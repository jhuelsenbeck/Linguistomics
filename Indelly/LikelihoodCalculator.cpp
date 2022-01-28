#include <cmath>
#include <iomanip>
#include <iostream>
#include "EigenSystem.hpp"
#include "IntVector.hpp"
#include "LikelihoodCalculator.hpp"
#include "Model.hpp"
#include "Msg.hpp"
#include "Node.hpp"
#include "ParameterAlignment.hpp"
#include "TransitionProbabilities.hpp"
#include "Tree.hpp"

#undef DEBUG_LIKELIHOOD

int LikelihoodCalculator::maxUnalignableDimension = 10;



LikelihoodCalculator::LikelihoodCalculator(ParameterAlignment* a, Model* m) {

    // memory addresses of important objects
    data = a;
    model = m;
    TransitionProbabilities& tp = TransitionProbabilities::transitionProbabilties();
    transitionProbabilityFactory = &tp;
    sequences = data->getRawSequenceMatrix();

    numStates = data->getNumStates();
    taxonMask = data->getTaxonMask();
    numTaxa = (int)taxonMask.getNumberSetBits();
    numNodes = 2 * numTaxa - 1;
    
    beta.resize(numNodes, 0.0);
    birthProbability.resize(numNodes, 0.0);
    extinctionProbability.resize(numNodes, 0.0);
    homologousProbability.resize(numNodes, 0.0);
    nonHomologousProbability.resize(numNodes, 0.0);

    fH = new double*[numNodes];
    fH[0] = new double[numNodes * numStates];
    for (int i=1; i<numNodes; i++)
        fH[i] = fH[i-1] + numStates;
    
    fI = new double*[numNodes];
    fI[0] = new double[numNodes * (numStates + 1)];
    for (int i=1; i<numNodes; i++)
        fI[i] = fI[i-1] + (numStates + 1);

    zeroH = new double[numStates];
    for (int i=0; i<numStates; i++)
        zeroH[i] = 0.0;
        
    zeroI = new double[numStates + 1];
    for (int i=0; i<numStates+1; i++)
        zeroI[i] = 0.0;
        
    possibleVectorIndices.resize(maxUnalignableDimension);
    nodeHomology.resize(numNodes);
    numHomologousEmissions.resize(numNodes);
    numHomologousEmissionsForClass.resize(maxUnalignableDimension + 1);
}

LikelihoodCalculator::~LikelihoodCalculator(void) {

    drainPool();
    delete [] fH[0];
    delete [] fH;
    delete [] fI[0];
    delete [] fI;
    delete [] zeroH;
    delete [] zeroI;
}

std::string LikelihoodCalculator::alignmentName(void) {

    return data->getName();
}

void LikelihoodCalculator::clearPpTable(void) {

    for (PartialProbabilitiesLookup::iterator it = partialProbabilities.begin(); it != partialProbabilities.end(); it++)
        returnToPool(it->first);
    partialProbabilities.clear();
}


void LikelihoodCalculator::drainPool(void) {

    for (std::vector<IntVector*>::iterator v=pool.begin(); v != pool.end(); v++)
        {
        allocated.erase(*v);
        delete (*v);
        }
}

IntVector* LikelihoodCalculator::getVector(void) {

    if (pool.empty() == true)
        {
        /* If the vector pool is empty, we allocate a new vector and return it. We
           do not need to add it to the vector pool. */
        IntVector* v = new IntVector(numTaxa);
        allocated.insert(v);
#       if defined (TRACK_ALLOCS)
        numAllocs++;
#       endif
        return v;
        }
    
    // Return a vector from the vector pool, remembering to remove it from the pool.
    IntVector* v = pool.back();
    pool.pop_back();
    return v;
}

IntVector* LikelihoodCalculator::getVector(IntVector& vec) {

    if (pool.empty() == true)
        {
        IntVector* v = new IntVector(numTaxa);
        allocated.insert(v);
        for (int i=0; i<vec.size(); ++i)
            (*v)[i] = vec[i];
#       if defined (TRACK_ALLOCS)
        numAllocs++;
#       endif
        return v;
        }
    IntVector* v = pool.back();
    pool.pop_back();
    for (int i=0; i<vec.size(); ++i)
        (*v)[i] = vec[i];
    return v;
}

void LikelihoodCalculator::initialize(void) {

    alignment = data->getIndelMatrix();
    unalignableRegionSize = 0;
    tree = model->getTree(taxonMask);

    RbBitSet bs = data->getTaxonMask();
    transitionProbabilities = transitionProbabilityFactory->getTransitionProbabilities(bs);
    std::vector<Node*>& downPassSequence = tree->getDownPassSequence();
    for (int n=0; n<downPassSequence.size(); n++)
        {
        Node* p = downPassSequence[n];
        p->setTransitionProbability( transitionProbabilities[p->getIndex()] );
        }
    stateEquilibriumFrequencies = transitionProbabilityFactory->getStationaryFrequencies();
    setBirthDeathProbabilities();
}

double LikelihoodCalculator::lnLikelihood(void) {

    initialize();
    
    // null emissions probability
	IntVector* pos = getVector();
    double nullEmissionFactor = partialProbability(pos, pos);
    double f = immortalProbability / nullEmissionFactor;
    partialProbabilities.insert( std::make_pair(pos, f) );

    // array of possible vector indices, used in inner loop
    //std::vector<int> possibleVectorIndices(maxUnalignableDimension, 0);
    for (int i=0; i<maxUnalignableDimension; i++)
        possibleVectorIndices[i] = 0;
    int firstNotUsed = 0;
    int len = (int)alignment.size();
    int currentAlignmentColumn = 0;
    //std::vector<int> state(len);
    state.resize(len);
    do
        {
        // Find all possible vectors from current position, pos
        int ptr;
        int numPossible = 0;
        IntVector* mask = getVector();
	    for (ptr=firstNotUsed; mask->zeroEntry() && ptr<len; ptr++)
            {
            if (state[ ptr ] != used)
                {
                if (mask->innerProduct( alignment[ptr] ) == 0)
                    {
                    state[ ptr ] = possible;
                    if (numPossible == maxUnalignableDimension)
                        {
                        // This gets too hairy - bail out.
                        unalignableRegionSize++;
                        Msg::warning("We bailed out because things became too hairy: numPossible = " + std::to_string(numPossible));
                        return -std::numeric_limits<double>::infinity();
                        }
                    possibleVectorIndices[numPossible++] = ptr;
                    }
                mask->add( alignment[ptr] );
                }
            }
        returnToPool(mask);

        // Loop over all combinations of possible vectors, which define edges from
        // pos to another possible position, by ordinary binary counting.
	    //IntVector newPos(*pos);
	    //IntVector signature(pos->size());
        IntVector* newPos = getVector(*pos);
        IntVector* signature = getVector();
	    bool unusedPos, foundNonZero;
        do
            {
            // Find next combination
            foundNonZero = false;
			for (int posPtr=numPossible-1; posPtr>=0; --posPtr)
                {
			    int curPtr = possibleVectorIndices[ posPtr ];
                currentAlignmentColumn = curPtr;
                if (state[ curPtr ] == possible)
                    {
                    state[ curPtr ] = edgeUsed;
					newPos->add( alignment[ curPtr ] );
                    // Compute signature vector
					signature->addMultiple( alignment[ curPtr ], posPtr + 1 );
                    // Signal: non-zero combination found, and stop
                    foundNonZero = true;
                    posPtr = 0;
                    }
                else
                    {
                    // It was edgeUsed (i.e., digit == 1), so reset digit and continue
					state[ curPtr ] = possible;
					newPos->subtract( alignment[ curPtr ] );
					signature->addMultiple( alignment[ curPtr ], -posPtr - 1 );
                    }
                }

            if (foundNonZero)
                {
                PartialProbabilitiesLookup::iterator it = partialProbabilities.find(pos);
                if (it == partialProbabilities.end())
                    Msg::error("Could not find pos vector in dpTable map.");
                double lft = it->second;
                it = partialProbabilities.find(newPos);
                double rht = 0.0;
                if (it == partialProbabilities.end())
                    {
                    unusedPos = true;
                    rht = 0.0;
                    }
                else
                    {
                    rht = it->second;
                    unusedPos = false;
                    }
                    
                double transFac = (-partialProbability(signature, pos)) / nullEmissionFactor;
                lft *= transFac;
                rht += lft;

                // If we are storing a value at a previously unused position, make sure we use a fresh key object
                if (unusedPos)
                    {
                    IntVector* newPosCopy = getVector(*newPos);
                    partialProbabilities.insert( std::make_pair(newPosCopy, rht) );
                    }
                else
                    {
                    PartialProbabilitiesLookup::iterator it2 = partialProbabilities.find(newPos);
                    if (it2 == partialProbabilities.end())
                        Msg::error("Should have found newPos in table. What happened?");
                    it2->second = rht;
                    //delete rht; // need to remember to delete rht if it's not being inserted back into the table
                    }
                }

            } while (foundNonZero);

    returnToPool(signature);

    // Now find next entry in partial probability table.  Use farthest unused vector
    --ptr;
    while (ptr >= 0 && state[ptr] != possible)
        {
        // Undo any possible used vector that we encounter
        if (state[ptr] == used)
            {
            pos->subtract( alignment[ptr] );
            state[ptr] = free;
            }
        --ptr;
        }
    
    if (ptr == -1)
        {
        // No more unused vectors, so we also fell through the edge loop above,
        // hence newPos contains the final position
        PartialProbabilitiesLookup::iterator it = partialProbabilities.find(newPos);
        if (it == partialProbabilities.end())
            {
            printTable();
            Msg::error("Could not find newPos in partialProbabilities");
            }
        double lnL = log(it->second);
        if (isinf(lnL) == true)
            {
            std::cout << "res = " << it->second << std::endl;
            data->print();
            printTable();
            }
        clearPpTable();
        returnToPool(newPos);
        if (pool.size() - allocated.size() != 0)
            Msg::warning("Some IntVectors were not returned to the pool");
        return lnL;
        }

    returnToPool(newPos);

    // Now use this farthest-out possible vector
    state[ptr] = used;
    pos->add( alignment[ptr] );
    if (ptr <= firstNotUsed)
        firstNotUsed++;
            
    } while (true);
}

#if 0

double LikelihoodCalculator::partialProbability(IntVector* signature, IntVector* pos) {

    // make certain that the conditional probabilities for all
    // nodes are zero before we start
    for (int i=0; i<numNodes; i++)
        memcpy( fH[i], zeroH, numStates*sizeof(double) );
    for (int i=0; i<numNodes; i++)
        memcpy( fI[i], zeroI, (numStates+1)*sizeof(double) );
        
    std::vector<int> nodeHomology(numNodes);
    std::vector<int> numHomologousEmissions(numNodes);
    std::vector<int> numHomologousEmissionsForClass(maxUnalignableDimension + 1);
    for (int i=0; i<numTaxa; i++)
        {
        numHomologousEmissionsForClass[ (*signature)[ i ] ]++;
        nodeHomology[ i ] = (*signature)[ i ];
        if ( (*signature)[i] == 0)
            numHomologousEmissions[ i ] = 0;
        else
            numHomologousEmissions[ i ] = 1;
        }
        
    // Loop over all nodes except root, and find out which nodes need carry nucleotides of what
    // homology class.
    bool isHomologyClashing = false;
    std::vector<Node*>& dpSequence = tree->getDownPassSequence();
    for (int n=0; n<dpSequence.size(); n++)
        {
        Node* p = dpSequence[n];
        if (p->getAncestor() != NULL)
            {
            int pIdx = p->getIndex();
            int pAncIdx = p->getAncestor()->getIndex();
            
            if ((numHomologousEmissions[pIdx] == numHomologousEmissionsForClass[ nodeHomology[pIdx] ]) || (nodeHomology[pIdx] == 0))
                {
                // This node is the MRCA of the homologous nucleotides iHom[i], or a gap --- do nothing
                }
            else
                {
                if (nodeHomology[pAncIdx] == 0 || nodeHomology[pAncIdx] == nodeHomology[pIdx])
                    {
                    // This node is not yet MRCA, so node above must be homologous
                    nodeHomology[pAncIdx] = nodeHomology[pIdx];
                    // If iHom[i]==0, iHomNum[i]==0.
                    numHomologousEmissions[pAncIdx] += numHomologousEmissions[pIdx];
                    }
                else
                    {
                    // Clashing homology; signal
                    isHomologyClashing = true;
                    }
                }
            }
        }
    if (isHomologyClashing)
        return 0.0;

    // pruning algorithm
    for (int n=0; n<dpSequence.size(); n++)
        {
        Node* p = dpSequence[n];
        int pIdx = p->getIndex();
        double* fI_p = &fI[pIdx][0];
        double* fH_p = &fH[pIdx][0];
        
        if (p->getIsLeaf() == true)
            {
            // initialize conditional probabilties for the tip
            if ( (*signature)[pIdx] == 0 )
                {
                // gap
                fI_p[numStates] = 1.0;
                }
            else
                {
                // nucleotide
                fH_p[ sequences[pIdx][(*pos)[pIdx]] ] = 1.0;
                }
            }
        else
            {
            // calculate conditional probabilities for interior nodes
            //std::vector<Node*> des = p->getDescendantsVector();
            p->getDescendantsVector(des);
            if (des.size() != 2)
                Msg::error("Expecting two descendants");
            Node* pLft = des[0];
            Node* pRht = des[1];
            int lftChildIdx = pLft->getIndex();
            int rhtChildIdx = pRht->getIndex();
            StateMatrix_t& tpLft = *(pLft->getTransitionProbability());
            StateMatrix_t& tpRht = *(pRht->getTransitionProbability());
            
            if ( (nodeHomology[pIdx] != 0) && (nodeHomology[pIdx] == nodeHomology[lftChildIdx]) && (nodeHomology[pIdx] == nodeHomology[rhtChildIdx]) )
                {
                // Case 1: One homology family spanning tree intersects both child edges, i.e.
                //         a 'homologous' nucleotide should travel down both edges.
                double* fH_lft = &fH[lftChildIdx][0];
                double* fH_rht = &fH[rhtChildIdx][0];
                double homologousProbability_lft = homologousProbability[lftChildIdx];
                double homologousProbability_rht = homologousProbability[rhtChildIdx];
                for (int i=0; i<numStates; ++i)
                    {
                    double sumLft = 0.0, sumRht = 0.0;
                    for (int j=0; j<numStates; ++j)
                        {
                        sumLft += fH_lft[j] * homologousProbability_lft * tpLft(i,j);
                        sumRht += fH_rht[j] * homologousProbability_rht * tpRht(i,j);
                        }
                    fH_p[i] = sumLft * sumRht;
                    }
                }
            else if (nodeHomology[pIdx] != 0)
                {
                // Case 2: One homology family Spanning tree intersects one of the child edges, i.e.
                //        'homologous' nucleotide must travel down a specific edge.
                if (nodeHomology[pIdx] == nodeHomology[lftChildIdx])
                    {
                    int homologousChildIdx   = lftChildIdx;
                    int inhomologousChildIdx = rhtChildIdx;
                    StateMatrix_t& tpHomologous   = *(pLft->getTransitionProbability());
                    StateMatrix_t& tpInhomologous = *(pRht->getTransitionProbability());
                    double* fH_homo   = &fH[homologousChildIdx][0];
                    double* fH_inhomo = &fH[inhomologousChildIdx][0];
                    double* fI_inhomo = &fI[inhomologousChildIdx][0];
                    double homologousProbability_homo   = homologousProbability[homologousChildIdx];
                    double homologousProbability_inhomo = homologousProbability[inhomologousChildIdx];
                    double nonHomologousProbability_inhomo = nonHomologousProbability[inhomologousChildIdx];
                    double extinctionProbability_inhomo = extinctionProbability[inhomologousChildIdx];
                    double birthProbability_inhomo = birthProbability[inhomologousChildIdx];
                    for (int i=0; i<numStates; i++)
                        {
                        double sumLft = 0.0;
                        double sumRht = extinctionProbability_inhomo * fI_inhomo[numStates];
                        for (int j=0; j<numStates; j++)
                            {
                            sumLft += fH_homo[j] * homologousProbability_homo * tpHomologous(i,j);
                            sumRht +=
                                (fH_inhomo[j] + fI_inhomo[j]) * (nonHomologousProbability_inhomo - extinctionProbability_inhomo * birthProbability_inhomo) * stateEquilibriumFrequencies[j] +
                                fI_inhomo[j] * homologousProbability_inhomo * tpInhomologous(i,j);
                            }
                        fH_p[i] = sumLft * sumRht;
                        }
                    }
                else
                    {
                    int homologousChildIdx = rhtChildIdx;
                    int inhomologousChildIdx = lftChildIdx;
                    StateMatrix_t& tpHomologous = *(pRht->getTransitionProbability());
                    StateMatrix_t& tpInhomologous = *(pLft->getTransitionProbability());
                    double* fH_homo   = &fH[homologousChildIdx][0];
                    double* fH_inhomo = &fH[inhomologousChildIdx][0];
                    double* fI_inhomo = &fI[inhomologousChildIdx][0];
                    double homologousProbability_homo   = homologousProbability[homologousChildIdx];
                    double homologousProbability_inhomo = homologousProbability[inhomologousChildIdx];
                    double nonHomologousProbability_inhomo = nonHomologousProbability[inhomologousChildIdx];
                    double extinctionProbability_inhomo = extinctionProbability[inhomologousChildIdx];
                    double birthProbability_inhomo = birthProbability[inhomologousChildIdx];
                    for (int i=0; i<numStates; i++)
                        {
                        double sumLft = 0.0;
                        double sumRht = extinctionProbability_inhomo * fI_inhomo[numStates];
                        for (int j=0; j<numStates; j++)
                            {
                            sumLft += fH_homo[j] * homologousProbability_homo * tpHomologous(i,j);
                            sumRht +=
                                (fH_inhomo[j] + fI_inhomo[j]) * (nonHomologousProbability_inhomo - extinctionProbability_inhomo * birthProbability_inhomo) * stateEquilibriumFrequencies[j] +
                                fI_inhomo[j] * homologousProbability_inhomo * tpInhomologous(i,j);
                            }
                        fH_p[i] = sumLft * sumRht;
                        }
                    }

                }
            else
                {
                // Case 3: No homology family's spanning tree intersect either child edge, i.e.
                //         a 'homologous' nucleotide may travel down either edge and pop out at a leaf,
                //         or may travel down both edges but has to die in at least one subtree, or
                //         an 'inhomologous' nucleotide may do whatever it likes here.
                double* fH_lft                      = &fH[lftChildIdx][0];
                double* fH_rht                      = &fH[rhtChildIdx][0];
                double* fI_lft                      = &fI[lftChildIdx][0];
                double* fI_rht                      = &fI[rhtChildIdx][0];
                double birthProbability_lft         = birthProbability[lftChildIdx];
                double birthProbability_rht         = birthProbability[rhtChildIdx];
                double homologousProbability_lft    = homologousProbability[lftChildIdx];
                double homologousProbability_rht    = homologousProbability[rhtChildIdx];
                double nonHomologousProbability_lft = nonHomologousProbability[lftChildIdx];
                double nonHomologousProbability_rht = nonHomologousProbability[rhtChildIdx];
                double extinctionProbability_lft    = extinctionProbability[lftChildIdx];
                double extinctionProbability_rht    = extinctionProbability[rhtChildIdx];
                for (int i=0; i<numStates; ++i)
                    {
                    double sumLft1 = 0.0, sumLft2 = 0.0;
                    double sumRht1 = extinctionProbability_lft * fI[lftChildIdx][numStates];
                    double sumRht2 = extinctionProbability_rht * fI[rhtChildIdx][numStates];
                    for (int j=0; j<numStates; ++j)
                        {
                        sumLft1 += fH_lft[j] * homologousProbability_lft * tpLft(i,j);
                        sumLft2 += fH_rht[j] * homologousProbability_rht * tpRht(i,j);
                        sumRht1 += (fH_lft[j] + fI_lft[j]) * (nonHomologousProbability_lft - extinctionProbability_lft * birthProbability_lft) * stateEquilibriumFrequencies[j] + fI_lft[j] * homologousProbability_lft * tpLft(i,j);
                        sumRht2 += (fH_rht[j] + fI_rht[j]) * (nonHomologousProbability_rht - extinctionProbability_rht * birthProbability_rht) * stateEquilibriumFrequencies[j] + fI_rht[j] * homologousProbability_rht * tpRht(i,j);
                        }
                    fH_p[i] = sumLft1 * sumRht2 + sumLft2 * sumRht1;
                    fI_p[i] = sumRht1 * sumRht2;
                    }
                double sumLft = fI_lft[numStates];
                double sumRht = fI_rht[numStates];
                for (int j=0; j<numStates; ++j)
                    {
                    sumLft -= birthProbability_lft * (fI_lft[j] + fH_lft[j]) * stateEquilibriumFrequencies[j];
                    sumRht -= birthProbability_rht * (fI_rht[j] + fH_rht[j]) * stateEquilibriumFrequencies[j];
                    }
                fI_p[numStates] = sumLft * sumRht;
               }
            
            }
        }

    // average probability over states at the root
    int rootIdx = tree->getRoot()->getIndex();
    double* fI_root = &fI[rootIdx][0];
    double* fH_root = &fH[rootIdx][0];
    double birthProbability_root = birthProbability[rootIdx];
    double res = fI_root[numStates];
    for (int i=0; i<numStates; i++)
        res -= (fI_root[i] + fH_root[i]) * birthProbability_root * stateEquilibriumFrequencies[i];
    return res;
}

#else

double LikelihoodCalculator::partialProbability(IntVector* signature, IntVector* pos) {

    // make certain that the conditional probabilities for all
    // nodes are zero before we start
    for (int i=0; i<numNodes; i++)
        memcpy( fH[i], zeroH, numStates*sizeof(double) );
    for (int i=0; i<numNodes; i++)
        memcpy( fI[i], zeroI, (numStates+1)*sizeof(double) );
        
//    std::vector<int> nodeHomology(numNodes);
//    std::vector<int> numHomologousEmissions(numNodes);
//    std::vector<int> numHomologousEmissionsForClass(maxUnalignableDimension + 1);
    for (int i=0; i<numNodes; i++)
        {
        nodeHomology[i] = 0;
        numHomologousEmissions[i] = 0;
        }
    for (int i=0; i<maxUnalignableDimension + 1; i++)
        numHomologousEmissionsForClass[i] = 0;
    for (int i=0; i<numTaxa; i++)
        {
        numHomologousEmissionsForClass[ (*signature)[ i ] ]++;
        nodeHomology[ i ] = (*signature)[ i ];
        if ( (*signature)[i] == 0)
            numHomologousEmissions[ i ] = 0;
        else
            numHomologousEmissions[ i ] = 1;
        }
        
    // Loop over all nodes except root, and find out which nodes need carry nucleotides of what
    // homology class.
    bool isHomologyClashing = false;
    std::vector<Node*>& dpSequence = tree->getDownPassSequence();
    for (int n=0; n<dpSequence.size(); n++)
        {
        Node* p = dpSequence[n];
        if (p->getAncestor() != NULL)
            {
            int pIdx = p->getIndex();
            int pAncIdx = p->getAncestor()->getIndex();
            
            if ((numHomologousEmissions[pIdx] == numHomologousEmissionsForClass[ nodeHomology[pIdx] ]) || (nodeHomology[pIdx] == 0))
                {
                // This node is the MRCA of the homologous nucleotides iHom[i], or a gap --- do nothing
                }
            else
                {
                if (nodeHomology[pAncIdx] == 0 || nodeHomology[pAncIdx] == nodeHomology[pIdx])
                    {
                    // This node is not yet MRCA, so node above must be homologous
                    nodeHomology[pAncIdx] = nodeHomology[pIdx];
                    // If iHom[i]==0, iHomNum[i]==0.
                    numHomologousEmissions[pAncIdx] += numHomologousEmissions[pIdx];
                    }
                else
                    {
                    // Clashing homology; signal
                    isHomologyClashing = true;
                    }
                }
            }
        }
    if (isHomologyClashing)
        return 0.0;
 

    // pruning algorithm
    for (int n=0; n<dpSequence.size(); n++)
        {
        Node* p = dpSequence[n];
        int pIdx = p->getIndex();
        if (p->getIsLeaf() == true)
            {
            // initialize conditional probabilties for the tip
            if ( (*signature)[pIdx] == 0 )
                {
                // gap
                fI[pIdx][numStates] = 1.0;
                }
            else
                {
                // nucleotide
                fH[pIdx][ sequences[pIdx][(*pos)[pIdx]] ] = 1.0;
                }
            }
        else
            {
            // calculate conditional probabilities for interior nodes
            p->getDescendantsVector(des);
            //std::vector<Node*> des = p->getDescendantsVector();
            if (des.size() != 2)
                Msg::error("Expecting two descendants");
            Node* pLft = des[0];
            Node* pRht = des[1];
            int lftChildIdx = pLft->getIndex();
            int rhtChildIdx = pRht->getIndex();
            StateMatrix_t& tpLft = *(pLft->getTransitionProbability());
            StateMatrix_t& tpRht = *(pRht->getTransitionProbability());
            
            if ( (nodeHomology[pIdx] != 0) && (nodeHomology[pIdx] == nodeHomology[lftChildIdx]) && (nodeHomology[pIdx] == nodeHomology[rhtChildIdx]) )
                {
                // Case 1: One homology family spanning tree intersects both child edges, i.e.
                //         a 'homologous' nucleotide should travel down both edges.
                for (int i=0; i<numStates; i++)
                    {
                    double lft = 0.0;
                    double rht = 0.0;
                    for (int j=0; j<numStates; j++)
                        {
                        lft += fH[lftChildIdx][j] * homologousProbability[lftChildIdx] * tpLft(i,j);
                        rht += fH[rhtChildIdx][j] * homologousProbability[rhtChildIdx] * tpRht(i,j);
                        }
                    fH[pIdx][i] = lft * rht;
                    }
                }
            else if (nodeHomology[pIdx] != 0)
                {
                // Case 2: One homology family Spanning tree intersects one of the child edges, i.e.
                //        'homologous' nucleotide must travel down a specific edge.
                int homologousChildIdx, inhomologousChildIdx;
                StateMatrix_t* tpHomologous = NULL;
                StateMatrix_t* tpInhomologous = NULL;
                if (nodeHomology[pIdx] == nodeHomology[lftChildIdx])
                    {
                    homologousChildIdx = lftChildIdx;
                    inhomologousChildIdx = rhtChildIdx;
                    tpHomologous = pLft->getTransitionProbability();
                    tpInhomologous = pRht->getTransitionProbability();
                    }
                else
                    {
                    homologousChildIdx = rhtChildIdx;
                    inhomologousChildIdx = lftChildIdx;
                    tpHomologous = pRht->getTransitionProbability();
                    tpInhomologous = pLft->getTransitionProbability();
                    }

                for (int i=0; i<numStates; i++)
                    {
                    double lft = 0.0;
                    double rht = extinctionProbability[inhomologousChildIdx] * fI[inhomologousChildIdx][numStates];
                    for (int j=0; j<numStates; j++)
                        {
                        lft += fH[homologousChildIdx][j] * homologousProbability[homologousChildIdx] * (*tpHomologous)(i,j);
                        rht +=
                            (fH[inhomologousChildIdx][j] + fI[inhomologousChildIdx][j]) * (nonHomologousProbability[inhomologousChildIdx] - extinctionProbability[inhomologousChildIdx] * birthProbability[inhomologousChildIdx]) * stateEquilibriumFrequencies[j] +
                            fI[inhomologousChildIdx][j] * homologousProbability[inhomologousChildIdx] * (*tpInhomologous)(i,j);
                        }
                    fH[pIdx][i] = lft * rht;
                    }
                }
            else
                {
                // Case 3: No homology family's spanning tree intersect either child edge, i.e.
                //         a 'homologous' nucleotide may travel down either edge and pop out at a leaf,
                //         or may travel down both edges but has to die in at least one subtree, or
                //         an 'inhomologous' nucleotide may do whatever it likes here.
                for (int i=0; i<numStates; i++)
                    {
                    double lft1 = 0.0, lft2 = 0.0;
                    double rht1 = extinctionProbability[lftChildIdx] * fI[lftChildIdx][numStates];
                    double rht2 = extinctionProbability[rhtChildIdx] * fI[rhtChildIdx][numStates];
                    for (int j=0; j<numStates; j++)
                        {
                        lft1 += fH[lftChildIdx][j] * homologousProbability[lftChildIdx] * tpLft(i,j);
                        lft2 += fH[rhtChildIdx][j] * homologousProbability[rhtChildIdx] * tpRht(i,j);
                        rht1 += (fH[lftChildIdx][j] + fI[lftChildIdx][j]) * (nonHomologousProbability[lftChildIdx] - extinctionProbability[lftChildIdx] * birthProbability[lftChildIdx]) * stateEquilibriumFrequencies[j] + fI[lftChildIdx][j] * homologousProbability[lftChildIdx] * tpLft(i,j);
                        rht2 += (fH[rhtChildIdx][j] + fI[rhtChildIdx][j]) * (nonHomologousProbability[rhtChildIdx] - extinctionProbability[rhtChildIdx] * birthProbability[rhtChildIdx]) * stateEquilibriumFrequencies[j] + fI[rhtChildIdx][j] * homologousProbability[rhtChildIdx] * tpRht(i,j);
                        }
                    fH[pIdx][i] = lft1 * rht2 + lft2 * rht1;
                    fI[pIdx][i] = rht1 * rht2;
                    }
                double lft = fI[lftChildIdx][numStates];
                double rht = fI[rhtChildIdx][numStates];
                for (int j=0; j<numStates; j++)
                    {
                    lft -= birthProbability[lftChildIdx] * (fI[lftChildIdx][j] + fH[lftChildIdx][j]) * stateEquilibriumFrequencies[j];
                    rht -= birthProbability[rhtChildIdx] * (fI[rhtChildIdx][j] + fH[rhtChildIdx][j]) * stateEquilibriumFrequencies[j];
                    }
                fI[pIdx][numStates] = lft * rht;
               }
            
            }
        }

    // average probability over states at the root
    int rootIdx = tree->getRoot()->getIndex();
    double res = fI[rootIdx][numStates];
    for (int i=0; i<numStates; i++)
        res -= (fI[rootIdx][i] + fH[rootIdx][i]) * birthProbability[rootIdx] * stateEquilibriumFrequencies[i];
    return res;
}
#endif

void LikelihoodCalculator::printTable(void) {

    std::cout << "Table" << std::endl;
    std::cout << std::fixed << std::setprecision(10) << std::scientific;
    int i = 0;
    for (PartialProbabilitiesLookup::iterator it = partialProbabilities.begin(); it != partialProbabilities.end(); it++)
        {
        std::cout << i++ << " -- " << *(it->first) << " -- " << it->second << std::endl;
        }
}

void LikelihoodCalculator::returnToPool(IntVector* v) {

    v->clean();
    pool.push_back(v);
}

void LikelihoodCalculator::setBirthDeathProbabilities(void) {

    insertionRate = model->getInsertionRate();
    deletionRate  = model->getDeletionRate();

    immortalProbability= 1.0;
    std::vector<Node*>& dpSequence = tree->getDownPassSequence();
    for (int n=0; n<dpSequence.size(); n++)
        {
        Node* p = dpSequence[n];
        int pIdx = p->getIndex();
        double v = p->getBranchLength();
        if (p->getAncestor() == NULL)
            {
            // root
            beta[pIdx] = 1.0 / deletionRate;
            homologousProbability[pIdx] = 0.0;
            }
        else
            {
            // internal node or tip
            beta[pIdx] = exp( (insertionRate - deletionRate) * v );
            beta[pIdx] = (1.0 - beta[pIdx]) / (deletionRate - insertionRate * beta[pIdx]);
            homologousProbability[pIdx] = exp( -deletionRate * v) * (1.0 - insertionRate * beta[pIdx] );
            }
        birthProbability[pIdx]         = insertionRate * beta[pIdx];
        extinctionProbability[pIdx]    = deletionRate * beta[pIdx];
        nonHomologousProbability[pIdx] = (1.0 - deletionRate * beta[pIdx]) * (1.0 - birthProbability[pIdx]) - homologousProbability[pIdx];
        immortalProbability *= (1.0 - birthProbability[pIdx]);
        }
        
#   if defined(DEBUG_LIKELIHOOD)
    std::cout << "TKF91 event probabilties" << std::endl;
    std::cout << std::fixed << std::setprecision(10);
    for (int i=0; i<numNodes; i++)
        {
        std::cout << i << " -- " << birthProbability[i] << " " << extinctionProbability[i] << " " << homologousProbability[i] << " " << nonHomologousProbability[i] << std::endl;
        }
#   endif
}
