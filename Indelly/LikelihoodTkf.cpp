#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include "EigenSystem.hpp"
#include "IntVector.hpp"
#include "Msg.hpp"
#include "LikelihoodTkf.hpp"
#include "Model.hpp"
#include "Node.hpp"
#include "ParameterAlignment.hpp"
#include "TransitionProbabilities.hpp"
#include "Tree.hpp"
#include "UserSettings.hpp"

#include <gmpxx.h>

#undef DEBUG_TKF91
int LikelihoodTkf::unalignableRegionSize = 0;
int LikelihoodTkf::maxUnalignableDimension = 10;

/* This is code from Beast (now defunct in Beast) for calculating the likelihood
   under the TKF91 model. The main changes were made in the hash map called iTable
   in the original code and called dpTable here. I use a map, sorting on the vector
   contents for the key. The original code hashes the vectors by the vector elements,
   so both should be equivalent, here. I also use a high performance arbitrary
   precision double in place of the roll-it-yourself class (BFloat) used in Beast. */



LikelihoodTkf::LikelihoodTkf(ParameterAlignment* a, Tree* t, Model* m) {

//    mpf_class x;

    data = a;
    tree = t;
    model = m;
    init();
}

LikelihoodTkf::~LikelihoodTkf(void) {

    clearDpTable();
}

void LikelihoodTkf::clearDpTable(void) {

//    for (std::map<IntVector*,mpfr::mpreal,CompIntVector>::iterator it = dpTable.begin(); it != dpTable.end(); it++)
//        delete it->first;
    for (std::map<IntVector*,double,CompIntVector>::iterator it = dpTable.begin(); it != dpTable.end(); it++)
        delete it->first;
    dpTable.clear();
}

void LikelihoodTkf::init(void) {

    // initialize some useful variables
    numStates = data->getNumStates();

    // initialize the iParent and iTau arrays based on the given tree.
    initTree(tree);
      
    // initialize the iAlignment array from the given alignment.
    initAlignment();

    // initialize the iSequences array from the given alignment.
    initSequences();

    // initialize the iTrans array from the substitution model -- must be called after populating tree!
    initTransitionProbabilities();

    // Initialise TKF91 coefficients in iB, iH, iN, iE, and iInitial
    initTKF91();
}

void LikelihoodTkf::initAlignment(void) {
            
    alignment = data->getIndelMatrix();
            
//    printVector("Alignment", alignment);
#   if 0
    std::cout << "Variable = alignment" << std::endl;
    std::cout << "name = " << data->getName() << std::endl;
    for (int i=0; i<alignment.size(); i++)
        {
        std::cout << std::setw(3) << i << " -- ";
        for (int j=0; j<alignment[i].size(); j++)
            std::cout << alignment[i][j] << " ";
        std::cout << std::endl;
        }
#   endif
}

void LikelihoodTkf::initSequences(void) {
    
    sequences = data->getRawSequenceMatrix();

#   if 0
    std::cout << "Variable = sequences" << std::endl;
    for (int i=0; i<sequences.size(); i++)
        {
        std::cout << std::setw(3) << i << " -- ";
        for (int j=0; j<sequences[i].size(); j++)
            std::cout << std::setw(2) << sequences[i][j] << " ";
        std::cout << std::endl;
        }
#   endif
}

void LikelihoodTkf::initTransitionProbabilities(void) {
        
    TransitionProbabilities& tip = TransitionProbabilities::transitionProbabilties();
    transitionProbabilities = tip.getTransitionProbabilities();
    stateEquilibriumFrequencies = tip.getStationaryFrequencies();

#   if 0
    tip.print();
#   endif
}

void LikelihoodTkf::initTKF91(void) {

    // initialize lambda and mu
    insertionRate = model->getInsertionRate();
    deletionRate  = model->getDeletionRate();

    // initialize parameters of TKF91 model
    UserSettings& settings = UserSettings::userSettings();
    numIndelCategories = settings.getNumIndelCategories();

    int numNodes = (int)parents.size();
    std::vector<double> beta(numNodes);
    birthProbability.resize(numNodes);
    extinctionProbability.resize(numNodes);
    homologousProbability.resize(numNodes);
    nonHomologousProbability.resize(numNodes);

    if (numIndelCategories == 1)
        {
        
        immortalProbability = 1.0;
        for (int i=0; i<numNodes; i++)
            {
            if (i == numNodes-1)
                {
                // root
                beta[ i ] = 1.0 / deletionRate;
                homologousProbability[i] = 0.0;
                }
            else
                {
                // internal node or tip
                beta[i] = exp( (insertionRate - deletionRate) * tau[i] );
                beta[i] = (1.0 - beta[i]) / (deletionRate - insertionRate * beta[i]);
                homologousProbability[i] = exp( -deletionRate * tau[i]) * (1.0 - insertionRate * beta[i] );
                }
            birthProbability[i]         = insertionRate * beta[i];
            extinctionProbability[i]    = deletionRate * beta[i];
            nonHomologousProbability[i] = (1.0 - deletionRate * beta[i]) * (1.0 - birthProbability[i]) - homologousProbability[i];
            immortalProbability *= (1.0 - birthProbability[i]);
            }
        }
    else
        {
        
        }
        
#   if defined(DEBUG_TKF91)
    std::cout << "TKF91 event probabilties" << std::endl;
    std::cout << std::fixed << std::setprecision(10);
    for (int i=0; i<numNodes; i++)
        {
        std::cout << i << " -- " << birthProbability[i] << " " << extinctionProbability[i] << " " << homologousProbability[i] << " " << nonHomologousProbability[i] << std::endl;
        }
#   endif
}

void LikelihoodTkf::initTree(Tree* t) {

    parents = t->getAncestorIndices();
    tau = t->getBranchLengthVector();
    
    //printTreeInfo();

#   if defined(DEBUG_TKF91)
    std::cout << "Tree" << std::endl;
    for (int i=0; i<parents.size(); i++)
        {
        std::cout << std::setw(3) << i << " -- " << std::setw(4) << parents[i] << " ";
        std::cout << std::fixed << std::setprecision(4) << tau[i] << " ";
        std::cout << std::endl;
        }
#   endif
}

void LikelihoodTkf::printTable(void) {

    std::cout << "dpTable" << std::endl;
    std::cout << std::fixed << std::setprecision(10) << std::scientific;
    int i = 0;
//    for (std::map<IntVector*,mpfr::mpreal,CompIntVector>::iterator it = dpTable.begin(); it != dpTable.end(); it++)
//        {
//        std::cout << i++ << " -- " << *(it->first) << " -- " << it->second.toString() << std::endl;
//        }
    for (std::map<IntVector*,double,CompIntVector>::iterator it = dpTable.begin(); it != dpTable.end(); it++)
        {
        std::cout << i++ << " -- " << *(it->first) << " -- " << it->second << std::endl;
        }
}

void LikelihoodTkf::printTreeInfo(void) {

    for (int i=0; i<parents.size(); i++)
        {
        std::cout << "iParent[" << i << "] = " << parents[i] << ";" << std::endl;
        }
        
    for (int i=0; i<tau.size(); i++)
        {
        std::cout << std::fixed << std::setprecision(10);
        std::cout << "iTau[" << i << "] = " << tau[i] << ";" << std::endl;
        }
}

void LikelihoodTkf::printVector(std::string header, std::vector<int>& v) {

    if (header != "")
        std::cout << header << std::endl;
    for (int i=0; i<v.size(); i++)
        std::cout << v[i] << " ";
    std::cout << std::endl;
}

void LikelihoodTkf::printVector(std::string header, std::vector<double>& v) {

    std::cout << std::fixed << std::setprecision(10);
    
    if (header != "")
        std::cout << header << std::endl;
    for (int i=0; i<v.size(); i++)
        std::cout << v[i] << " ";
    std::cout << std::endl;
}

void LikelihoodTkf::printVector(std::string header, std::vector< std::vector<double> >& v) {

    std::cout << std::fixed << std::setprecision(10);

    if (header != "")
        std::cout << header << std::endl;
    for (int i=0; i<v.size(); i++)
        {
        for (int j=0; j<v[i].size(); j++)
            std::cout << v[i][j] << " ";
        std::cout << std::endl;
        }
}

void LikelihoodTkf::printVector(std::string header, std::vector< std::vector<int> >& v) {

    if (header != "")
        std::cout << header << std::endl;
    for (int i=0; i<v.size(); i++)
        {
        for (int j=0; j<v[i].size(); j++)
            std::cout << std::setw(2) << v[i][j] << " ";
        std::cout << std::endl;
        }
}

double LikelihoodTkf::tkfLike(void) {
    
    int len = (int)alignment.size();
    int numLeaves = (int)alignment[0].size();
    int firstNotUsed = 0;                                    // First not-'used' alignment vector (for efficiency)
    std::vector<int> state(len);                             // Helper array, to traverse the region in the DP table corresp. to the alignment
	IntVector pos(numLeaves);                                // Current position; sum of all used vectors
    int currentAlignmentColumn = 0;                          // Keep track of which alignment column is being calculated
    
    // Calculate correction factor for null emissions ("wing folding", or linear equation solving.)
    double nullEmissionFactor = treeRecursion( &pos, &pos, -1 );
    double f = immortalProbability / nullEmissionFactor;
    dpTable.insert( std::make_pair(new IntVector(pos), f) );   // insert the first probability into dpTable
    //printTable();
 
    // Array of possible vector indices, used in inner loop
    std::vector<int> possibleVectorIndices(maxUnalignableDimension, 0);

    do
        {
        // Find all possible vectors from current position, pos
	    IntVector mask(numLeaves);
        int ptr;
        int numPossible = 0;
	    for (ptr=firstNotUsed; mask.zeroEntry() && ptr<len; ptr++)
            {
            if (state[ ptr ] != used)
                {
                if (mask.innerProduct( alignment[ptr] ) == 0)
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
                mask.add( alignment[ptr] );
                }
            }

        // Loop over all combinations of possible vectors, which define edges from
        // pos to another possible position, by ordinary binary counting.
	    IntVector newPos(pos);
	    IntVector signature(pos.size());
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
					newPos.add( alignment[ curPtr ] );
                    // Compute signature vector
					signature.addMultiple( alignment[ curPtr ], posPtr+1 );
                    // Signal: non-zero combination found, and stop
                    foundNonZero = true;
                    posPtr = 0;
                    }
                else
                    {
                    // It was edgeUsed (i.e., digit == 1), so reset digit and continue
					state[ curPtr ] = possible;
					newPos.subtract( alignment[ curPtr ] );
					signature.addMultiple( alignment[ curPtr ], -posPtr-1 );
                    }
                }

            if (foundNonZero)
                {
                std::map<IntVector*,double,CompIntVector>::iterator it = dpTable.find(&pos);
                if (it == dpTable.end())
                    Msg::error("Could not find pos vector in dpTable map.");
                double lft = it->second;
                it = dpTable.find(&newPos);
                double rht;
                if (it == dpTable.end())
                    {
                    unusedPos = true;
                    rht = 0.0;
                    }
                else
                    {
                    rht = it->second;
                    unusedPos = false;
                    }
                    
//                std::cout << "currentAlignmentColumn = " << currentAlignmentColumn << std::endl;
                double transFac = (-treeRecursion( &signature, &pos, currentAlignmentColumn )) / nullEmissionFactor;
                lft *= transFac;
                rht += lft;

                // If we are storing a value at a previously unused position, make sure we use a fresh key object
                if (unusedPos)
                    {
                    dpTable.insert( std::make_pair(new IntVector(newPos), rht) );
                    //printTable();
                    }
                else
                    {
                    std::map<IntVector*,double,CompIntVector>::iterator it2 = dpTable.find(&newPos);
                    if (it2 == dpTable.end())
                        Msg::error("Should have found newPos in table. What happened?");
                    it2->second = rht;
                    //printTable();
                    }
                }

            } while (foundNonZero);

    // Now find next entry in DP table.  Use farthest unused vector
    --ptr;
    while (ptr >= 0 && state[ptr] != possible)
        {
        // Undo any possible used vector that we encounter
        if (state[ptr] == used)
            {
            pos.subtract( alignment[ptr] );
            state[ptr] = free;
            }
        --ptr;
        }
    
    if (ptr == -1)
        {
        // No more unused vectors, so we also fell through the edge loop above,
        // hence newPos contains the final position
        std::map<IntVector*,double,CompIntVector>::iterator it = dpTable.find(&newPos);
        if (it == dpTable.end())
            Msg::error("Could not find newPos in dpTable");
        //std::cout << "found " << it->first << " " << it->second << " " << log(it->second) << std::endl;
        double lnL = log(it->second);
        clearDpTable();
        return lnL;
        }
    
    // Now use this farthest-out possible vector
    state[ptr] = used;
    pos.add( alignment[ptr] );
    if (ptr <= firstNotUsed)
        firstNotUsed++;
            
    } while (true);
}

void LikelihoodTkf::debugPrint(void) {

    std::cout << "Node information:" << std::endl;
    for (int i=0; i<parents.size(); i++)
        {
        std::cout << std::setw(3) << i << " : " << std::setw(3) << parents[i] << " " << std::fixed << std::setprecision(5) << tau[i] << " ";
        std::cout << birthProbability[i] << " " << extinctionProbability[i] << " " << homologousProbability[i] << " " << nonHomologousProbability[i];
        std::cout << std::endl;
        }
        
    std::cout << "Parameters:" << std::endl;
    std::cout << "  Lambda = " << insertionRate << std::endl;
    std::cout << "  Mu     = " << deletionRate << std::endl;
    std::cout << "  E(N)   = " << insertionRate / (deletionRate - insertionRate) << std::endl;
    
    for (int i=0; i<sequences.size(); i++)
        {
        for (int j=0; j<sequences[i].size(); j++)
            std::cout << sequences[i][j] << " ";
        std::cout << std::endl;
        }
        
    for (int i=0; i<transitionProbabilities.size(); i++)
        {
        std::cout << "iTrans[" << i << "] (" << tau[i] << "):" << std::endl;
        for (int j=0; j<numStates; j++)
            {
            std::cout << "  ";
            for (int k=0; k<numStates; k++)
                std::cout << std::fixed << std::setprecision(3) << transitionProbabilities[i][j][k] << " ";
            std::cout << std::endl;
            }
        }
}

double LikelihoodTkf::treeRecursion(IntVector* signature, IntVector* pos, int siteColumn) {

    int numLeaves = (int)signature->size();                                       // Dimension of alignment columns, i.e. number of leaves
    int numNodes = (int)parents.size();                                           // Number of internal nodes
    std::vector<int> nodeHomology(numNodes);                                      // Homology for every node; 0 if node need not be homologous to an emitted nucleotide
    std::vector<int> numHomologousEmissions(numNodes);                            // Number of homologous emissions accounted for by homologous nucleotide at this node
    std::vector<int> numHomologousEmissionsForClass(maxUnalignableDimension + 1); // Number of emissions for each class of homologous nucleotides
    std::vector<int> lftChildrenIndices(numNodes);                                // Left children
    std::vector<int> rhtChildrenIndices(numNodes);                                // Right children

    std::vector< std::vector<double> > fH(numNodes);
    std::vector< std::vector<double> > fI(numNodes);
	for (int i=0; i<numNodes; i++)
        {
	    fH[i].resize(numStates, 0.0);
	    fI[i].resize(numStates + 1, 0.0);   // Extra position for 'gap' entry
        }
    for (int i=0; i<numLeaves; i++)
        {
        numHomologousEmissionsForClass[ (*signature)[ i ] ]++;
        nodeHomology[ i ] = (*signature)[ i ];
        if ( (*signature)[i] == 0)
            numHomologousEmissions[ i ] = 0;
        else
            numHomologousEmissions[ i ] = 1;
        }
        
    // Loop over all nodes except root, and find out which nodes need carry nucleotides of what
    // homology class.  Also fill iChild* arrays, which point to the two children of every parent.
    bool isHomologyClashing = false;
    for (int i=0; i<numNodes-1; i++)
        {
        if ((numHomologousEmissions[i] == numHomologousEmissionsForClass[ nodeHomology[i] ]) || (nodeHomology[i] == 0))
            {
            // This node is the MRCA of the homologous nucleotides iHom[i], or a gap - do nothing
            }
        else
            {
            if (nodeHomology[ parents[ i ]] == 0 || nodeHomology[ parents[ i ]] == nodeHomology[i])
                {
                // This node is not yet MRCA, so node above must be homologous
                nodeHomology[parents[i]] = nodeHomology[i];
                // If iHom[i]==0, iHomNum[i]==0.
                numHomologousEmissions[parents[i]] += numHomologousEmissions[i];
                }
            else
                {
                // Clashing homology; signal
                isHomologyClashing = true;
                }
            }
        if (lftChildrenIndices[ parents[ i ] ] == 0)
            lftChildrenIndices[ parents[ i ] ] = i;
        else
            rhtChildrenIndices[ parents[ i ] ] = i;
        }

    // Bail out - cheaper than do this implicitly in the recursion below
    if (isHomologyClashing)
        {
        //std::cout << "Homology clashes!" << std::endl;
        return 0.0;
        }

//    std::cout << *signature << std::endl;
    
    // Start recursion.  First initialise the leaves
    for (int i=0; i<numLeaves; i++)
        {
        if ( (*signature)[i] == 0)
            {
            // gap
            fI[i][numStates] = 1.0;
            }
        else
            {
            // nucleotide
            fH[i][ sequences[i][(*pos)[i]] ] = 1.0;
            }
        }

    // Now do the recursion, bottom-up, on all internal nodes
    for (int i=numLeaves; i<numNodes; i++)
        {
        int lftChildIdx = lftChildrenIndices[i];
        int rhtChildIdx = rhtChildrenIndices[i];
        
        // Find out whether:
        // 1- One homology family spanning tree intersects both child edges, i.e.
        //    a 'homologous' nucleotide should travel down both edges,
        // 2- One homology family Spanning tree intersects one of the child edges, i.e.
        //    'homologous' nucleotide must travel down a specific edge
        // 3- No homology family's spanning tree intersect either child edge, i.e.
        //    a 'homologous' nucleotide may travel down either edge and pop out at a leaf,
        //    or may travel down both edges but has to die in at least one subtree, or
        //    an 'inhomologous' nucleotide may do whatever it likes here.
        if ((nodeHomology[i] != 0) && (nodeHomology[i] == nodeHomology[lftChildIdx]) && (nodeHomology[i] == nodeHomology[rhtChildIdx]))
            {
            // case 1
            for (int j=0; j<numStates; j++)
                {
                double lft = 0.0;
                double rht = 0.0;
                for (int k=0; k<numStates; k++)
                    {
                    lft += fH[lftChildIdx][k] * homologousProbability[lftChildIdx] * transitionProbabilities[lftChildIdx][j][k];
                    rht += fH[rhtChildIdx][k] * homologousProbability[rhtChildIdx] * transitionProbabilities[rhtChildIdx][j][k];
                    }
                fH[i][j] = lft * rht;
                }
                // Others: 0.0
            }
        else if (nodeHomology[i] != 0)
            {
            // case 2.  Figure out which is the homologous child, and which the inhomologous one.
            int homologousChildIdx, inhomologousChildIdx;
            if (nodeHomology[i] == nodeHomology[lftChildIdx])
                {
                homologousChildIdx = lftChildIdx;
                inhomologousChildIdx = rhtChildIdx;
                }
            else
                {
                homologousChildIdx = rhtChildIdx;
                inhomologousChildIdx = lftChildIdx;
                }

            for (int j=0; j<numStates; j++)
                {
                double lft = 0.0;
                double rht = extinctionProbability[inhomologousChildIdx] * fI[inhomologousChildIdx][numStates];
                for (int k=0; k<numStates; k++)
                    {
                    lft += fH[homologousChildIdx][k] * homologousProbability[homologousChildIdx] * transitionProbabilities[homologousChildIdx][j][k];
                    rht +=
                        (fH[inhomologousChildIdx][k] + fI[inhomologousChildIdx][k]) * (nonHomologousProbability[inhomologousChildIdx] - extinctionProbability[inhomologousChildIdx] * birthProbability[inhomologousChildIdx]) * stateEquilibriumFrequencies[k] +
                        fI[inhomologousChildIdx][k] * homologousProbability[inhomologousChildIdx] * transitionProbabilities[inhomologousChildIdx][j][k];
                    }
                fH[i][j] = lft * rht;
                }
            // Others: 0.0
            }
        else
            {
            // case 3
            int child1Idx = lftChildIdx;
            int child2Idx = rhtChildIdx;
            for (int j=0; j<numStates; j++)
                {
                double lft1 = 0.0, lft2 = 0.0;
                double rht1 = extinctionProbability[child1Idx] * fI[child1Idx][numStates];
                double rht2 = extinctionProbability[child2Idx] * fI[child2Idx][numStates];
                for (int k=0; k<numStates; k++)
                    {
                    lft1 += fH[child1Idx][k] * homologousProbability[child1Idx] * transitionProbabilities[child1Idx][j][k];
                    lft2 += fH[child2Idx][k] * homologousProbability[child2Idx] * transitionProbabilities[child2Idx][j][k];
                    rht1 += (fH[child1Idx][k] + fI[child1Idx][k]) * (nonHomologousProbability[child1Idx] - extinctionProbability[child1Idx] * birthProbability[child1Idx]) * stateEquilibriumFrequencies[k] + fI[child1Idx][k] * homologousProbability[child1Idx] * transitionProbabilities[child1Idx][j][k];
                    rht2 += (fH[child2Idx][k] + fI[child2Idx][k]) * (nonHomologousProbability[child2Idx] - extinctionProbability[child2Idx] * birthProbability[child2Idx]) * stateEquilibriumFrequencies[k] + fI[child2Idx][k] * homologousProbability[child2Idx] * transitionProbabilities[child2Idx][j][k];
                    }
                fH[i][j] = lft1 * rht2 + lft2 * rht1;  // homology pops out below iC1 + homology pops out below iC2
                fI[i][j] = rht1 * rht2;                // no homology with j below i.
                }
            double lft = fI[child1Idx][numStates];
            double rht = fI[child2Idx][numStates];
            for (int j=0; j<numStates; j++)
                {
                lft -= birthProbability[child1Idx] * (fI[child1Idx][j] + fH[child1Idx][j]) * stateEquilibriumFrequencies[j];
                rht -= birthProbability[child2Idx] * (fI[child2Idx][j] + fH[child2Idx][j]) * stateEquilibriumFrequencies[j];
                }
            fI[i][numStates] = lft * rht;
            }
            
        }

    // Now calculate the final result
    int rootIdx = numNodes - 1;
    double res = fI[rootIdx][numStates];
    for (int i=0; i<numStates; i++)
        res -= (fI[rootIdx][i] + fH[rootIdx][i]) * birthProbability[rootIdx] * stateEquilibriumFrequencies[i];
    return res;
}

