#include <cmath>
#include <iomanip>
#include <iostream>
#include "AlignmentProposal.hpp"
#include "Container.hpp"
#include "IntVector.hpp"
#include "Model.hpp"
#include "Msg.hpp"
#include "ParameterAlignment.hpp"
#include "RandomVariable.hpp"
#include "TransitionProbabilities.hpp"
#include "Tree.hpp"

/* This code proposes a new alignment and is taken from Beast. This code was originally
   in Java and is defunct in Beast, as far as I can tell. */
   
int AlignmentProposal::bigUnalignableRegion = 0;



AlignmentProposal::AlignmentProposal(ParameterAlignment* a, Tree* t, RandomVariable* r, Model* m, double expnt, double gp) {

    if (expnt < 1.0)
        Msg::error("Exponenet parameter must be greater than 1");
    if (gp >= 0.0)
        Msg::error("Gap penalty parameter must be less than 0");
        
    // addresses of important objects
    alignmentParm = a;
    tree = t;
    rv = r;
    model = m;
    
    // initialize instance variables
    numTaxa = t->getNumTaxa();
    numNodes = (int)tree->getAncestorIndices().size();
    taxonMask = alignmentParm->getTaxonMask();
    gap = gp;
    basis = expnt;
    gapCode = alignmentParm->getGapCode();
    maxlength = 20; // was 1000
    maxUnalignDimension = 17;
        
    // allocate the profile for this alignment
    profile = new int**[numNodes];

    for (int i=0; i<numNodes; i++)
        {
        profile[i] = new int*[numNodes];
        profile[i][0] = new int[numNodes*maxlength];
        for (int j=1; j<numNodes; j++)
            profile[i][j] = profile[i][j-1] + maxlength;
        for (int j = 0; j < numNodes; j++)
            memset(profile[i][j], 0, maxlength * sizeof(int));
    }
                
    possibles.resize(maxUnalignDimension);
    state.resize(3*maxUnalignDimension);
    numStates = alignmentParm->getNumStates();

    dp = new int*[maxlength];
    dp[0] = new int[maxlength*maxlength];
    for (int i=1; i<maxlength; i++)
        dp[i] = dp[i-1] + maxlength;
    for (int i=0; i<maxlength; i++)
        for (int j=0; j<maxlength; j++)
            dp[i][j] = 0;
        
    scoring = new double*[numStates];
    scoring[0] = new double[numStates*numStates];
    for (int i=1; i<numStates; i++)
        scoring[i] = scoring[i-1] + numStates;
    for (int i=0; i<numStates; i++)
        for (int j=0; j<numStates; j++)
            scoring[i][j] = 0.0;
            
    xProfile = new int[numNodes];
    yProfile = new int[numNodes];
    for (int i=0; i<numNodes; i++)
        {
        xProfile[i] = 0;
        yProfile[i] = 0;
        }
}

AlignmentProposal::~AlignmentProposal(void) {

    freeProfile(profile, numNodes);
    
    delete [] scoring[0];
    delete [] scoring;
    delete [] xProfile;
    delete [] yProfile;
}

void AlignmentProposal::cleanTable(std::map<IntVector*, int, CompIntVector>& m) {

    for (std::map<IntVector*,int,CompIntVector>::iterator it = m.begin(); it != m.end(); it++)
        returnToPool(it->first);
    m.clear();
}

int AlignmentProposal::countPaths(std::vector<std::vector<int> >& inputAlignment, int startCol, int endCol) {
		
    // size of alignment
    int len = endCol - startCol + 1;

    // get a numSites X numTaxa matrix containing the pattern of indels
    alignmentParm->getIndelMatrix(inputAlignment, alignment);

    // Enter first count into DP table
    std::map<IntVector*,int,CompIntVector> dpTable;
    IntVector* pos = getVector();
    IntVector* newVec = getVector(*pos);
    dpTable.insert( std::make_pair( newVec, 1) );  // copy of the pos object is inserted
		
    // Array of possible vector indices, used in inner loop
    for (int i=0; i<maxUnalignDimension; i++)
        possibles[i] = 0;
    state.resize(len);
    for (int i=0; i<len; i++)
        state[i] = 0;
    do
        {
        // Find all possible vectors from current position, iPos
        //IntVector mask(numLeaves);
        IntVector* mask = getVector();
        int ptr, numPossible = 0, firstNotUsed = 0;
        for (ptr = firstNotUsed; mask->zeroEntry() && ptr<len; ptr++)
            {
            if (state[ ptr ] != used)
                {
                if (mask->innerProduct( alignment[ptr] ) == 0)
                    {
                    state[ ptr ] = possible;
                    if (numPossible == maxUnalignDimension)
                        {
                        // This gets too hairy - bail out.
                        bigUnalignableRegion++;
                        returnToPool(pos);
                        cleanTable(dpTable);
                        return -1;
                        }
                    possibles[numPossible++] = ptr;
                    }
                mask->add( alignment[ptr] );
                }
            }
        returnToPool(mask);
        
        // Loop over all combinations of possible vectors, which define edges from
        // iPos to another possible position, by ordinary binary counting.
        IntVector* newPos = getVector(*pos);
        IntVector* signature = getVector();
        int posPtr;
        bool unusedPos;
        bool foundNonZero;
        do
            {
            // Find next combination
            foundNonZero = false;
            for (posPtr = numPossible - 1; posPtr >= 0; --posPtr)
                {
                int curPtr = possibles[ posPtr ];
                if (state[ curPtr ] == possible)
                    {
                    state[ curPtr ] = edgeUsed;
                    newPos->add( alignment[ curPtr ] );
                    // Compute signature vector
                    signature->addMultiple( alignment[ curPtr ], posPtr+1 );
                    // Signal: non-zero combination found, and stop
                    foundNonZero = true;
                    posPtr = 0;
                    }
                else
                    {
                    // It was eEdgeUsed (i.e., digit == 1), so reset digit and continue
                    state[ curPtr ] = possible;
                    newPos->subtract( alignment[ curPtr ] );
                    signature->addMultiple( alignment[ curPtr ], -posPtr-1 );
                    }
                }
        
            if (foundNonZero)
                {
                std::map<IntVector*, int, CompIntVector>::iterator it = dpTable.find(pos);
                if (it == dpTable.end())
                    Msg::error("Could not find pos in dpTable map.");
                int left = it->second;
                int right;
                it = dpTable.find(newPos);
                if (it == dpTable.end())
                    {
                    right = 0;
                    unusedPos = true;
                    }
                else
                    {
                    right = it->second;
                    unusedPos = false;
                    }
                right += left;

                // If we are storing a value at a previously unused position, make sure we use a fresh key object
                if (unusedPos)
                    {
                    IntVector* v = getVector(*newPos);
                    dpTable.insert( std::make_pair(v, right) );
                    }
                else
                    {
                    std::map<IntVector*, int, CompIntVector>::iterator it2 = dpTable.find(newPos);
                    if (it2 == dpTable.end())
                        Msg::error("We should have found newPos in dpTable for the alignment proposal!");
                    it2->second = right;
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
                pos->subtract( alignment[ptr] );
                state[ptr] = free;
                }
            --ptr;
            }
        
        if (ptr == -1)
            {
            // No more unused vectors, so we also fell through the edge loop above,
            // hence iNewPos contains the final position
            std::map<IntVector*, int, CompIntVector>::iterator it = dpTable.find(newPos);
            if (it == dpTable.end())
                Msg::error("Could not find newPos in dpTable map.");
            int result = it->second;
            returnToPool(pos);
            returnToPool(newPos);
            returnToPool(signature);
            cleanTable(dpTable);
            return result;
            }
        else
            {
            returnToPool(newPos);
            returnToPool(signature);
            }
        
        // Now use this farthest-out possible vector
        state[ptr] = used;
        pos->add( alignment[ptr] );
        if (ptr <= firstNotUsed)
            firstNotUsed++;
        
        } while (true);

    //returnToPool(pos);
    //cleanTable(dpTable);
    return 0;
}

void AlignmentProposal::drainPool(void) {

    for (std::vector<IntVector*>::iterator v=pool.begin(); v != pool.end(); v++)
        {
        allocated.erase(*v);
        delete (*v);
        }
}

IntVector* AlignmentProposal::getVector(void) {

    if (pool.empty() == true)
        {
        /* If the vector pool is empty, we allocate a new vector and return it. We
           do not need to add it to the vector pool. */
        IntVector* v = new IntVector(numTaxa);
        allocated.insert(v);
        return v;
        }
    
    // Return a vector from the vector pool, remembering to remove it from the pool.
    IntVector* v = pool.back();
    pool.pop_back();
    return v;
}

IntVector* AlignmentProposal::getVector(IntVector& vec) {

    if (pool.empty() == true)
        {
        IntVector* v = new IntVector(numTaxa);
        allocated.insert(v);
        for (int i=0; i<vec.size(); ++i)
            (*v)[i] = vec[i];
        return v;
        }
    IntVector* v = pool.back();
    pool.pop_back();
    for (int i=0; i<vec.size(); ++i)
        (*v)[i] = vec[i];
    return v;
}

void AlignmentProposal::freeProfile(int*** x, int n) {

    for (int i=0; i<n; i++)
        {
        delete [] x[i][0];
        delete [] x[i];
        }
    delete [] x;
}

void AlignmentProposal::initialize(void) {

    // get a pointer to the tree
    tree = model->getTree(taxonMask);
    
    // initialize profile
    for (int i=0; i<numNodes; i++)
        for (int j=0; j<numNodes; j++)
            for (int k=0; k<maxlength; k++)
                profile[i][j][k] = 0;
                
    // check that we have no stray IntVectors
    if (allocated.size() != pool.size())
        {
        std::cout << "allocated.size() = " << allocated.size() << std::endl;
        std::cout << "pool.size()      = " << pool.size() << std::endl;
        Msg::error("Expecting all IntVectors to have been returned to the pool");
        }
}

void AlignmentProposal::print(std::string s, std::vector<int>& x) {

    std::cout << s << std::endl;
    for (int i=0; i<x.size(); i++)
        std::cout << x[i] << " ";
    std::cout << std::endl;
}

void AlignmentProposal::print(std::string s, std::vector<std::vector<int> >& x) {

    std::cout << s << std::endl;
    for (int i=0; i<x.size(); i++)
        {
        for (int j=0; j<x[i].size(); j++)
            std::cout << x[i][j] << " ";
        std::cout << std::endl;
        }
}

void AlignmentProposal::print(std::string s, int*** x, int a, int b, int c) {

    std::cout << s << std::endl;
    for (int i=0; i<a; i++)
        {
        std::cout << i << ":" << std::endl;
        for (int j=0; j<b; j++)
            {
            for (int k=0; k<c; k++)
                std::cout << x[i][j][k] << " ";
            std::cout << std::endl;
            }
        }
}

double AlignmentProposal::propose(std::vector<std::vector<int> >& newAlignment, double extensionProb) {
    
    // check that the extenstion probability is between 0 and 1
    if (extensionProb <= 0.0 || extensionProb > 1.0)
        Msg::error("Extension parameter must be in the range (0,1]");
        
    // initialize variables that may have changed since last update
    initialize();
    std::vector<std::vector<int> >& curAlignment = alignmentParm->getAlignment();

    // initialize tree
    std::vector<int> parents = tree->getAncestorIndices();
    std::vector<int> lftChildrenIndices(numNodes);
    std::vector<int> rhtChildrenIndices(numNodes);
    for (int i=0; i<numNodes-1; i++)
        {
        if (lftChildrenIndices[ parents[ i ] ] == 0)
            lftChildrenIndices[ parents[ i ] ] = i;
        else
            rhtChildrenIndices[ parents[ i ] ] = i;
        }

    // sort children
    std::vector<int> sorter(numNodes);
    for (int i=0; i<curAlignment.size(); i++)
        sorter[i] = i;
    for( int i=(int)curAlignment.size(); i<parents.size(); i++)
        {
        if (sorter[lftChildrenIndices[i]] > sorter[rhtChildrenIndices[i]])
            sorter[i] = sorter[rhtChildrenIndices[i]];
        else
            sorter[i] = sorter[lftChildrenIndices[i]];
        // for internal nodes, sorter gets the lower value
        if(sorter[lftChildrenIndices[i]] > sorter[rhtChildrenIndices[i]])
            {
            // lftChildrenIndices[i] has the smaller sorter value
            // which guarantees that the root's profile contains the sequences in the proper order
            int temp = lftChildrenIndices[i];
            lftChildrenIndices[i] = rhtChildrenIndices[i];
            rhtChildrenIndices[i] = temp;
            }
        }

    // cutting a part of the alignment (we cut out the entire alignment)
#   if 0
    int len = (int)curAlignment[0].size();
    int pos = 0;
#   else
    int len = 0;
    for (len=1; len < curAlignment[0].size() && rv->uniformRv() < extensionProb; len++)
        ;
    int pos = rv->uniformRvInt(curAlignment[0].size() - len + 1);
#   endif

    //std::vector<std::vector<int> > profileNumber;
    profileNumber.resize(numNodes);
    for (int i=0; i<numNodes; i++)
        {
        profileNumber[i].resize(curAlignment.size());
        xProfile[i] = 0;
        yProfile[i] = 0;
        }
    for (int i=0; i<numNodes; i++)
        {
        xProfile[i] = 0;
        yProfile[i] = 0;
        }
    for (int i=0; i<curAlignment.size(); i++)
        {
        yProfile[i] = 1;
        for (int j=pos; j<pos+len; j++)
            {
            if (curAlignment[i][j] != gapCode)
                {
                profile[i][0][xProfile[i]] = curAlignment[i][j];
                xProfile[i]++;
                }
            }
        // profile number
        profileNumber[i][0] = i;
        }
    
    /* dynamic programming algorithm for the proposal */
    TransitionProbabilities& tip = TransitionProbabilities::transitionProbabilties();
    DoubleMatrix** tiProbs = tip.getTransitionProbabilities(taxonMask);
    std::vector<double>& stationaryFrequencies = tip.getStationaryFrequencies();
    for (int i=0; i<numStates; i++)
        for (int j=0; j<numStates; j++)
            scoring[i][j] = 0.0;

    double proposalLoglikelihood = 0.0;
    //std::vector<std::vector<int> > tempProfile;
    tempProfile.resize(numNodes);
    for (int i=0; i<numNodes; i++)
        tempProfile[i].resize(maxlength);
    for (int k=(int)curAlignment.size(); k<numNodes; k++)
        {
        int lftChild = lftChildrenIndices[k];
        int rhtChild = rhtChildrenIndices[k];
        
        // obtaining the scoring matrix
        for (int i=0; i<numStates; i++)
            {
            for (int j=0; j<numStates; j++)
                {
                scoring[i][j] = 0.0;
                for (int l=0; l<numStates; l++)
                    scoring[i][j] += (*tiProbs[lftChild])(i,l) * (*tiProbs[rhtChild])(l,j);
                }
            }
        for (int i=0; i<numStates; i++)
            {
            for (int j=0; j<numStates; j++)
                {
                scoring[i][j] = log(scoring[i][j] / stationaryFrequencies[j]);
                }
            }

        // initialising the table
        dp[0][0] = 0.0;
        for (int i=1; i<=xProfile[lftChild]; i++)
            {
            dp[i][0] = dp[i - 1][0] + gap;
            }
        for (int j=1; j<=xProfile[rhtChild]; j++)
            {
            dp[0][j] = dp[0][j - 1] + gap;
            }

        // main dynamic programming
        double iTempdouble3 = 0.0;
        for (int i=1; i<=xProfile[lftChild]; i++)
            for (int j=1; j<=xProfile[rhtChild]; j++)
                {
                double iTempdouble1 = dp[i - 1][j] + gap;
                double iTempdouble2 = dp[i][j - 1] + gap;
//                { //JPH why restrict the scope?
                int iCounter = 0;
                iTempdouble3 = 0.0;
                for (int l=0; l<yProfile[lftChild]; l++)
                    for (int m=0; m<yProfile[rhtChild]; m++)
                        {
                        if (profile[lftChild][l][i - 1] != gapCode && profile[rhtChild][m][j - 1] != gapCode)
                            {
                            iTempdouble3 += scoring[profile[lftChild][l][i - 1]][profile[rhtChild][m][j - 1]];
                            iCounter++;
                            }
                        }
                    iTempdouble3 = dp[i - 1][j - 1] + iTempdouble3/iCounter;
//                }

                if (iTempdouble1 > iTempdouble2)
                    dp[i][j] = iTempdouble1;
                else
                    dp[i][j] = iTempdouble2;
                if (iTempdouble3 > dp[i][j])
                    dp[i][j] = iTempdouble3;
                }

        // Stochastic traceback, generating the new profile
        int iI = xProfile[lftChild];
        int jJ = xProfile[rhtChild];
        double iTempdouble4;
        while (iI > 0 || jJ > 0)
            {
            double iTempdouble1 = 0.0, iTempdouble2 = 0.0;
            if (iI > 0)
                {
                iTempdouble1 = dp[iI - 1][jJ] + gap;
                }
            if (jJ > 0)
                {
                iTempdouble2 = dp[iI][jJ - 1] + gap;
                }
            if (iI > 0 && jJ > 0)
                {
                int iCounter = 0;
                double iTempdouble3 = 0.0;
                for (int l=0; l<yProfile[lftChild]; l++)
                    for (int m=0; m<yProfile[rhtChild]; m++)
                        {
                        if (profile[lftChild][l][iI - 1] != gapCode && profile[rhtChild][m][jJ - 1] != gapCode)
                            {
                            iTempdouble3 += scoring[profile[lftChild][l][iI - 1]][profile[rhtChild][m][jJ - 1]];
                            iCounter++;
                            }
                        }
                iTempdouble3 = dp[iI - 1][jJ - 1] + iTempdouble3/iCounter;
                }

            // normalising
            iTempdouble4 = 0.0;
            if (iI > 0)
                {
                iTempdouble4 += exp(iTempdouble1 * log(basis));
                }
            if (jJ > 0)
                {
                iTempdouble4 += exp(iTempdouble2 * log(basis));
                }
            if (iI > 0 && jJ > 0)
                {
                iTempdouble4 += exp(iTempdouble3 * log(basis));
                }

            if(iI > 0)
                iTempdouble1 = exp(iTempdouble1 * log(basis))/iTempdouble4;
            else
                iTempdouble1 = 0.0;

            if(jJ > 0)
                iTempdouble2 = exp(iTempdouble2 * log(basis))/iTempdouble4;
            else
                iTempdouble2 = 0.0;

            if(iI > 0 && jJ > 0)
                iTempdouble3 = exp(iTempdouble3 * log(basis))/iTempdouble4;
            else
                iTempdouble3 = 0.0;

            // stochastic step
            iTempdouble4 = rv->uniformRv();
            if (iTempdouble4 < iTempdouble1)
                {
                // the probability is in the denominator, that's why we substract
                proposalLoglikelihood -= log(iTempdouble1);
                for (int l=0; l<yProfile[lftChild]; l++)
                    {
                    tempProfile[l][xProfile[k]] = profile[lftChild][l][iI - 1];
                    }
                for (int l=yProfile[lftChild]; l<yProfile[lftChild] + yProfile[rhtChild]; l++)
                    {
                    tempProfile[l][xProfile[k]] = gapCode;
                    }
                xProfile[k]++;
                iI--;
                }
            else if (iTempdouble4 < iTempdouble1 + iTempdouble2)
                {
                // the probability is in the denominator, that's why we substract
                proposalLoglikelihood -= log(iTempdouble2);
                for (int l=0; l<yProfile[lftChild]; l++)
                    {
                    tempProfile[l][xProfile[k]] = gapCode;
                    }
                for (int l=yProfile[lftChild]; l < yProfile[lftChild] + yProfile[rhtChild]; l++)
                    {
                    tempProfile[l][xProfile[k]] = profile[rhtChild][l - yProfile[lftChild]][jJ - 1];
                    }
                xProfile[k]++;
                jJ--;
                }
            else
                {
                if (iTempdouble3 <= 0.0)
                    {
                    freeProfile(profile, numNodes);
                    Msg::error("Problem proposing new alignment");
                    }
                // the probability is in the denominator, that's why we substract
                proposalLoglikelihood -= log(iTempdouble3);
                for (int l=0; l<yProfile[lftChild]; l++)
                    tempProfile[l][xProfile[k]] = profile[lftChild][l][iI - 1];
                for (int l=yProfile[lftChild]; l < yProfile[lftChild] + yProfile[rhtChild]; l++)
                    tempProfile[l][xProfile[k]] = profile[rhtChild][l - yProfile[lftChild]][jJ - 1];
                xProfile[k]++;
                jJ--;
                iI--;
                }

            } // finished the traceback, iI, jJ

        yProfile[k] = yProfile[lftChild] + yProfile[rhtChild];

        //Making new profileNumbers
        for (int i=0; i<yProfile[lftChild]; i++)
            profileNumber[k][i] = profileNumber[lftChild][i];
        for (int i=yProfile[lftChild]; i<yProfile[lftChild] + yProfile[rhtChild]; i++)
            profileNumber[k][i] = profileNumber[rhtChild][i - yProfile[lftChild]];

        // now reversing the profile, loading from the iTempProfile into iProfile
        for (int i=0; i<xProfile[k]; i++)
            {
            for (int j=0; j<yProfile[k]; j++)
                {
                profile[k][j][i] = tempProfile[j][xProfile[k] - i - 1];
                }
            }

        } // k, index for the internal nodes, profile[root] is the obtained multiple alignment
    
    int iL2 = xProfile[numNodes - 1];

    // putting the new sub-alignment into the new proposal
//    std::vector<std::vector<int> > returnedAlignment(curAlignment.size());
    newAlignment.resize( curAlignment.size() );
    for (int i=0; i < curAlignment.size(); i++)
        {
        newAlignment[i].resize( curAlignment[0].size() + iL2 - len );
        }

    for (int i=0; i<pos; i++)
        {
        for (int j=0; j<curAlignment.size(); j++)
            {
            newAlignment[j][i] = curAlignment[j][i];
            }
        }
    for (int i=pos; i<pos + iL2; i++)
        {
        for (int j=0; j<curAlignment.size(); j++)
            {
            newAlignment[profileNumber[numNodes - 1][j]][i] = profile[numNodes - 1][j][i - pos];
            }
        }

    for (int i=pos + len; i<curAlignment[0].size(); i++)
        {
        for (int j=0; j<curAlignment.size(); j++)
            {
            newAlignment[j][i + iL2 - len] = curAlignment[j][i];
            }
        }
        
//    print("curAlignment", curAlignment);
//    print("newAlignment", newAlignment);
    
    // I involve the combinatorical factor into the Hastings ratio
    // countPaths(int[][] iInputAlignment, int iStartCol, int iEndCol)
//    int iCF = countPaths(returnedAlignment, pos, pos+iL2-1);
    int iCF = countPaths(newAlignment, pos, pos+iL2-1);
    if (iCF == -1)
        {
        freeProfile(profile, numNodes);
        Msg::error("Something is seriously wrong with the alignment");
        proposalLoglikelihood -= exp(100.0);
        iCF = 1;
        }
    proposalLoglikelihood -= log((double)iCF);

    iCF = countPaths(curAlignment, pos, pos+len-1);
    if (iCF == -1)
        {
        freeProfile(profile, numNodes);
        Msg::error("Something is seriously wrong with the alignment");
        proposalLoglikelihood -= exp(100.0);
        iCF = 1;
        }
    proposalLoglikelihood += log((double)iCF);

    /*
    Now I'm going to calculate the Hastings ratio.
    To do that, I must make the profile of each node
    I copy the multiple alignment into profiles, and then delete the all-gaps columns
    The traceback is easy: based on checking all-gap subcolumns
    */

    // collecting the sequences from the alignment -- now without cut out!
    // they will be the first profiles
    for (int i=0; i<curAlignment.size(); i++)
        {
        yProfile[i] = 1;
        for (int j = pos; j <pos + len; j++)
            profile[i][0][j - pos] = curAlignment[i][j];
        }

    // profiles at internal nodes are obtained by merging profiles of children
    for (int i=(int)curAlignment.size(); i<numNodes; i++)
        {
        yProfile[i] = yProfile[lftChildrenIndices[i]] + yProfile[rhtChildrenIndices[i]];
        for(int j=0; j<len; j++)
            {
            for (int k = 0; k < yProfile[lftChildrenIndices[i]]; k++)
                profile[i][k][j] = profile[lftChildrenIndices[i]][k][j];
            for (int k = yProfile[lftChildrenIndices[i]]; k < yProfile[lftChildrenIndices[i]] + yProfile[rhtChildrenIndices[i]]; k++)
                profile[i][k][j] = profile[rhtChildrenIndices[i]][k - yProfile[lftChildrenIndices[i]]][j];
            }
        }

    // and now -- cutting out all-gap columns
    for (int i=0; i<numNodes; i++)
        {
        // first i load it into the temporary profile and then back...
        for (int j=0; j<len; j++)
            for (int k=0; k < yProfile[i]; k++)
                tempProfile[k][j] = profile[i][k][j];
        xProfile[i] = 0;
        for (int j = 0; j < len; j++)
            {
            // checking for all-gap;
            bool notAllGaps = false;
            for (int k = 0; k < yProfile[i]; k++)
                {
                if (tempProfile[k][j] != gapCode)
                    {
                    notAllGaps = true;
                    break;
                    }
                }
            if (notAllGaps)
                {
                // if not all-gap then
                for (int k = 0; k < yProfile[i]; k++)
                    profile[i][k][xProfile[i]] =  tempProfile[k][j];
                xProfile[i]++;
                }
            }
        }

    // profiles are ready, so now the DP
    for (int k = (int)curAlignment.size(); k < numNodes; k++)
        {
        // obtaining the scoring matrix
        for (int i = 0; i < numStates; i++)
            for (int j = 0; j < numStates; j++)
                {
                scoring[i][j] = 0.0;
                for (int l = 0; l < numStates; l++)
                    scoring[i][j] += (*tiProbs[lftChildrenIndices[k]])(i,l) * (*tiProbs[rhtChildrenIndices[k]])(l,j);
                }
        for (int i = 0; i < numStates; i++)
            for (int j = 0; j < numStates; j++)
                scoring[i][j] = log(scoring[i][j] / stationaryFrequencies[j]);

        // initialising the table
        //std::vector<std::vector<double> > dp;
        //dp.resize(maxlength);
        //for (int i=0; i<maxlength; i++)
        //    dp[i].resize(maxlength);
        dp[0][0] = 0.0;
        double iTempdouble2 = 0.0;
        for (int i = 1; i <= xProfile[lftChildrenIndices[k]]; i++)
            {
            dp[i][0] = dp[i - 1][0] + gap;
            }
        for (int j = 1; j <= xProfile[rhtChildrenIndices[k]]; j++)
            {
            iTempdouble2 = 0.0;
            for (int m = 0; m < yProfile[rhtChildrenIndices[k]]; m++)
                if(profile[rhtChildrenIndices[k]][m][j - 1] != gapCode)
                    iTempdouble2 += gap;
            dp[0][j] = dp[0][j - 1] + gap;
            }

        // main dynamic programming
        double iTempdouble1 = 0.0, iTempdouble3 = 0.0;
        for (int i = 1; i <= xProfile[lftChildrenIndices[k]]; i++)
            for (int j = 1; j <= xProfile[rhtChildrenIndices[k]]; j++)
                {
                iTempdouble1 = dp[i - 1][j] + gap;

                iTempdouble2 = 0.0;
                for (int m = 0; m < yProfile[rhtChildrenIndices[k]]; m++)
                    if (profile[rhtChildrenIndices[k]][m][j - 1] != gapCode)
                        iTempdouble2 += gap;
                iTempdouble2 = dp[i][j - 1] + gap;

                {
                int iCounter = 0;
                iTempdouble3 = 0.0;
                for (int l = 0; l < yProfile[lftChildrenIndices[k]]; l++)
                    for (int m = 0; m < yProfile[rhtChildrenIndices[k]]; m++)
                        {
                        if (profile[lftChildrenIndices[k]][l][i - 1] != gapCode && profile[rhtChildrenIndices[k]][m][j - 1] != gapCode)
                            {
                            iTempdouble3 += scoring[profile[lftChildrenIndices[k]][l][i - 1]][profile[rhtChildrenIndices[k]][m][j - 1]];
                            iCounter++;
                            }
                        }
                iTempdouble3 = dp[i - 1][j - 1] + iTempdouble3/iCounter;
                }

            if (iTempdouble1 > iTempdouble2)
                dp[i][j] = iTempdouble1;
            else
                dp[i][j] = iTempdouble2;
            if (iTempdouble3 > dp[i][j])
                dp[i][j] = iTempdouble3;
            }

        // traceback, now not stochastic
        int i = xProfile[lftChildrenIndices[k]];
        int j = xProfile[rhtChildrenIndices[k]];
        int n = xProfile[k];
        double iTempdouble4;
        while (n > 0)
            {
            if(i > 0)
                {
                iTempdouble1 = dp[i - 1][j] + gap;
                }
            if (j > 0)
                {
                iTempdouble2 = dp[i][j - 1] + gap;
                }
            if (i > 0 && j > 0)
                {
                int iCounter = 0;
                iTempdouble3 = 0.0;
                for (int l = 0; l < yProfile[lftChildrenIndices[k]]; l++)
                    for (int m = 0; m < yProfile[rhtChildrenIndices[k]]; m++)
                        {
                        if (profile[lftChildrenIndices[k]][l][i - 1] != gapCode && profile[rhtChildrenIndices[k]][m][j - 1] != gapCode)
                            {
                            iTempdouble3 += scoring[profile[lftChildrenIndices[k]][l][i - 1]][profile[rhtChildrenIndices[k]][m][j - 1]];
                            iCounter++;
                            }
                        }
                iTempdouble3 = dp[i - 1][j - 1] + iTempdouble3/iCounter;
                }

            // normalising
            iTempdouble4 = 0.0;
            if (i > 0)
                iTempdouble4 += exp(iTempdouble1 * log(basis));
            if (j > 0)
                iTempdouble4 += exp(iTempdouble2 * log(basis));
            if (i > 0 && j > 0)
                iTempdouble4 += exp(iTempdouble3 * log(basis));

            if (i > 0)
                iTempdouble1 = exp(iTempdouble1 * log(basis))/iTempdouble4;
            else
                iTempdouble1 = 0.0;

            if (j > 0)
                iTempdouble2 = exp(iTempdouble2 * log(basis))/iTempdouble4;
            else
                iTempdouble2 = 0.0;

            if (i > 0 && j > 0)
                iTempdouble3 = exp(iTempdouble3 * log(basis))/iTempdouble4;
            else
                iTempdouble3 = 0.0;

            // finding the proper step (all-gap subcolumns tell us the way
            bool iTempBoole = true;
            for (int l = 0; l < yProfile[lftChildrenIndices[k]] && iTempBoole; l++)
                iTempBoole = (profile[k][l][n-1] == gapCode);
            if (iTempBoole)
                {
                // all the first part is gap, so we decrease j
                j--;
                proposalLoglikelihood += log(iTempdouble2);
                n--;
                }
            else
                {
                iTempBoole = true;
                for (int l = yProfile[lftChildrenIndices[k]]; l < yProfile[lftChildrenIndices[k]] + yProfile[rhtChildrenIndices[k]] && iTempBoole; l++)
                    iTempBoole = (profile[k][l][n-1] == gapCode);
                if (iTempBoole)
                    {
                    // all the second part is gap, so we decrease i
                    i--;
                    proposalLoglikelihood += log(iTempdouble1);
                    n--;
                    }
            else
                {
                // no all-gap subcolumn, decreasing both i and j;
                i--;
                j--;
                proposalLoglikelihood += log(iTempdouble3);
                n--;
                }
            }

        } // end of traceback

    } // end of k

#   if 0
    // debug print
    std::cout << alignmentParm->getName() << std::endl;
    std::cout << "original" << std::endl;
    for (int i=0; i<curAlignment.size(); i++)
        {
        for (int j=0; j<curAlignment[i].size(); j++)
            std::cout << std::setw(2) << curAlignment[i][j] << " ";
        std::cout << std::endl;
        }

    std::cout << "proposed" << std::endl;
    for (int i=0; i<newAlignment.size(); i++)
        {
        for (int j=0; j<newAlignment[i].size(); j++)
            std::cout << std::setw(2) << newAlignment[i][j] << " ";
        std::cout << std::endl;
        }
#   endif
        
    return proposalLoglikelihood + (len - iL2) * log(1.0-extensionProb);
}

void AlignmentProposal::returnToPool(IntVector* v) {

    v->clean();
    pool.push_back(v);
}
