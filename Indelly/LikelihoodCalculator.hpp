#ifndef LikelihoodCalculator_hpp
#define LikelihoodCalculator_hpp

#include <set>
#include <string>
#include <vector>
#include "IntVector.hpp"
#include "JphMatrix.hpp"
#include "RbBitSet.h"
#include "TransitionProbabilities.hpp"
class IndelMatrix;
class Model;
class Node;
class ParameterAlignment;
class TransitionProbabilities;
class Tree;
typedef std::map<IntVector*, double, CompIntVector> PartialProbabilitiesLookup;

struct IndelProbabilities {

    double      insertionRate;
    double      deletionRate;
    double      immortalProbability;
    double*     beta;
    double*     birthProbability;
    double*     extinctionProbability;
    double*     homologousProbability;
    double*     nonHomologousProbability;
};

struct IndelCombinatorics {

    int         maximumSequenceLength;
    int*        state;
    int*        possibleVectorIndices;
    int*        nodeHomology;
    int*        numHomologousEmissions;
    int*        numHomologousEmissionsForClass;
};



class LikelihoodCalculator {

    public:
        LikelihoodCalculator(void) = delete;
                                        LikelihoodCalculator(ParameterAlignment* a, Model* m);
                                       ~LikelihoodCalculator(void);
        std::string                     alignmentName(void);
        double                          lnLikelihood(void);
    
    private:
        const int                       maxUnalignableDimension  = 10,
                                        maxUnalignableDimension1 = maxUnalignableDimension + 1;
        const double                    minBranchLength = 1e-6;
        void                            allcoateIndelCombinatorics(int nn, int maxSeqLen);
        void                            allocateIndelProbabilities(int nn);
        void                            clearPpTable(void);
        void                            drainPool(void);
        void                            freeIndelCombinatorics(void);
        void                            freeIndelProbabilities(void);
        IntVector*                      getVector(void);
        IntVector*                      getVector(IntVector& vec);
        int                             getNumAllocated(void) { return (int) allocated.size(); }
        int                             getNumInPool(void) { return (int)pool.size(); }
        void                            initialize(void);
        double                          partialProbability(IntVector* signature, IntVector* pos);
        void                            printTable(void);
        void                            returnToPool(IntVector* n);
        void                            setBirthDeathProbabilities(void);
        double                          prune(IntVector* signature, IntVector* pos, std::vector<Node*>& dpSequence);

        RbBitSet                        taxonMask;
        ParameterAlignment*             data;
        Model*                          model;
        Tree*                           tree;
        std::vector<Node*>              des;       // instantiated once on construction and modified via reference
        std::vector<IntVector*>         pool;
        std::set<IntVector*>            allocated;
        IndelMatrix*                    alignment;
        std::vector<std::vector<int> >  sequences;  // instantiated once on construction of object, never modified
        int                             numStates,
                                        numStates1;
        int                             numTaxa;
        int                             numNodes;
        int                             unalignableRegionSize;
        enum                            StateLabels { free, possible, edgeUsed, used };
        PartialProbabilitiesLookup      partialProbabilities;
        TransitionProbabilities*        transitionProbabilityFactory;
        DoubleMatrix*                   transitionProbabilities;
        std::vector<double>             stateEquilibriumFrequencies;
        IndelProbabilities              indelProbs;
        double**                        fH;
        double**                        fI;
        IndelCombinatorics              indelCombos;
};

#endif
