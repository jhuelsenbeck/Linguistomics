#ifndef LikelihoodCalculator_hpp
#define LikelihoodCalculator_hpp

#include <set>
#include <string>
#include <vector>
#include "IntVector.hpp"
#include "JphMatrix.hpp"
#include "RbBitSet.h"
#include "TransitionProbabilities.hpp"
class Model;
class Node;
class ParameterAlignment;
class TransitionProbabilities;
class Tree;
typedef std::map<IntVector*, double, CompIntVector> PartialProbabilitiesLookup;



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

        void                            clearPpTable(void);
        void                            drainPool(void);
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
        std::vector<IntVector*>         pool;
        std::set<IntVector*>            allocated;
        std::vector<std::vector<int> >  alignment;
        std::vector<std::vector<int> >  sequences;
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
        double                          insertionRate;
        double                          deletionRate;
        std::vector<double>             beta;
        std::vector<double>             birthProbability;
        std::vector<double>             extinctionProbability;
        std::vector<double>             homologousProbability;
        std::vector<double>             nonHomologousProbability;
        double**                        fH;
        double**                        fI;
        double                          immortalProbability;
        std::vector<int>                possibleVectorIndices;
        std::vector<int>                state;
        std::vector<int>                nodeHomology;
        std::vector<int>                numHomologousEmissions;
        std::vector<int>                numHomologousEmissionsForClass;
        std::vector<Node*>              des;
};

#endif
