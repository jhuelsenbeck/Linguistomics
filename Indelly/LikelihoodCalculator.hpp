#ifndef LikelihoodCalculator_hpp
#define LikelihoodCalculator_hpp

#include <map>
#include <set>
#include <vector>
#include "IntVector.hpp"
#include "RbBitSet.h"
#include "TransitionProbabilities.hpp"
class Model;
class ParameterAlignment;
class TransitionProbabilities;
class Tree;
typedef std::map<IntVector*, double, CompIntVector> PartialProbabilitiesLookup;



class LikelihoodCalculator {

    public:
                                        LikelihoodCalculator(void) = delete;
                                        LikelihoodCalculator(ParameterAlignment* a, Model* m);
                                       ~LikelihoodCalculator(void);
        double                          lnLikelihood(void);
    
    private:
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
        RbBitSet                        taxonMask;
        ParameterAlignment*             data;
        Model*                          model;
        Tree*                           tree;
        std::vector<IntVector*>         pool;
        std::set<IntVector*>            allocated;
        std::vector<std::vector<int> >  alignment;
        std::vector<std::vector<int> >  sequences;
        int                             numStates;
        int                             numTaxa;
        int                             numNodes;
        int                             unalignableRegionSize;
        static int                      maxUnalignableDimension;
        static constexpr double         minBranchLength = 1e-6;
        enum                            StateLabels { free, possible, edgeUsed, used };
        PartialProbabilitiesLookup      partialProbabilities;
        TransitionProbabilities*        transitionProbabilityFactory;
        std::vector<StateMatrix_t*>     transitionProbabilities;
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
        double*                         zeroH;
        double*                         zeroI;
        double                          immortalProbability;
};

#endif
