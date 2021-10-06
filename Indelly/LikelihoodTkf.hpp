#ifndef LikelihoodTkf_H
#define LikelihoodTkf_H

#include <map>
#include <vector>
#include "mpreal.h"
#include "IntVector.hpp"
#include "RbBitSet.h"
#include "TransitionProbabilities.hpp"
class Model;
class Node;
class ParameterAlignment;
class SiteLikelihood;
class Tree;



class LikelihoodTkf {

    public:
                                                            LikelihoodTkf(void) = delete;
                                                            LikelihoodTkf(LikelihoodTkf& lkf) = delete;
                                                            LikelihoodTkf(ParameterAlignment* a, Model* m);
                                                           ~LikelihoodTkf(void);
        double                                              tkfLike(void);
        
    private:
        void                                                clearDpTable(void);
        void                                                init(void);
        void                                                initAlignment(void);
        void                                                initSequences(void);
        void                                                initTKF91(void);
        void                                                initTransitionProbabilities(void);
        void                                                initTree(void);
        void                                                printTable(void);
        void                                                printVector(std::string header, std::vector<int>& v);
        void                                                printVector(std::string header, std::vector<double>& v);
        void                                                printVector(std::string header, std::vector< std::vector<double> >& v);
        void                                                printVector(std::string header, std::vector< std::vector<int> >& v);
        double                                              treeRecursion(IntVector* signature, IntVector* pos, int siteColumn);
        ParameterAlignment*                                 data;
        Tree*                                               tree;
        Model*                                              model;
        std::vector<std::vector<int> >                      alignment;
        std::vector<std::vector<int> >                      sequences;
        std::vector<StateMatrix_t*>                         transitionProbabilities;
        std::vector<double>                                 stateEquilibriumFrequencies;
        int                                                 numIndelCategories;
        std::vector<double>                                 birthProbability;
        std::vector<double>                                 extinctionProbability;
        std::vector<double>                                 homologousProbability;
        std::vector<double>                                 nonHomologousProbability;
        SiteLikelihood*                                     siteProbs;
        double**                                            fH;
        double**                                            fI;
        double                                              immortalProbability;
        int                                                 numStates;
        double                                              insertionRate;
        double                                              deletionRate;
        RbBitSet                                            taxonMask;
//        std::map<IntVector*,double,CompIntVector>         dpTable;
        std::map<IntVector*,mpfr::mpreal*,CompIntVector>    dpTable;
        enum                                                StateLabels { free, possible, edgeUsed, used };
        static int                                          unalignableRegionSize;
        static int                                          maxUnalignableDimension;
        static constexpr double                             minEdgeLength = 1e-6;
};

#endif
