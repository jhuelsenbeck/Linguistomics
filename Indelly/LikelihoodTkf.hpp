#ifndef LikelihoodTkf_H
#define LikelihoodTkf_H

#include <map>
#include <vector>
#include "IntVector.hpp"
class Model;
class Node;
class ParameterAlignment;
class Tree;



class LikelihoodTkf {

    public:
                                                        LikelihoodTkf(void) = delete;
                                                        LikelihoodTkf(LikelihoodTkf& lkf) = delete;
                                                        LikelihoodTkf(ParameterAlignment* a, Tree* t, Model* m, std::string sm);
                                                       ~LikelihoodTkf(void);
        double                                          tkfLike(void);
        
    private:
        void                                            clearDpTable(void);
        void                                            debugPrint(void);
        void                                            init(void);
        void                                            initAlignment(void);
        void                                            initSequences(void);
        void                                            initTKF91(void);
        void                                            initTransitionProbabilities(void);
        void                                            initTree(Tree* t);
        void                                            printTable(void);
        void                                            printTreeInfo(void);
        void                                            printVector(std::string header, std::vector<int>& v);
        void                                            printVector(std::string header, std::vector<double>& v);
        void                                            printVector(std::string header, std::vector< std::vector<double> >& v);
        double                                          treeRecursion(IntVector* signature, IntVector* pos);
        ParameterAlignment*                             data;
        Tree*                                           tree;
        Model*                                          model;
        std::vector<int>                                parents;
        std::vector<double>                             tau;
        std::vector<std::vector<int> >                  alignment;
        std::vector<std::vector<int> >                  sequences;
        std::vector<double**>                           transitionProbabilities;
        std::vector<double>                             stateEquilibriumFrequencies;
        std::vector<double>                             birthProbability;
        std::vector<double>                             extinctionProbability;
        std::vector<double>                             homologousProbability;
        std::vector<double>                             nonHomologousProbability;
        double                                          immortalProbability;
        int                                             numStates;
        double                                          insertionRate;
        double                                          deletionRate;
        std::map<IntVector*,double,CompIntVector>       dpTable;
        enum                                            StateLabels { free, possible, edgeUsed, used };
        static int                                      unalignableRegionSize;
        static int                                      maxUnalignableDimension;
        static constexpr double                         minEdgeLength = 1e-6;
        std::string                                     substitutionModel;
};

#endif
