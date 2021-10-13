#ifndef ParameterTree_hpp
#define ParameterTree_hpp

#include <map>
#include <string>
#include "Parameter.hpp"
#include "RbBitSet.h"
class Alignment;
class RandomVariable;
class Tree;

struct TreePair {

    Tree*   trees[2];
};



class ParameterTree : public Parameter {

    public:
                                        ParameterTree(void) = delete;
                                        ParameterTree(const ParameterTree& pt) = delete;
                                        ParameterTree(RandomVariable* r, Model* m, std::string treeStr, std::vector<std::string> tNames, std::vector<Alignment*>& wordAlignments, double itl);
                                       ~ParameterTree(void);
        void                            accept(void);
        void                            clearSubtrees(void);
        Tree*                           getActiveTree(void) { return fullTree.trees[0]; }
        Tree*                           getActiveTree(RbBitSet& mask);
        std::map<RbBitSet,TreePair>&    getSubtrees(void) { return subTrees; }
        std::string                     getHeader(void) { return ""; }
        double                          lnPriorProbability(void);
        void                            print(void);
        void                            printNewick(void);
        std::string                     getString(void);
        void                            reject(void);
        double                          update(void);
        void                            updateSubtrees(void);
                
    private:
        bool                            checkSubtreeCompatibility(double tolerance);
        int                             countMaskBits(std::vector<bool>& m);
        void                            initializeSubtrees(std::vector<Alignment*>& alns);
        void                            nniArea(std::vector<Node*>& backbone, Node*& incidentNode);
        void                            normalize(std::vector<double>& vec, double minVal);
        double                          updateBrlenProportions(void);
        double                          updateBranchlengthsFromPrior(void);
        double                          updateTopologyFromPrior(void);
        double                          updateNni(void);
        double                          updateSpr(void);
        double                          updateTreeLength(void);
        double                          betaT;
        TreePair                        fullTree;
        std::map<RbBitSet,TreePair>     subTrees;
};

#endif
