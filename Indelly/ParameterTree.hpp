#ifndef ParameterTree_hpp
#define ParameterTree_hpp

#include <map>
#include <string>
#include "Parameter.hpp"
class Alignment;
class RandomVariable;
class Tree;



class ParameterTree : public Parameter {

    public:
                                        ParameterTree(void) = delete;
                                        ParameterTree(const ParameterTree& pt) = delete;
                                        ParameterTree(RandomVariable* r, Model* m, std::string treeStr, std::vector<std::string> tNames, double itl);
                                       ~ParameterTree(void);
        void                            accept(void);
        void                            addSubtrees(std::vector<Alignment*>& alns);
        void                            clearSubtrees(void);
        Tree*                           getActiveTree(void) { return trees[0]; }
        Tree*                           getActiveTree(std::string mask);
        std::string                     getHeader(void) { return ""; }
        double                          lnPriorProbability(void);
        void                            print(void);
        std::string                     getString(void);
        void                            reject(void);
        double                          update(void);
                
    private:
        int                             countMaskBits(std::vector<bool>& m);
        void                            nniArea(std::vector<Node*>& backbone, Node*& incidentNode);
        void                            normalize(std::vector<double>& vec, double minVal);
        double                          updateBrlenProportions(void);
        double                          updateBranchlengthsFromPrior(void);
        double                          updateTopologyFromPrior(void);
        double                          updateNni(void);
        double                          updateSpr(void);
        double                          updateTreeLength(void);
        double                          betaT;
        Tree*                           trees[2];
        std::map<std::string,Tree*>     subtrees;
};

#endif
