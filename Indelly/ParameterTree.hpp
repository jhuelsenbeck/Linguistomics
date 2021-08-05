#ifndef ParameterTree_hpp
#define ParameterTree_hpp

#include <string>
#include "Parameter.hpp"
class RandomVariable;
class Tree;



class ParameterTree : public Parameter {

    public:
                                        ParameterTree(void) = delete;
                                        ParameterTree(const ParameterTree& pt) = delete;
                                        ParameterTree(RandomVariable* r, Model* m, std::string treeStr, std::vector<std::string> tNames, double itl);
                                       ~ParameterTree(void);
        void                            accept(void);
        Tree*                           getActiveTree(void) { return trees[0]; }
        std::string                     getHeader(void) { return ""; }
        double                          lnPriorProbability(void);
        void                            print(void);
        std::string                     getString(void);
        void                            reject(void);
        double                          update(void);
                
    private:
        void                            normalize(std::vector<double>& vec, double minVal);
        double                          updateBrlenProportions(void);
        double                          updateSpr(void);
        double                          updateTreeLength(void);
        Tree*                           trees[2];
        double                          betaT;
};

#endif
