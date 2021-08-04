#ifndef Model_hpp
#define Model_hpp

#include <string>
#include <vector>
class Alignment;
class Parameter;
class ParameterAlignment;
class RandomVariable;
class Tree;



class Model {

    public:
                                            Model(void) = delete;
                                            Model(RandomVariable* r);
                                           ~Model(void);
        void                                accept(void);
        std::string                         getLastUpdate(void);
        ParameterAlignment*                 getAlignment(int idx);
        std::vector<ParameterAlignment*>    getAlignments(void);
        double                              getDeletionRate(void);
        std::vector<double>&                getEquilibriumFrequencies(void);
        std::vector<double>&                getExchangabilityRates(void);
        double                              getInsertionRate(void);
        int                                 getNumAlignments(void);
        std::string                         getParameterHeader(void);
        std::string                         getParameterString(void);
        Tree*                               getTree(void);
        std::string                         getUpdatedParameterName(void);
        double                              lnLikelihood(void);
        double                              lnPriorProbability(void);
        void                                reject(void);
        double                              update(void);
    
    private:
        std::vector<std::string>            readWords(std::string fp, std::vector<Alignment*>& wordAlignments);
        void                                wordLnLike(int i, ParameterAlignment* aln, Tree* t);
        RandomVariable*                     rv;
        double*                             threadLnL;
        std::vector<Parameter*>             parameters;
        std::vector<ParameterAlignment*>    wordParameterAlignments;
        int                                 updatedParameterIdx;
};

#endif
