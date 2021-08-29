#ifndef Model_hpp
#define Model_hpp

#include <map>
#include <set>
#include <string>
#include <vector>
#include "json.hpp"
class Alignment;
class Parameter;
class ParameterAlignment;
class ParameterTree;
class RandomVariable;
class Tree;
class UserSettings;



class Model {

    public:
                                                Model(void) = delete;
                                                Model(RandomVariable* rj);
                                               ~Model(void);
        void                                    accept(void);
        std::string                             getLastUpdate(void);
        ParameterAlignment*                     getAlignment(int idx);
        std::vector<ParameterAlignment*>        getAlignments(void);
        std::vector<std::string>&               getCanonicalTaxonList(void) { return canonicalTaxonList; }
        double                                  getDeletionRate(void);
        std::vector<double>&                    getEquilibriumFrequencies(void);
        std::vector<double>&                    getExchangabilityRates(void);
        std::vector<double>&                    getIndelGammaRates(void);
        double                                  getInsertionRate(void);
        int                                     getNumAlignments(void);
        std::string                             getParameterHeader(void);
        std::string                             getParameterString(void);
        ParameterTree*                          getParameterTree(void);
        std::string                             getStateSetsJsonString(void);
        Tree*                                   getTree(void);
        Tree*                                   getTree(std::string mask);
        std::string                             getUpdatedParameterName(void);
        double                                  lnLikelihood(void);
        double                                  lnPriorProbability(void);
        void                                    reject(void);
        double                                  update(void);
    
    private:
        std::vector<Alignment*>                 initializeAlignments(nlohmann::json& j);
        void                                    initializeParameters(std::vector<Alignment*>& wordAlignments, nlohmann::json& j);
        void                                    initializeStateSets(nlohmann::json& j);
        void                                    initializeTransitionProbabilities(int numStates, nlohmann::json& j);
        nlohmann::json                          parseJsonFile(void);
        void                                    wordLnLike(int i, ParameterAlignment* aln, Tree* t);
        RandomVariable*                         rv;
        double*                                 threadLnL;
        std::vector<Parameter*>                 parameters;
        std::vector<ParameterAlignment*>        wordParameterAlignments;
        int                                     updatedParameterIdx;
        int                                     substitutionModel;
        std::map<std::string,std::set<int> >    stateSets;
        std::vector<std::string>                canonicalTaxonList;
};

#endif
