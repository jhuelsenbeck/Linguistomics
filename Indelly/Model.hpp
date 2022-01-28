#ifndef Model_hpp
#define Model_hpp

#include <map>
#include <set>
#include <string>
#include <vector>
#include "json.hpp"
#include "RbBitSet.h"
#include "threads.hpp"
class Alignment;
class LikelihoodCalculator;
class Parameter;
class ParameterAlignment;
class ParameterTree;
class Partition;
class RandomVariable;
class Tree;
class UserSettings;



class Model {

    public:
                                                Model(void) = delete;
                                                Model(RandomVariable* r, ThreadPool* p);
                                               ~Model(void);
        void                                    accept(void);
        void                                    flipActiveLikelihood(void);
        void                                    flipActiveLikelihood(int idx);
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
        Tree*                                   getTree(RbBitSet& mask);
        Tree*                                   getTree(const RbBitSet& mask);
        std::string                             getUpdatedParameterName(void);
        double                                  lnLikelihood(void);
        double                                  lnPriorProbability(void);
        void                                    reject(void);
        void                                    setUpdateLikelihood(void);
        void                                    setUpdateLikelihood(int idx);
        double                                  update(void);
            
    private:
        std::vector<Alignment*>                 initializeAlignments(nlohmann::json& j);
        void                                    initializeParameters(std::vector<Alignment*>& wordAlignments, nlohmann::json& j);
        void                                    initializeStateSets(nlohmann::json& j);
        void                                    initializeTransitionProbabilities(std::vector<Alignment*>& wordAlignments);
        nlohmann::json                          parseJsonFile(void);
        void                                    wordLnLike(int i);
        RandomVariable*                         rv;
        ThreadPool*                             threadPool;
        double*                                 threadLnL;
        std::vector<bool>                       updateLikelihood;
        std::vector<int>                        activeLikelihood;
        std::vector<double>                     wordLnLikelihoods[2];
        std::vector<Parameter*>                 parameters;
        std::vector<ParameterAlignment*>        wordParameterAlignments;
        std::vector<LikelihoodCalculator*>      wordLikelihoodCalculators;
        int                                     updatedParameterIdx;
        int                                     substitutionModel;
        Partition*                              partitionInfo;
        std::map<std::string,std::set<int> >    stateSets;
        std::vector<std::string>                canonicalTaxonList;
};

#endif
