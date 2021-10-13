#ifndef TransitionProbabilities_hpp
#define TransitionProbabilities_hpp

#include <map>
#include <string>
#include <vector>
#include "EigenSystem.hpp"
#include "RbBitSet.h"
#include "UserSettings.hpp"
class Alignment;
class Model;
class Tree;

struct TransitionProbabilitiesPair {

    std::vector<StateMatrix_t*> probs[2];
};



class TransitionProbabilities {

    public:
        static TransitionProbabilities&                 transitionProbabilties(void) {
                                                            static TransitionProbabilities tp;
                                                            return tp;
                                                        }
        void                                            flipActive(void);
        int                                             getNumNodes(void) { return numNodes; }
        int                                             getNumStates(void) { return numStates; }
        std::vector<double>                             getStationaryFrequencies(void) { return stationaryFreqs[activeProbs]; }
        std::vector<StateMatrix_t*>                     getTransitionProbabilities(RbBitSet& bs);
        StateMatrix_t*                                  getTransitionProbabilities(RbBitSet& bs, int nodeIdx);
        void                                            initialize(Model* m, std::vector<Alignment*>& alns, int nn, int ns, int sm);
        void                                            print(void);
        void                                            setNeedsUpdate(bool tf) { needsUpdate = tf; }
        void                                            setTransitionProbabilities(void);
    
    private:
                                                        TransitionProbabilities(void);
                                                       ~TransitionProbabilities(void);
                                                        TransitionProbabilities(const TransitionProbabilities& tp) = delete;
        void                                            padeTransitionProbabilities(Tree* t, const StateMatrix_t& Q, const std::vector<StateMatrix_t*>& probs);
        void                                            setTransitionProbabilitiesJc69(void);
        void                                            setTransitionProbabilitiesUsingEigenSystem(void);
        void                                            setTransitionProbabilitiesUsingPadeMethod(void);
        bool                                            isInitialized;
        int                                             numNodes;
        int                                             numStates;
        int                                             activeProbs;
        std::map<RbBitSet,TransitionProbabilitiesPair>  transProbs;
        std::vector<double>                             stationaryFreqs[2];
        Model*                                          modelPtr;
        bool                                            needsUpdate;
        int                                             substitutionModel;
        int                                             numRateCategories;
};

#endif
