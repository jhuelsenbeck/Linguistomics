#ifndef TransitionProbabilities_hpp
#define TransitionProbabilities_hpp

#include <vector>
class Model;



class TransitionProbabilities {

    public:
        static TransitionProbabilities& transitionProbabilties(void) {
                                            static TransitionProbabilities tp;
                                            return tp;
                                        }
        void                            flipActive(void);
        int                             getNumNodes(void) { return numNodes; }
        int                             getNumStates(void) { return numStates; }
        std::vector<double**>           getTransitionProbabilities(void) { return probs[activeProbs]; }
        double**                        getTransitionProbabilities(int nodeIdx) { return probs[activeProbs][nodeIdx]; }
        void                            initialize(Model* m, int nn, int ns);
        void                            print(void);
        void                            setNeedsUpdate(bool tf) { needsUpdate = tf; }
        void                            setTransitionProbabilities(void);
    
    private:
                                        TransitionProbabilities(void);
                                       ~TransitionProbabilities(void);
                                        TransitionProbabilities(const TransitionProbabilities& tp) = delete;
        bool                            isInitialized;
        int                             numNodes;
        int                             numStates;
        int                             activeProbs;
        std::vector<double**>           probs[2];
        Model*                          modelPtr;
        bool                            needsUpdate;
};

#endif
