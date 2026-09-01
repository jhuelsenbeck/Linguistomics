#ifndef LikelihoodCalculator_hpp
#define LikelihoodCalculator_hpp

#include <cstddef>
#include <cstdint>
#include <unordered_map>
#include <vector>
#include "MathCache.hpp"
#include "Threads.hpp"
class Alignment;
class Node;
class ParameterAlignment;
class ParameterFrequencies;
class ParameterIndelRates;
class ParameterTree;
class TransitionProbabilities;
class Tree;


struct TKF91Probabilities {

    double                  insertionRate;          // lambda
    double                  deletionRate;           // mu
    double                  immortalProbability;    // prod_{n in T} (1 - B_n), root included
    std::vector<double>     beta;
    std::vector<double>     birthProbability;       // B_n = lambda * beta
    std::vector<double>     extinctionProbability;  // E_n = mu * beta
    std::vector<double>     homologousProbability;  // H_n = e^{-mu t} (1 - lambda beta)
    std::vector<double>     nonHomologousProbability; // N_n = (1 - e^{-mu t} - mu beta)(1 - lambda beta)
};


class LikelihoodCalculator : public LikelihoodTask {

    public:
                                        // finite stand-in for -infinity, which is undefined behavior
                                        // under -ffast-math; small enough to sum over every cognate
        static constexpr double         impossibleLnLikelihood = -1.0e100;

                                        LikelihoodCalculator(void) = delete;
                                        LikelihoodCalculator(TransitionProbabilities* tpc, ParameterAlignment* a, ParameterTree* t, ParameterIndelRates* r, ParameterFrequencies* f);
                                        LikelihoodCalculator(const LikelihoodCalculator&) = delete;
        LikelihoodCalculator&           operator = (const LikelihoodCalculator&) = delete;
                                       ~LikelihoodCalculator(void);

    protected:
        double                          computeLnLikelihood(void) override;

    private:
        struct TKF91Event {

            uint64_t    birthBit;       // single bit for r, tested against the state
            uint64_t    aliveMask;      // nodes labelled other than "-"
            uint64_t    clearMask;      // subtree(r) | disjointLower(r)
            double      weight;
        };

        typedef std::unordered_map<uint64_t,double> StateMap;

        void                            initializeTopology(void);
        void                            cacheTransitionMatrices(void);
        void                            setBirthDeathProbabilities(void);
        void                            setStationaryFrequencies(void);
        void                            buildColumns(void);
        void                            enumerateAlivePatterns(int r, uint64_t emitMask, std::vector<uint64_t>& out);
        double                          patternWeight(int r, uint64_t aliveMask, size_t column);
        void                            buildEventTable(uint64_t emitMask, size_t column, std::vector<TKF91Event>& out);
        void                            applyEvents(const std::vector<TKF91Event>& events, const StateMap& in, StateMap& out);
        void                            silentClosure(StateMap& a);

                                        // sentinel column index for a silent (no emission) event
        static constexpr size_t         silentColumn = static_cast<size_t>(-1);
                                        // a state is a bitmask over post-order positions
        static constexpr size_t         maxNodes = 64;
                                        // silent-event geometric series controls
        static constexpr size_t         maxClosureIterations = 1000;
        static constexpr double         closureTolerance = 1.0e-15;

                                        // pointers and sizes
        Tree*                           tree;
        Alignment*                      alignment;
        double*                         equilibriumFrequencies;
        size_t                          numTaxa;
        size_t                          numNodes;
        size_t                          numStates;
        size_t                          numSegments;

                                        // parameter access
        TransitionProbabilities*        tiProbs;
        ParameterAlignment*             myAlignment;
        ParameterTree*                  myTree;
        ParameterIndelRates*            myIndelRates;
        ParameterFrequencies*           myFrequencies;

        std::vector<Node*>              poNode;
        std::vector<int>                poLeft;          // -1 when absent
        std::vector<int>                poRight;
        std::vector<int>                poSeqRow;        // alignment row for a leaf, -1 otherwise
        std::vector<char>               poIsLeaf;
        std::vector<size_t>             poSubtreeSize;
        std::vector<uint64_t>           poSubtree;       // bitmask of subtree(n)
        std::vector<uint64_t>           poDisjointLower; // disjoint from n and ordered below it
        std::vector<double*>            poTiMatrix;      // transition matrix for the branch above n

                                        // columns of the fixed alignment
        std::vector<uint64_t>           columnEmit;      // emitting leaves, as a post-order bitmask

                                        // scratch, reused across calls
        std::unordered_map<const Node*,int> nodePosition;   // Node -> post-order position
        std::vector<double>             fWork;           // numNodes * numStates
        std::vector<uint64_t>           patternWork;
        std::vector<TKF91Event>         silentEvents;
        std::vector<TKF91Event>         columnEvents;
        StateMap                        alphaCurrent;
        StateMap                        alphaNext;
        StateMap                        closureCurrent;
        StateMap                        closureNext;

        TKF91Probabilities              tkf91Probs;
        unsigned                        taxonMask;
};

#endif
