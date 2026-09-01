#ifndef TKF91StochasticMapping_hpp
#define TKF91StochasticMapping_hpp

#include <cstddef>
#include <cstdint>
#include <string>
#include <unordered_map>
#include <vector>
#include "json.hpp"
#include "Threads.hpp"
#include "UserSettings.hpp"
class Alignment;
class Node;
class ParameterAlignment;
class ParameterFrequencies;
class ParameterIndelRates;
class ParameterTree;
class RandomVariable;
class RateMatrix;
class TransitionProbabilities;
class Tree;


struct MappedEvent {

    enum Kind { Substitution, Insertion, Deletion };

    Kind        kind;
    double      time;                   // 0 = parent end, branchLength = child end
    int         lineage;                // index into MappedHistory::lineages
    int         fromState;              // substitution only, else -1
    int         toState;                // substitution / insertion, else -1
};

struct MappedLineage {

    int                 column;         // alignment column, or -1 if silent
    std::vector<char>   presence;       // presence[nodePostOrderIndex] = 0/1
    std::vector<int>    nodeState;      // state at each alive node, -1 elsewhere
};

struct MappedHistory {

    std::vector<MappedLineage>            lineages;
    std::vector<std::vector<MappedEvent>> events;   // events[branchPostOrderIndex]
    size_t                                numSilentLineages = 0;
};


class TKF91StochasticMapping {

    public:
                                        TKF91StochasticMapping(void) = delete;
                                        TKF91StochasticMapping(TransitionProbabilities* tpc, ParameterAlignment* a, ParameterTree* t, ParameterIndelRates* r, ParameterFrequencies* f, RateMatrix* rateMatrix, SubstitutionModel mt, RandomVariable* rng);
                                        TKF91StochasticMapping(const TKF91StochasticMapping&) = delete;
        TKF91StochasticMapping&         operator = (const TKF91StochasticMapping&) = delete;
        virtual                        ~TKF91StochasticMapping(void);

                                        // MAIN ENTRY: draw one full history from the posterior
        MappedHistory                   sampleHistory(void);
        std::vector<MappedHistory>      sampleHistories(size_t n);

        double                          validateExpectations(size_t n, std::string* report);

        void                            appendHistoryToJsonFile(const MappedHistory& hist, const std::string& path, int cognateIndex, long generation);
        nlohmann::json                  historyToJson(const MappedHistory& hist, int cognateIndex, long generation) const;
        std::string                     newickString(void) const;
        int                             getCognateIndex(void);

    protected:
        virtual void                    drawEndpointConditionedSubstitutions(double segmentLength, int startState, int endState, int lineage, std::vector<MappedEvent>& out);

        int                             drawFreeEndSubstitutions(double segmentLength, int startState, int lineage, std::vector<MappedEvent>& out);
        int                             drawFreshLinkSubstitutions(double segmentLength, int endState, int lineage, std::vector<MappedEvent>& out);

    private:
        struct TKF91Event {

            uint64_t    birthBit;       // bit of the birth node r
            uint64_t    aliveMask;      // nodes the lineage occupies
            uint64_t    clearMask;      // subtree(r) | disjointLower(r)
            double      weight;
            int         column;         // alignment column, or -1 for a silent event
        };

        typedef std::unordered_map<uint64_t,double> StateMap;

        void                            initializeTopology(void);
        void                            setBirthDeathProbabilities(void);
        void                            setStationaryFrequencies(void);
        void                            cacheTransitionMatrices(void);
        void                            cacheRateMatrix(void);
        void                            buildF81RateMatrix(void);
        void                            buildColumns(void);

                                        // event tables: silent (emitMask 0) and one per column
        void                            buildAllEventTables(void);
        void                            enumerateAlivePatterns(int r, uint64_t emitMask, std::vector<uint64_t>& out);
        double                          patternWeight(int r, uint64_t aliveMask, size_t column);
        void                            buildEventTable(uint64_t emitMask, size_t column, std::vector<TKF91Event>& out);
        static uint64_t                 applyEvent(uint64_t state, const TKF91Event& e)
            {
            return (state & ~e.clearMask) | e.aliveMask;
            }

                                        // layer 1: reachable states, backward cost-to-go, forward draw
        void                            enumerateReachableStates(void);
        void                            computeBackward(void);
        void                            sampleEventSequence(std::vector<const TKF91Event*>& drawn);

                                        // layer 2 + 3: turn one lineage's alive-set into node states,
                                        // H/N/E resolution and timed events
        void                            realizeLineage(const TKF91Event& ev, int lineageIndex, MappedHistory& hist);

        int                             simulateForward(double segmentLength, int startState, int lineage, std::vector<MappedEvent>& out);
        void                            emitSubstitutions(std::vector<MappedEvent>& branchEvents, const std::vector<MappedEvent>& subs, double offset) const;
        double                          drawDeathTime(double branchLength);
        void                            requireRateMatrix(void) const;

                                        // rng helpers (default to the injected RandomVariable)
        void                            buildNewick(int i, std::string& s) const;
        double                          uniform01(void);
        double                          exponential(double rate);
        int                             drawFromEquilibrium(void);

        static constexpr double         neverHappens = 1.0e300;   // a waiting time no branch can reach
        static constexpr double         notComputed  = -1.0;      // validateExpectations: not run

        static constexpr size_t         maxNodes = 64;
        static constexpr size_t         maxClosureIterations = 100000;
        static constexpr size_t         maxRejectionAttempts = 1000000;
        static constexpr double         closureTolerance = 1.0e-15;
        static constexpr double         reversibilityTolerance = 1.0e-10;
        static constexpr size_t         silentColumn = static_cast<size_t>(-1);

                                        // injected dependencies
        TransitionProbabilities*        tiProbs;
        ParameterAlignment*             myAlignment;
        ParameterTree*                  myTree;
        ParameterIndelRates*            myIndelRates;
        ParameterFrequencies*           myFrequencies;
        RateMatrix*                     myRateMatrix;
        SubstitutionModel               modelType;
        RandomVariable*                 randomVariable;

                                        // resolved state
        Tree*                           tree;
        Alignment*                      alignment;
        size_t                          numTaxa;
        size_t                          numNodes;
        size_t                          numStates;
        size_t                          numSegments;

                                        // TKF91 per-node factors, post-order indexed
        double                          insertionRate;
        double                          deletionRate;
        double                          immortalProbability;
        std::vector<double>             beta;
        std::vector<double>             birthProbability;
        std::vector<double>             extinctionProbability;
        std::vector<double>             homologousProbability;
        std::vector<double>             nonHomologousProbability;
        std::vector<double>             equilibriumFrequencies;
        std::vector<double>             rateMatrix;        // Q, row-major from x to

                                        // topology in post-order coordinates
        std::vector<Node*>              poNode;
        std::vector<int>                poLeft;
        std::vector<int>                poRight;
        std::vector<int>                poParent;
        std::vector<int>                poSeqRow;
        std::vector<char>               poIsLeaf;
        std::vector<size_t>             poSubtreeSize;
        std::vector<uint64_t>           poSubtree;
        std::vector<uint64_t>           poDisjointLower;
        std::vector<double>             poBranchLength;
        std::vector<double*>            poTiMatrix;
        std::unordered_map<const Node*,int> nodePosition;

                                        // columns and events
        std::vector<uint64_t>           columnEmit;
        std::vector<TKF91Event>         silentEvents;
        std::vector<std::vector<TKF91Event>> columnEventTables;

                                        // layer 1 working set
        std::vector<uint64_t>           states;
        std::unordered_map<uint64_t,int> stateIndex;
        std::vector<std::vector<double>> backward;         // backward[c][stateIdx], c = 0..numSegments+1

                                        // scratch
        std::vector<double>             fWork;
        std::vector<uint64_t>           patternWork;
        std::vector<MappedEvent>        substitutionWork;
        std::vector<MappedEvent>        reversalWork;

        bool                            rateMatrixIsReversible;

        unsigned                        taxonMask;
};

#endif
