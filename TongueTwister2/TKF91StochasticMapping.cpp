#include <algorithm>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <system_error>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <limits>
#include "Alignment.hpp"
#include "Container.hpp"
#include "Msg.hpp"
#include "Node.hpp"
#include "ParameterAlignment.hpp"
#include "ParameterFrequencies.hpp"
#include "ParameterIndelRates.hpp"
#include "ParameterTree.hpp"
#include "RandomVariable.hpp"
#include "RateMatrix.hpp"
#include "TKF91StochasticMapping.hpp"
#include "TransitionProbabilitiesGpu.hpp"
#include "Tree.hpp"


// isnan and !(x > 0.0) can both be folded away under -ffast-math, so look at the bits
static inline bool isNotFinite(double x) {

    uint64_t bits;
    std::memcpy(&bits, &x, sizeof(bits));
    return ((bits >> 52) & 0x7FFULL) == 0x7FFULL;      // all-ones exponent: Inf or NaN
}

static inline bool isPositiveFinite(double x) {

    return (isNotFinite(x) == false) && (x > 0.0);
}


TKF91StochasticMapping::TKF91StochasticMapping(TransitionProbabilities* tpc, ParameterAlignment* a, ParameterTree* t, ParameterIndelRates* r, ParameterFrequencies* f, RateMatrix* rateMatrix, SubstitutionModel mt, RandomVariable* rng) :
    tiProbs(tpc), myAlignment(a), myTree(t), myIndelRates(r), myFrequencies(f),
    myRateMatrix(rateMatrix), modelType(mt), randomVariable(rng) {

    alignment = myAlignment->getAlignment();
    if (alignment == nullptr)
        Msg::error("Null alignment when initializing TKF91 stochastic mapping");
    taxonMask = myAlignment->getTaxonMask();
    numStates = myAlignment->getNumStates();
    if (numStates == 0)
        Msg::error("Expecting at least one character state in TKF91 stochastic mapping");
    equilibriumFrequencies.assign(numStates, 0.0);

    tree = myTree->getTree(taxonMask);
    if (tree == nullptr)
        Msg::error("Could not find tree when initializing TKF91 stochastic mapping");
    numSegments = alignment->getNumSegments();
    initializeTopology();
}

TKF91StochasticMapping::~TKF91StochasticMapping(void) {
}


// --------------------------------------------------------------- rng helpers

double TKF91StochasticMapping::uniform01(void) {

    if (randomVariable == nullptr)
        Msg::error("TKF91 stochastic mapping has no RandomVariable");
    return randomVariable->uniformRv();
}

double TKF91StochasticMapping::exponential(double rate) {

    if (isPositiveFinite(rate) == false)
        return neverHappens;
    return -std::log(1.0 - uniform01()) / rate;
}

int TKF91StochasticMapping::drawFromEquilibrium(void) {

    double u = uniform01(), acc = 0.0;
    for (size_t s=0; s<numStates; s++)
        {
        acc += equilibriumFrequencies[s];
        if (u <= acc)
            return (int)s;
        }
    return (int)numStates - 1;
}

int TKF91StochasticMapping::getCognateIndex(void) {

    return (int)myAlignment->getCognateIndex();
}

// ----------------------------------------------------------- shared setup
// (topology / factors / columns identical in meaning to LikelihoodCalculator)

void TKF91StochasticMapping::initializeTopology(void) {

    const std::vector<Node*>& postOrder = tree->getPostOrder();
    numNodes = postOrder.size();
    if (numNodes == 0)
        Msg::error("Empty post-order traversal in TKF91 stochastic mapping");
    if (numNodes > maxNodes)
        Msg::error("TKF91 stochastic mapping handles at most " + std::to_string(maxNodes) + " nodes");

    poNode.assign(numNodes, nullptr);
    poLeft.assign(numNodes, -1);
    poRight.assign(numNodes, -1);
    poParent.assign(numNodes, -1);
    poSeqRow.assign(numNodes, -1);
    poIsLeaf.assign(numNodes, 0);
    poSubtreeSize.assign(numNodes, 0);
    poSubtree.assign(numNodes, 0);
    poDisjointLower.assign(numNodes, 0);
    poBranchLength.assign(numNodes, 0.0);
    poTiMatrix.assign(numNodes, nullptr);

    std::unordered_map<const Node*,int>& position = nodePosition;
    position.clear();
    position.reserve(numNodes * 2);
    for (size_t i=0; i<numNodes; i++)
        {
        if (postOrder[i] == nullptr)
            Msg::error("Null node in the post-order traversal of TKF91 stochastic mapping");
        poNode[i] = postOrder[i];
        position[postOrder[i]] = (int)i;
        }

    size_t numLeaves = 0;
    for (size_t i=0; i<numNodes; i++)
        {
        Node* p = poNode[i];
        poIsLeaf[i] = (p->getIsLeaf() == true) ? 1 : 0;
        poBranchLength[i] = p->getBranchLength();
        if (p->getAncestor() != nullptr)
            poParent[i] = position[p->getAncestor()];
        if (poIsLeaf[i])
            {
            numLeaves++;
            poSeqRow[i] = p->getIndex();
            }
        else
            {
            Node* lft = p->getDescendant(0);
            Node* rht = p->getDescendant(1);
            if (lft == nullptr || rht == nullptr)
                Msg::error("Interior node without two descendants in TKF91 stochastic mapping");
            poLeft[i] = position[lft];
            poRight[i] = position[rht];
            }
        uint64_t m = ((uint64_t)1 << i);
        size_t sz = 1;
        if (poLeft[i] >= 0)
            {
            m |= poSubtree[poLeft[i]];
            sz += poSubtreeSize[poLeft[i]];
            }
        if (poRight[i] >= 0)
            {
            m |= poSubtree[poRight[i]];
            sz += poSubtreeSize[poRight[i]];
            }
        poSubtree[i] = m;
        poSubtreeSize[i] = sz;
        }

    numTaxa = numLeaves;
    if (alignment->getNumTaxa() != numTaxa)
        Msg::error("Tree has " + std::to_string(numTaxa) + " leaves but the alignment has " +
                   std::to_string(alignment->getNumTaxa()) + " rows in TKF91 stochastic mapping");

    for (size_t r=0; r<numNodes; r++)
        {
        uint64_t m = 0;
        for (size_t k=0; k<r; k++)
            {
            const bool kBelowR = ((poSubtree[r] >> k) & 1) != 0;
            const bool rBelowK = ((poSubtree[k] >> r) & 1) != 0;
            if (kBelowR == false && rBelowK == false)
                m |= ((uint64_t)1 << k);
            }
        poDisjointLower[r] = m;
        }

    beta.resize(numNodes);
    birthProbability.resize(numNodes);
    extinctionProbability.resize(numNodes);
    homologousProbability.resize(numNodes);
    nonHomologousProbability.resize(numNodes);
    fWork.resize(numNodes * numStates);
}

void TKF91StochasticMapping::setBirthDeathProbabilities(void) {

    insertionRate = myIndelRates->getInsertionRate();
    deletionRate  = myIndelRates->getDeletionRate();
    const double lambda = insertionRate, mu = deletionRate;
    if (isPositiveFinite(lambda) == false || isNotFinite(mu) || !(mu > lambda))
        Msg::error("TKF91 requires 0 < lambda < mu");

    immortalProbability = 1.0;
    for (size_t i=0; i<numNodes; i++)
        {
        Node* p = poNode[i];
        const double v = p->getBranchLength();
        if (p->getAncestor() == nullptr)
            {
            beta[i] = 1.0 / mu;
            homologousProbability[i] = 0.0;
            }
        else
            {
            const double expPart = std::exp((lambda - mu) * v);
            beta[i] = (1.0 - expPart) / (mu - lambda * expPart);
            homologousProbability[i] = std::exp(-mu * v) * (1.0 - lambda * beta[i]);
            }
        birthProbability[i]      = lambda * beta[i];
        extinctionProbability[i] = mu * beta[i];
        nonHomologousProbability[i] = (1.0 - mu * beta[i]) * (1.0 - birthProbability[i]) - homologousProbability[i];
        immortalProbability *= (1.0 - birthProbability[i]);
        }
}

void TKF91StochasticMapping::setStationaryFrequencies(void) {

    if (myFrequencies == nullptr)
        {
        const double x = 1.0 / numStates;
        for (size_t i=0; i<numStates; i++)
            equilibriumFrequencies[i] = x;
        }
    else
        {
        std::vector<double>& x = myFrequencies->getFrequencies();
        if (x.size() != numStates)
            Msg::error("Stationary frequency vector has the wrong dimension in TKF91 stochastic mapping");
        for (size_t i=0; i<numStates; i++)
            equilibriumFrequencies[i] = x[i];
        }
}

void TKF91StochasticMapping::cacheTransitionMatrices(void) {

    for (size_t i=0; i<numNodes; i++)
        {
        Node* p = poNode[i];
        poTiMatrix[i] = (p->getAncestor() != nullptr)
                      ? tiProbs->getTransitionProbability(p->getBranchLength()).begin()
                      : nullptr;
        }
}

void TKF91StochasticMapping::buildF81RateMatrix(void) {

    double sumSquared = 0.0;
    for (size_t i=0; i<numStates; i++)
        sumSquared += equilibriumFrequencies[i] * equilibriumFrequencies[i];

    const double denominator = 1.0 - sumSquared;
    if (isPositiveFinite(denominator) == false)
        Msg::error("Degenerate stationary frequencies in TKF91 stochastic mapping: cannot normalize the rate matrix");
    const double c = 1.0 / denominator;

    rateMatrix.assign(numStates * numStates, 0.0);
    for (size_t i=0; i<numStates; i++)
        {
        double rowSum = 0.0;
        for (size_t j=0; j<numStates; j++)
            {
            if (i == j)
                continue;
            rateMatrix[i * numStates + j] = c * equilibriumFrequencies[j];
            rowSum += rateMatrix[i * numStates + j];
            }
        rateMatrix[i * numStates + i] = -rowSum;
        }

    rateMatrixIsReversible = true;      // pi(i) c pi(j) == pi(j) c pi(i), identically
}

void TKF91StochasticMapping::cacheRateMatrix(void) {

    if (myRateMatrix == nullptr)
        {
        if (modelType == jc69 || modelType == f81)
            {
            buildF81RateMatrix();
            return;
            }
        Msg::error("TKF91 stochastic mapping was constructed without a rate matrix, so substitution histories cannot be recorded");
        }
    DoubleMatrix& Q = myRateMatrix->getQ();
    if (Q.getNumRows() != numStates || Q.getNumCols() != numStates)
        Msg::error("Rate matrix dimension does not match the number of states in TKF91 stochastic mapping");
    rateMatrix.assign(numStates * numStates, 0.0);
    for (size_t i=0; i<numStates; i++)
        for (size_t j=0; j<numStates; j++)
            rateMatrix[i * numStates + j] = Q(i, j);

    rateMatrixIsReversible = true;
    for (size_t i=0; i<numStates && rateMatrixIsReversible == true; i++)
        {
        for (size_t k=i+1; k<numStates; k++)
            {
            const double forward = equilibriumFrequencies[i] * rateMatrix[i * numStates + k];
            const double reverse = equilibriumFrequencies[k] * rateMatrix[k * numStates + i];
            const double scale = std::fabs(forward) + std::fabs(reverse);
            if (std::fabs(forward - reverse) > reversibilityTolerance * ((scale > 0.0) ? scale : 1.0))
                {
                rateMatrixIsReversible = false;
                break;
                }
            }
        }
}

void TKF91StochasticMapping::buildColumns(void) {

    columnEmit.assign(numSegments, 0);
    for (size_t segment=0; segment<numSegments; segment++)
        {
        uint64_t emit = 0;
        for (size_t i=0; i<numNodes; i++)
            {
            if (poIsLeaf[i] == 0)
                continue;
            const size_t charState = (size_t)(*alignment)(poSeqRow[i], segment);
            if (charState == numStates)
                continue;
            if (charState > numStates)
                Msg::error("Unrepresentable character state in TKF91 stochastic mapping");
            emit |= ((uint64_t)1 << i);
            }
        if (emit == 0)
            Msg::error("All-gap column " + std::to_string(segment) + " in TKF91 stochastic mapping");
        columnEmit[segment] = emit;
        }
}


// ------------------------------------------ alive-pattern / event machinery

void TKF91StochasticMapping::enumerateAlivePatterns(int r, uint64_t emitMask, std::vector<uint64_t>& out) {

    out.clear();
    if (poIsLeaf[r])
        {
        if ((emitMask >> r) & 1)
            out.push_back((uint64_t)1 << r);
        return;
        }
    std::vector<uint64_t> lftOptions, rhtOptions;
    for (int pass=0; pass<2; pass++)
        {
        const int c = (pass == 0) ? poLeft[r] : poRight[r];
        std::vector<uint64_t>& options = (pass == 0) ? lftOptions : rhtOptions;
        options.clear();
        const uint64_t needed = emitMask & poSubtree[c];
        if (needed == 0)
            options.push_back(0);
        if (poIsLeaf[c])
            {
            if (needed != 0)
                options.push_back((uint64_t)1 << c);
            }
        else
            {
            std::vector<uint64_t> sub;
            enumerateAlivePatterns(c, emitMask, sub);
            for (uint64_t s : sub)
                options.push_back(s);
            }
        }
    const uint64_t rBit = ((uint64_t)1 << r);
    for (uint64_t l : lftOptions)
        for (uint64_t h : rhtOptions)
            out.push_back(rBit | l | h);
}

double TKF91StochasticMapping::patternWeight(int r, uint64_t aliveMask, size_t column) {

    const size_t n = numStates;
    const size_t lo = (size_t)r + 1 - poSubtreeSize[r];
    for (size_t i=lo; i<=(size_t)r; i++)
        {
        if (((aliveMask >> i) & 1) == 0)
            continue;
        double* const fi = fWork.data() + i * n;
        if (poIsLeaf[i])
            {
            std::memset(fi, 0, n * sizeof(double));
            const size_t charState = (size_t)(*alignment)(poSeqRow[i], column);
            fi[charState] = 1.0;
            continue;
            }
        for (size_t a=0; a<n; a++)
            fi[a] = 1.0;
        for (int pass=0; pass<2; pass++)
            {
            const int c = (pass == 0) ? poLeft[i] : poRight[i];
            if (((aliveMask >> c) & 1) == 0)
                {
                const double e = extinctionProbability[c];
                for (size_t a=0; a<n; a++)
                    fi[a] *= e;
                continue;
                }
            const double* const fc = fWork.data() + (size_t)c * n;
            const double* const tp = poTiMatrix[c];
            const double hc = homologousProbability[c];
            const double nc = nonHomologousProbability[c];
            double piDotF = 0.0;
            for (size_t g=0; g<n; g++)
                piDotF += equilibriumFrequencies[g] * fc[g];
            const double newTerm = nc * piDotF;
            for (size_t a=0; a<n; a++)
                {
                const double* const row = tp + a * n;
                double dot = 0.0;
                for (size_t g=0; g<n; g++)
                    dot += row[g] * fc[g];
                fi[a] *= (hc * dot + newTerm);
                }
            }
        }
    const double* const fr = fWork.data() + (size_t)r * n;
    double total = 0.0;
    for (size_t a=0; a<n; a++)
        total += equilibriumFrequencies[a] * fr[a];
    return birthProbability[r] * total;
}

void TKF91StochasticMapping::buildEventTable(uint64_t emitMask, size_t column, std::vector<TKF91Event>& out) {

    out.clear();
    for (size_t r=0; r<numNodes; r++)
        {
        if ((emitMask & ~poSubtree[r]) != 0)
            continue;
        enumerateAlivePatterns((int)r, emitMask, patternWork);
        for (uint64_t alive : patternWork)
            {
            const double w = patternWeight((int)r, alive, column);
            if (w == 0.0)
                continue;
            TKF91Event e;
            e.birthBit  = ((uint64_t)1 << r);
            e.aliveMask = alive;
            e.clearMask = poSubtree[r] | poDisjointLower[r];
            e.weight    = w;
            e.column    = (column == silentColumn) ? -1 : (int)column;
            out.push_back(e);
            }
        }
}

void TKF91StochasticMapping::buildAllEventTables(void) {

    // silent events do not depend on a column; column c uses site c for its leaf
    // characters (its aliveMask leaves are the emitting leaves of column c)
    buildEventTable(0, silentColumn, silentEvents);
    columnEventTables.assign(numSegments, std::vector<TKF91Event>());
    for (size_t c=0; c<numSegments; c++)
        buildEventTable(columnEmit[c], c, columnEventTables[c]);
}


// ------------------------------------------------- layer 1: skeleton sampler

void TKF91StochasticMapping::enumerateReachableStates(void) {

    // BFS from the all-alive state under silent events and every column emission
    states.clear();
    stateIndex.clear();
    const uint64_t full = (numNodes == 64) ? ~(uint64_t)0 : (((uint64_t)1 << numNodes) - 1);

    std::vector<uint64_t> frontier;
    auto push = [&](uint64_t s)
        {
        if (stateIndex.find(s) == stateIndex.end())
            {
            stateIndex[s] = (int)states.size();
            states.push_back(s);
            frontier.push_back(s);
            }
        };
    push(full);

    while (!frontier.empty())
        {
        std::vector<uint64_t> next;
        for (uint64_t s : frontier)
            {
            for (const TKF91Event& e : silentEvents)
                {
                if (s & e.birthBit)
                    {
                    uint64_t s2 = applyEvent(s, e);
                    if (stateIndex.find(s2) == stateIndex.end())
                        {
                        stateIndex[s2] = (int)states.size();
                        states.push_back(s2);
                        next.push_back(s2);
                        }
                    }
                }
            for (size_t c=0; c<numSegments; c++)
                {
                for (const TKF91Event& e : columnEventTables[c])
                    {
                    if (s & e.birthBit)
                        {
                        uint64_t s2 = applyEvent(s, e);
                        if (stateIndex.find(s2) == stateIndex.end())
                            {
                            stateIndex[s2] = (int)states.size();
                            states.push_back(s2);
                            next.push_back(s2);
                            }
                        }
                    }
                }
            }
        frontier.swap(next);
        }
}

void TKF91StochasticMapping::computeBackward(void) {

    // B[c][S] for c = 1..numSegments+1 (phase 1..C+1), plus B[0] unused.
    // B[C+1] = (I-Z)^{-1} 1 ; B[c] = (I-Z)^{-1} (E_c B[c+1]).
    const size_t C = numSegments;
    const size_t ns = states.size();
    backward.assign(C + 2, std::vector<double>(ns, 0.0));

    auto closure = [&](std::vector<double>& base, std::vector<double>& out) {
        // solve x = base + Z x  by Neumann iteration (cost-to-go form)
        out = base;
        for (size_t it=0; it<maxClosureIterations; it++)
            {
            double maxDelta = 0.0, scale = 1.0;
            std::vector<double> nx = base;
            for (size_t si=0; si<ns; si++)
                {
                const uint64_t s = states[si];
                double add = 0.0;
                for (const TKF91Event& e : silentEvents)
                    if (s & e.birthBit)
                        add += e.weight * out[stateIndex[applyEvent(s, e)]];
                nx[si] += add;
                }
            for (size_t si=0; si<ns; si++)
                {
                maxDelta = std::max(maxDelta, std::fabs(nx[si] - out[si]));
                scale = std::max(scale, std::fabs(nx[si]));
                }
            out.swap(nx);
            if (maxDelta <= closureTolerance * scale)
                return;
            }
        Msg::error("Backward silent-event series failed to converge in TKF91 stochastic mapping");
    };

    // phase C+1: base = 1
    {
        std::vector<double> base(ns, 1.0);
        closure(base, backward[C + 1]);
    }
    for (size_t c=C; c>=1; c--)
        {
        std::vector<double> base(ns, 0.0);
        for (size_t si=0; si<ns; si++)
            {
            const uint64_t s = states[si];
            double v = 0.0;
            for (const TKF91Event& e : columnEventTables[c - 1])
                if (s & e.birthBit)
                    v += e.weight * backward[c + 1][stateIndex[applyEvent(s, e)]];
            base[si] = v;
            }
        closure(base, backward[c]);
        if (c == 1)
            break;
        }
}

void TKF91StochasticMapping::sampleEventSequence(std::vector<const TKF91Event*>& drawn) {

    drawn.clear();
    const size_t C = numSegments;
    const uint64_t full = (numNodes == 64) ? ~(uint64_t)0 : (((uint64_t)1 << numNodes) - 1);
    uint64_t S = full;
    size_t c = 1;

    std::vector<const TKF91Event*> cand;
    std::vector<int> candPhase;         // phase after the event (for column advance)
    std::vector<double> w;

    for (size_t guard=0; ; guard++)
        {
        if (guard > 10ull * (C + 1) + 1000000ull)
            Msg::error("Runaway history in TKF91 stochastic mapping");

        cand.clear();
        candPhase.clear();
        w.clear();

        // silent events (phase unchanged)
        for (const TKF91Event& e : silentEvents)
            if (S & e.birthBit)
                {
                const uint64_t s2 = applyEvent(S, e);
                cand.push_back(&e);
                candPhase.push_back((int)c);
                w.push_back(e.weight * backward[c][stateIndex[s2]]);
                }
        // emit column c (phase -> c+1)
        if (c <= C)
            for (const TKF91Event& e : columnEventTables[c - 1])
                if (S & e.birthBit)
                    {
                    const uint64_t s2 = applyEvent(S, e);
                    cand.push_back(&e);
                    candPhase.push_back((int)c + 1);
                    w.push_back(e.weight * backward[c + 1][stateIndex[s2]]);
                    }
        // terminate (only in phase C+1)
        bool canTerminate = (c == C + 1);
        double terminateW = canTerminate ? 1.0 : 0.0;

        double norm = terminateW;
        for (double x : w) norm += x;
        if (isPositiveFinite(norm) == false)
            Msg::error("Zero-mass state in TKF91 stochastic mapping forward sampler");

        double u = uniform01() * norm, acc = 0.0;
        int chosen = -1;
        for (size_t k=0; k<w.size(); k++)
            {
            acc += w[k];
            if (u <= acc)
                {
                chosen = (int)k;
                break;
                }
            }
        if (chosen < 0)
            {
            // terminate
            if (!canTerminate)
                Msg::error("Sampling fell through without terminating in TKF91 stochastic mapping");
            return;
            }

        drawn.push_back(cand[chosen]);
        S = applyEvent(S, *cand[chosen]);
        c = (size_t)candPhase[chosen];
        }
}


// ------------------------------------ layers 2 + 3: realize one link lineage

void TKF91StochasticMapping::realizeLineage(const TKF91Event& ev, int lineageIndex, MappedHistory& hist) {

    const size_t n = numStates;
    const uint64_t A = ev.aliveMask;
    const int r = __builtin_ctzll(ev.birthBit);   // birth node = lowest set bit of birthBit

    MappedLineage lin;
    lin.column = ev.column;
    lin.presence.assign(numNodes, 0);
    lin.nodeState.assign(numNodes, -1);
    for (size_t i=0; i<numNodes; i++)
        if ((A >> i) & 1)
            lin.presence[i] = 1;

    // ---- layer 2a: peel this lineage (fill fWork for its alive nodes) ----
    // recompute the same f used in patternWeight so we can FFBS the node states
    const size_t lo = (size_t)r + 1 - poSubtreeSize[r];
    for (size_t i=lo; i<=(size_t)r; i++)
        {
        if (((A >> i) & 1) == 0)
            continue;
        double* const fi = fWork.data() + i * n;
        if (poIsLeaf[i])
            {
            std::memset(fi, 0, n * sizeof(double));
            const int col = (ev.column >= 0) ? ev.column : 0;      // silent lineages have no leaves
            if (poIsLeaf[i] && ((A >> i) & 1))
                {
                const size_t cs = (size_t)(*alignment)(poSeqRow[i], col);
                fi[cs] = 1.0;
                }
            continue;
            }
        for (size_t a=0; a<n; a++) fi[a] = 1.0;
        for (int pass=0; pass<2; pass++)
            {
            const int c = (pass == 0) ? poLeft[i] : poRight[i];
            if (((A >> c) & 1) == 0)
                {
                const double e = extinctionProbability[c];
                for (size_t a=0;a<n;a++) fi[a]*=e;
                continue;
                }
            const double* const fc = fWork.data() + (size_t)c * n;
            const double* const tp = poTiMatrix[c];
            const double hc = homologousProbability[c];
            const double nc = nonHomologousProbability[c];
            double piDotF = 0.0;
            for (size_t g=0; g<n; g++) piDotF += equilibriumFrequencies[g]*fc[g];
            const double newTerm = nc * piDotF;
            for (size_t a=0; a<n; a++)
                {
                const double* const row = tp + a*n;
                double dot = 0.0;
                for (size_t g=0; g<n; g++) dot += row[g]*fc[g];
                fi[a] *= (hc*dot + newTerm);
                }
            }
        }

    // ---- layer 2b: FFBS node states, top (r) down ----
    // birth-node state ~ pi(a) f_r(a)
    {
        const double* const fr = fWork.data() + (size_t)r * n;
        double norm = 0.0;
        for (size_t a=0; a<n; a++) norm += equilibriumFrequencies[a]*fr[a];
        double u = uniform01()*norm, acc = 0.0;
        int drawn = (int)n-1;
        for (size_t a=0; a<n; a++){ acc += equilibriumFrequencies[a]*fr[a]; if (u<=acc){drawn=(int)a;break;} }
        lin.nodeState[r] = drawn;
    }
    // children in reverse post-order (parent before child)
    for (size_t ii=0; ii<numNodes; ii++)
        {
        const size_t i = numNodes - 1 - ii;
        if (((A >> i) & 1) == 0 || (int)i == r)
            continue;
        const int par = poParent[i];
        if (par < 0 || ((A >> par) & 1) == 0 || lin.nodeState[par] < 0)
            continue;
        if (poIsLeaf[i])
            {
            lin.nodeState[i] = (ev.column >= 0) ? (int)(*alignment)(poSeqRow[i], ev.column) : drawFromEquilibrium();
            continue;
            }

        // child state g ~ (H_i p(a->g) + N_i pi(g)) f_i(g)
        const int a = lin.nodeState[par];
        const double* const fi = fWork.data() + i*n;
        const double* const tp = poTiMatrix[i];
        const double hi = homologousProbability[i], ni = nonHomologousProbability[i];
        const double* const row = tp + (size_t)a*n;
        double norm = 0.0;
        for (size_t g=0; g<n; g++) norm += (hi*row[g] + ni*equilibriumFrequencies[g]) * fi[g];
        double u = uniform01()*norm, acc = 0.0;
        int drawn = (int)n-1;
        for (size_t g=0; g<n; g++){ acc += (hi*row[g]+ni*equilibriumFrequencies[g])*fi[g]; if (u<=acc){drawn=(int)g;break;} }
        lin.nodeState[i] = drawn;
        }

    if (poParent[r] >= 0)
        {
        const int born = drawFreshLinkSubstitutions(poBranchLength[r], lin.nodeState[r], lineageIndex, substitutionWork);
        MappedEvent b;
        b.kind      = MappedEvent::Insertion;
        b.time      = 0.0;
        b.lineage   = lineageIndex;
        b.fromState = -1;
        b.toState   = born;
        hist.events[r].push_back(b);
        emitSubstitutions(hist.events[r], substitutionWork, 0.0);
        }

    for (size_t v=0; v<numNodes; v++)
        {
        const int u = poParent[v];
        if (u < 0)
            continue;
        const bool aliveU = ((A >> u) & 1) != 0;
        const bool aliveV = ((A >> v) & 1) != 0;
        const double L = poBranchLength[v];

        if (aliveU && aliveV)
            {
            // split H vs N given the two endpoint states, exactly the two terms
            // of the alive-child factor
            const int a = lin.nodeState[u], g = lin.nodeState[v];
            const double* const tp = poTiMatrix[v];
            const double hTerm = homologousProbability[v] * tp[(size_t)a*n + (size_t)g];
            const double nTerm = nonHomologousProbability[v] * equilibriumFrequencies[g];
            const bool survived = (uniform01() * (hTerm + nTerm) <= hTerm);
            if (survived)
                {
                // link persisted: substitutions a->g along the branch
                drawEndpointConditionedSubstitutions(L, a, g, lineageIndex, substitutionWork);
                emitSubstitutions(hist.events[v], substitutionWork, 0.0);
                }
            else
                {
                const double tDeath = drawDeathTime(L);
                const int died = drawFreeEndSubstitutions(tDeath, a, lineageIndex, substitutionWork);
                emitSubstitutions(hist.events[v], substitutionWork, 0.0);
                MappedEvent d;
                d.kind      = MappedEvent::Deletion;
                d.time      = tDeath;
                d.lineage   = lineageIndex;
                d.fromState = died;
                d.toState   = -1;
                hist.events[v].push_back(d);

                const double tBirth = tDeath + uniform01()*(L - tDeath);
                const int born = drawFreshLinkSubstitutions(L - tBirth, g, lineageIndex, substitutionWork);
                MappedEvent b;
                b.kind      = MappedEvent::Insertion;
                b.time      = tBirth;
                b.lineage   = lineageIndex;
                b.fromState = -1;
                b.toState   = born;
                hist.events[v].push_back(b);
                emitSubstitutions(hist.events[v], substitutionWork, tBirth);
                }
            }
        else if (aliveU && !aliveV)
            {
            const double tDeath = drawDeathTime(L);
            const int died = drawFreeEndSubstitutions(tDeath, lin.nodeState[u], lineageIndex, substitutionWork);
            emitSubstitutions(hist.events[v], substitutionWork, 0.0);
            MappedEvent d;
            d.kind      = MappedEvent::Deletion;
            d.time      = tDeath;
            d.lineage   = lineageIndex;
            d.fromState = died;
            d.toState   = -1;
            hist.events[v].push_back(d);
            }
        }

    hist.lineages[lineageIndex] = std::move(lin);
}


// ------------------------------------------------------------- default layer 3

void TKF91StochasticMapping::requireRateMatrix(void) const {

    if (rateMatrix.empty())
        Msg::error("No rate matrix is cached in TKF91 stochastic mapping, so substitutions cannot be drawn");
}

double TKF91StochasticMapping::drawDeathTime(double branchLength) {

    // time of a deletion known to fall within the branch: exponential with rate mu,
    // truncated to [0, branchLength]. Degenerates to uniform as mu*L goes to zero.
    if (deletionRate * branchLength < 1e-12)
        return uniform01() * branchLength;
    return -std::log(1.0 - uniform01() * (1.0 - std::exp(-deletionRate * branchLength))) / deletionRate;
}

void TKF91StochasticMapping::emitSubstitutions(std::vector<MappedEvent>& branchEvents, const std::vector<MappedEvent>& subs, double offset) const {

    // the samplers work in segment-local time; branch time is what gets recorded
    branchEvents.reserve(branchEvents.size() + subs.size());
    for (size_t i=0; i<subs.size(); i++)
        {
        MappedEvent e = subs[i];
        e.time += offset;
        branchEvents.push_back(e);
        }
}

int TKF91StochasticMapping::simulateForward(double segmentLength, int startState, int lineage, std::vector<MappedEvent>& out) {

    out.clear();
    const size_t n = numStates;
    int cur = startState;
    double t = 0.0;
    for (;;)
        {
        const double qii = -rateMatrix[(size_t)cur*n + (size_t)cur];
        if (isPositiveFinite(qii) == false)
            break;                                    // absorbing state: nothing further happens
        t += exponential(qii);
        if (t >= segmentLength)
            break;
        double u = uniform01()*qii, acc = 0.0;
        int nxt = cur;
        for (size_t j=0; j<n; j++)
            {
            if ((int)j == cur)
                continue;
            acc += rateMatrix[(size_t)cur*n + j];
            if (u <= acc)
                {
                nxt = (int)j;
                break;
                }
            }
        MappedEvent e;
        e.kind      = MappedEvent::Substitution;
        e.time      = t;
        e.lineage   = lineage;
        e.fromState = cur;
        e.toState   = nxt;
        out.push_back(e);
        cur = nxt;
        }
    return cur;
}

/* Both ends pinned: the link is present at both nodes and its states there are
   already drawn. Nielsen rejection sampling. */
void TKF91StochasticMapping::drawEndpointConditionedSubstitutions(double segmentLength, int startState, int endState, int lineage, std::vector<MappedEvent>& out) {

    requireRateMatrix();
    for (size_t attempt=0; attempt<maxRejectionAttempts; attempt++)
        {
        if (simulateForward(segmentLength, startState, lineage, out) == endState)
            return;
        }
    Msg::error("Endpoint-conditioned substitution sampling failed to converge in TKF91 stochastic mapping");
}

int TKF91StochasticMapping::drawFreeEndSubstitutions(double segmentLength, int startState, int lineage, std::vector<MappedEvent>& out) {

    requireRateMatrix();
    return simulateForward(segmentLength, startState, lineage, out);
}

int TKF91StochasticMapping::drawFreshLinkSubstitutions(double segmentLength, int endState, int lineage, std::vector<MappedEvent>& out) {

    requireRateMatrix();

    if (rateMatrixIsReversible == true)
        {
        const int born = simulateForward(segmentLength, endState, lineage, reversalWork);
        out.clear();
        out.reserve(reversalWork.size());
        for (size_t i=reversalWork.size(); i>0; i--)
            {
            const MappedEvent& f = reversalWork[i-1];
            MappedEvent e;
            e.kind      = MappedEvent::Substitution;
            e.time      = segmentLength - f.time;     // later in the forward pass is earlier here
            e.lineage   = lineage;
            e.fromState = f.toState;                  // and the jump runs the other way
            e.toState   = f.fromState;
            out.push_back(e);
            }
        return born;
        }

    for (size_t attempt=0; attempt<maxRejectionAttempts; attempt++)
        {
        const int b = drawFromEquilibrium();
        if (simulateForward(segmentLength, b, lineage, out) == endState)
            return b;
        }
    Msg::error("Fresh-link substitution sampling failed to converge in TKF91 stochastic mapping");
    return endState;
}


// ------------------------------------------------------------------ entry point

MappedHistory TKF91StochasticMapping::sampleHistory(void) {

    tree = myTree->getTree(taxonMask);
    if (tree == nullptr)
        Msg::error("Could not find tree in TKF91 stochastic mapping");
    alignment = myAlignment->getAlignment();
    numSegments = alignment->getNumSegments();

    initializeTopology();
    setBirthDeathProbabilities();
    setStationaryFrequencies();
    cacheTransitionMatrices();
    cacheRateMatrix();
    buildColumns();
    buildAllEventTables();

    enumerateReachableStates();
    computeBackward();

    if (const char* probe = std::getenv("TKF91_PROBE_NORMALIZER"))
        {
        (void)probe;
        const uint64_t full = (numNodes==64)?~(uint64_t)0:(((uint64_t)1<<numNodes)-1);
        double Z = immortalProbability * backward[1][stateIndex[full]];
        std::fprintf(stderr, "NORMALIZER %.15e\n", Z);
        }

    std::vector<const TKF91Event*> drawn;
    sampleEventSequence(drawn);

    MappedHistory hist;
    hist.events.assign(numNodes, std::vector<MappedEvent>());
    hist.lineages.assign(drawn.size(), MappedLineage());
    hist.numSilentLineages = 0;

    for (size_t k=0; k<drawn.size(); k++)
        {
        realizeLineage(*drawn[k], (int)k, hist);
        if (drawn[k]->column < 0)
            hist.numSilentLineages++;
        }

    for (size_t v=0; v<numNodes; v++)
        std::sort(hist.events[v].begin(), hist.events[v].end(),
                  [](const MappedEvent& x, const MappedEvent& y){ return x.time < y.time; });

    return hist;
}

std::vector<MappedHistory> TKF91StochasticMapping::sampleHistories(size_t n) {

    std::vector<MappedHistory> out;
    out.reserve(n);
    for (size_t i=0; i<n; i++)
        out.push_back(sampleHistory());
    return out;
}

double TKF91StochasticMapping::validateExpectations(size_t n, std::string* report) {

    (void)n;
    if (report != nullptr)
        *report = "validateExpectations() is not yet implemented; layer 3 (continuous times) is unvalidated. Layers 1-2 are exact.";
    // negative means "not computed": a relative discrepancy is non-negative, and a
    // NaN sentinel would be undefined behaviour under -ffinite-math-only
    return notComputed;
}


namespace {

    // full round-trip precision; std::to_string would silently truncate branch
    // lengths to six decimals and make the recorded tree not quite the one used
    std::string exact(double x) {

        std::ostringstream os;
        os << std::setprecision(17) << x;
        return os.str();
    }

    const char* kindName(MappedEvent::Kind k) {

        switch (k)
            {
            case MappedEvent::Substitution: return "substitution";
            case MappedEvent::Insertion:    return "insertion";
            case MappedEvent::Deletion:     return "deletion";
            }
        return "unknown";
    }
}

void TKF91StochasticMapping::buildNewick(int i, std::string& s) const {

    if (poIsLeaf[i])
        {
        s += poNode[i]->getName();
        }
    else
        {
        s += "(";
        buildNewick(poLeft[i], s);
        s += ",";
        buildNewick(poRight[i], s);
        s += ")";
        }
    if (poParent[i] >= 0)
        s += ":" + exact(poBranchLength[i]);
}

std::string TKF91StochasticMapping::newickString(void) const {

    std::string s;
    if (numNodes == 0)
        return s;
    buildNewick((int)numNodes - 1, s);       // post-order puts the root last
    s += ";";
    return s;
}

nlohmann::json TKF91StochasticMapping::historyToJson(const MappedHistory& hist, int cognateIndex, long generation) const {

    using nlohmann::json;
    json j;

    j["format"]     = "tkf91-stochastic-map";
    j["version"]    = 2;                            // 2: records carry substitution events
    j["cognate"]    = cognateIndex;
    j["generation"] = generation;

    {
        json t;
        t["numNodes"] = (uint64_t)numNodes;
        t["numTaxa"]  = (uint64_t)numTaxa;

        json parent = json::array(), left = json::array(), right = json::array();
        json isLeaf = json::array(), brlen = json::array();
        json leafRow = json::array(), names = json::array();
        for (size_t i=0; i<numNodes; i++)
            {
            parent.push_back(poParent[i]);
            left.push_back(poLeft[i]);
            right.push_back(poRight[i]);
            isLeaf.push_back(poIsLeaf[i] ? 1 : 0);
            brlen.push_back(poBranchLength[i]);
            leafRow.push_back(poSeqRow[i]);                 // alignment row, -1 for interior
            names.push_back(poIsLeaf[i] ? std::string(poNode[i]->getName()) : std::string());
            }
        t["parent"]       = parent;
        t["left"]         = left;
        t["right"]        = right;
        t["isLeaf"]       = isLeaf;
        t["branchLength"] = brlen;
        t["leafRow"]      = leafRow;
        t["taxonName"]    = names;
        t["newick"]       = newickString();                 // convenience for other tools
        j["tree"] = t;
    }

    // ---- the parameter values the draw was conditioned on ----
    {
        json p;
        p["insertionRate"] = insertionRate;
        p["deletionRate"]  = deletionRate;
        p["equilibriumFrequencies"] = equilibriumFrequencies;
        if (rateMatrix.empty() == false)
            {
            json Q = json::array();
            for (size_t i=0; i<numStates; i++)
                {
                json row = json::array();
                for (size_t k=0; k<numStates; k++)
                    row.push_back(rateMatrix[i * numStates + k]);
                Q.push_back(row);
                }
            p["rateMatrix"] = Q;                            // as used: average rate 1
            }
        j["parameters"] = p;
    }

    j["alignment"] = { {"numColumns", (uint64_t)numSegments}, {"numStates", (uint64_t)numStates} };

    {
        json lins = json::array();
        for (size_t k=0; k<hist.lineages.size(); k++)
            {
            const MappedLineage& L = hist.lineages[k];
            json e;
            e["index"]  = (uint64_t)k;
            e["column"] = L.column;
            e["silent"] = (L.column < 0);

            json pres = json::array();
            for (size_t i=0; i<L.presence.size(); i++)
                pres.push_back(L.presence[i] ? 1 : 0);      // vector<char> -> 0/1, not a string
            e["presence"]  = pres;
            e["nodeState"] = L.nodeState;                   // -1 where the lineage is absent
            lins.push_back(e);
            }
        j["lineages"] = lins;
        j["numSilentLineages"] = (uint64_t)hist.numSilentLineages;
    }

    {
        json evs = json::array();
        for (size_t b=0; b<hist.events.size(); b++)
            {
            for (size_t m=0; m<hist.events[b].size(); m++)
                {
                const MappedEvent& e = hist.events[b][m];
                json o;
                o["branch"]  = (uint64_t)b;
                o["kind"]    = kindName(e.kind);
                o["time"]    = e.time;
                o["lineage"] = e.lineage;
                if (e.fromState >= 0)
                    o["from"] = e.fromState;
                if (e.toState >= 0)
                    o["to"] = e.toState;
                evs.push_back(o);
                }
            }
        j["events"] = evs;
    }

    return j;
}

void TKF91StochasticMapping::appendHistoryToJsonFile(const MappedHistory& hist, const std::string& path, int cognateIndex, long generation) {

    const std::string text = historyToJson(hist, cognateIndex, generation).dump();

    std::error_code ec;
    const bool exists = std::filesystem::exists(path, ec) && std::filesystem::file_size(path, ec) > 0;

    if (exists == false)
        {
        std::ofstream out(path, std::ios::binary | std::ios::trunc);
        if (out.is_open() == false)
            Msg::error("Could not create stochastic-mapping file " + path);
        out << "[\n" << text << "\n]\n";
        if (out.good() == false)
            Msg::error("Failed writing stochastic-mapping file " + path);
        return;
        }

    // locate the closing bracket, which our own writes put within the last few bytes
    std::ifstream in(path, std::ios::binary);
    if (in.is_open() == false)
        Msg::error("Could not open stochastic-mapping file " + path);
    in.seekg(0, std::ios::end);
    const std::streamoff size = in.tellg();
    const std::streamoff window = (size < 64) ? size : 64;

    std::string tail((size_t)window, '\0');
    in.seekg(size - window);
    in.read(&tail[0], window);
    in.close();

    const size_t rel = tail.find_last_of(']');
    if (rel == std::string::npos)
        Msg::error("File " + path + " does not end in a JSON array; refusing to append");
    const std::streamoff cut = (size - window) + (std::streamoff)rel;

    std::filesystem::resize_file(path, (std::uintmax_t)cut, ec);
    if (ec)
        Msg::error("Could not truncate stochastic-mapping file " + path + ": " + ec.message());

    std::ofstream out(path, std::ios::binary | std::ios::app);
    if (out.is_open() == false)
        Msg::error("Could not reopen stochastic-mapping file " + path);
    out << ",\n" << text << "\n]\n";
    if (out.good() == false)
        Msg::error("Failed appending to stochastic-mapping file " + path);
}
