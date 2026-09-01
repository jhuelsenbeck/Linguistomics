#include <cmath>
#include <cstdint>
#include <cstring>
#include <limits>
#include <iomanip>
#include <iostream>
#include "Alignment.hpp"
#include "LikelihoodCalculator.hpp"
#include "Msg.hpp"
#include "Node.hpp"
#include "ParameterAlignment.hpp"
#include "ParameterFrequencies.hpp"
#include "ParameterIndelRates.hpp"
#include "ParameterTree.hpp"
#include "TransitionProbabilitiesGpu.hpp"
#include "Tree.hpp"


static inline bool isNotFinite(double x) {

    uint64_t bits;
    std::memcpy(&bits, &x, sizeof(bits));
    return ((bits >> 52) & 0x7FFULL) == 0x7FFULL;      // all-ones exponent: Inf or NaN
}

/* The value must be a real, strictly positive number. Replaces !(x > 0.0),
   which silently stops catching NaN under -ffinite-math-only. */
static inline bool isPositiveFinite(double x) {

    return (isNotFinite(x) == false) && (x > 0.0);
}

LikelihoodCalculator::LikelihoodCalculator(TransitionProbabilities* tpc, ParameterAlignment* a, ParameterTree* t, ParameterIndelRates* r, ParameterFrequencies* f) :
    tiProbs(tpc), myAlignment(a), myTree(t), myIndelRates(r), myFrequencies(f) {

    alignment = myAlignment->getAlignment();
    if (alignment == nullptr)
        Msg::error("Null alignment when initializing the TKF91 calculator");

    taxonMask = myAlignment->getTaxonMask();

    numStates = myAlignment->getNumStates();
    if (numStates == 0)
        Msg::error("Expecting at least one character state in the TKF91 calculator");

    equilibriumFrequencies = new double[numStates];

    tree = myTree->getTree(taxonMask);
    if (tree == nullptr)
        Msg::error("Could not find tree when initializing the TKF91 calculator");

    numSegments = alignment->getNumSegments();

    // initializeTopology establishes numNodes, numTaxa and every derived size
    initializeTopology();
}

LikelihoodCalculator::~LikelihoodCalculator(void) {

    delete [] equilibriumFrequencies;
}

void LikelihoodCalculator::initializeTopology(void) {

    const std::vector<Node*>& postOrder = tree->getPostOrder();

    numNodes = postOrder.size();
    if (numNodes == 0)
        Msg::error("Empty post-order traversal in the TKF91 calculator");

    // a state is a bitmask over post-order positions, so the tree must fit a uint64_t
    if (numNodes > maxNodes)
        Msg::error("The exact TKF91 calculator handles at most " + std::to_string(maxNodes) +
                   " nodes (" + std::to_string(maxNodes / 2) + " taxa); this tree has " +
                   std::to_string(numNodes));

    poNode.assign(numNodes, nullptr);
    poLeft.assign(numNodes, -1);
    poRight.assign(numNodes, -1);
    poSeqRow.assign(numNodes, -1);
    poIsLeaf.assign(numNodes, 0);
    poSubtreeSize.assign(numNodes, 0);
    poSubtree.assign(numNodes, 0);
    poDisjointLower.assign(numNodes, 0);
    poTiMatrix.assign(numNodes, nullptr);

    // map each Node to its post-order position. Held as a member so the buckets
    // survive between evaluations; this runs on every likelihood call.
    std::unordered_map<const Node*,int>& position = nodePosition;
    position.clear();
    position.reserve(numNodes * 2);
    for (size_t i=0; i<numNodes; i++)
        {
        if (postOrder[i] == nullptr)
            Msg::error("Null node at post-order position " + std::to_string(i) +
                       " in the TKF91 calculator");
        poNode[i] = postOrder[i];
        position[postOrder[i]] = (int)i;
        }

    size_t numRoots = 0;
    size_t numLeaves = 0;

    for (size_t i=0; i<numNodes; i++)
        {
        Node* p = poNode[i];

        if (p->getAncestor() == nullptr)
            numRoots++;
        else if (position.find(p->getAncestor()) == position.end())
            Msg::error("Node at post-order position " + std::to_string(i) +
                       " has an ancestor that is not in the traversal; the tree passed to the "
                       "TKF91 calculator is not a single connected tree");

        poIsLeaf[i] = (p->getIsLeaf() == true) ? 1 : 0;

        if (poIsLeaf[i])
            {
            numLeaves++;
            poSeqRow[i] = p->getIndex();
            }
        else
            {
            // children precede their parent in post-order, so they are already placed
            Node* lft = p->getDescendant(0);
            Node* rht = p->getDescendant(1);
            if (lft == nullptr || rht == nullptr)
                Msg::error("Interior node at post-order position " + std::to_string(i) +
                           " has " + std::to_string(p->numDescendants()) +
                           " descendants; the TKF91 calculator expects a rooted binary tree");
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

    if (numRoots != 1)
        Msg::error("Found " + std::to_string(numRoots) + " nodes without an ancestor in the tree "
                   "passed to the TKF91 calculator; expecting exactly one root");

    // The alignment supplies one row per leaf of THIS tree. If a subtree is being
    // scored, its alignment must have been reduced to match it.
    numTaxa = numLeaves;
    if (alignment->getNumTaxa() != numTaxa)
        Msg::error("The tree passed to the TKF91 calculator has " + std::to_string(numTaxa) +
                   " leaves but the alignment has " + std::to_string(alignment->getNumTaxa()) +
                   " rows (post-order length " + std::to_string(numNodes) +
                   ", Tree::getNumNodes() " + std::to_string(tree->getNumNodes()) + ")");

    std::vector<int> seen(numTaxa, 0);
    for (size_t i=0; i<numNodes; i++)
        {
        if (poIsLeaf[i] == 0)
            continue;
        const int row = poSeqRow[i];
        if (row < 0 || (size_t)row >= numTaxa)
            Msg::error("Leaf at post-order position " + std::to_string(i) + " has Node::index " +
                       std::to_string(row) + ", which is not a valid alignment row for a tree with " +
                       std::to_string(numTaxa) + " leaves. If this is a subtree, its leaves are "
                       "probably still carrying indices from the full taxon list.");
        if (seen[row] != 0)
            Msg::error("Alignment row " + std::to_string(row) +
                       " is claimed by more than one leaf in the TKF91 calculator");
        seen[row] = 1;
        }

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

    tkf91Probs.beta.resize(numNodes);
    tkf91Probs.birthProbability.resize(numNodes);
    tkf91Probs.extinctionProbability.resize(numNodes);
    tkf91Probs.homologousProbability.resize(numNodes);
    tkf91Probs.nonHomologousProbability.resize(numNodes);
    fWork.resize(numNodes * numStates);
}

void LikelihoodCalculator::cacheTransitionMatrices(void) {

    for (size_t i=0; i<numNodes; i++)
        {
        Node* p = poNode[i];
        if (p->getAncestor() != nullptr)
            {
            DoubleMatrix& tp = tiProbs->getTransitionProbability(p->getBranchLength());
            poTiMatrix[i] = tp.begin();
            }
        else
            {
            poTiMatrix[i] = nullptr;
            }
        }
}

void LikelihoodCalculator::setBirthDeathProbabilities(void) {

    tkf91Probs.insertionRate = myIndelRates->getInsertionRate();
    tkf91Probs.deletionRate = myIndelRates->getDeletionRate();

    const double lambda = tkf91Probs.insertionRate;
    const double mu     = tkf91Probs.deletionRate;

    // TKF91 requires lambda < mu for a finite equilibrium sequence length. Without
    // it beta and the root factors are meaningless rather than merely extreme.
    if (isPositiveFinite(lambda) == false || isNotFinite(mu) || !(mu > lambda))
        Msg::error("TKF91 requires 0 < lambda < mu");

    tkf91Probs.immortalProbability = 1.0;
    for (size_t i=0; i<numNodes; i++)
        {
        Node* p = poNode[i];
        const double v = p->getBranchLength();
        if (p->getAncestor() == nullptr)
            {
            // Root: tau = infinity by the stationarity assumption, so beta -> 1/mu,
            // giving B = lambda/mu, E = 1 and H = N = 0.
            tkf91Probs.beta[i] = 1.0 / mu;
            tkf91Probs.homologousProbability[i] = 0.0;
            }
        else
            {
            const double expPart = std::exp((lambda - mu) * v);
            tkf91Probs.beta[i] = (1.0 - expPart) / (mu - lambda * expPart);
            tkf91Probs.homologousProbability[i] = std::exp(-mu * v) * (1.0 - lambda * tkf91Probs.beta[i]);
            }
        tkf91Probs.birthProbability[i]      = lambda * tkf91Probs.beta[i];
        tkf91Probs.extinctionProbability[i] = mu * tkf91Probs.beta[i];
        // N = (1 - e^{-mu t} - mu beta)(1 - lambda beta), written here as
        // (1 - mu beta)(1 - B) - H, which is the same expression rearranged.
        tkf91Probs.nonHomologousProbability[i] = (1.0 - mu * tkf91Probs.beta[i]) *
                                                 (1.0 - tkf91Probs.birthProbability[i]) -
                                                 tkf91Probs.homologousProbability[i];
        // the immortal link prefactor of Theorem 1 runs over all nodes, root included
        tkf91Probs.immortalProbability *= (1.0 - tkf91Probs.birthProbability[i]);
        }
}

void LikelihoodCalculator::setStationaryFrequencies(void) {

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
            Msg::error("Stationary frequency vector has the wrong dimension in the TKF91 calculator");
        for (size_t i=0; i<numStates; i++)
            equilibriumFrequencies[i] = x[i];
        }
}

/* Turn the fixed alignment into one emission bitmask per column. Bits are
   post-order positions of the emitting leaves. */
void LikelihoodCalculator::buildColumns(void) {

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
                continue;                       // gap
            if (charState > numStates)
                {
                Msg::error("Unrepresentable character state " + std::to_string(charState) +
                           " at alignment row " + std::to_string(poSeqRow[i]) +
                           ", column " + std::to_string(segment) +
                           " (ambiguity codes are not supported by the TKF91 calculator)");
                }
            emit |= ((uint64_t)1 << i);
            }

        if (emit == 0)
            Msg::error("All-gap column " + std::to_string(segment) + " in the TKF91 alignment");

        columnEmit[segment] = emit;
        }
}

void LikelihoodCalculator::enumerateAlivePatterns(int r, uint64_t emitMask, std::vector<uint64_t>& out) {

    out.clear();

    if (poIsLeaf[r])
        {
        // a leaf birth node is alive by definition, so it must be the emitting leaf
        if ((emitMask >> r) & 1)
            out.push_back((uint64_t)1 << r);
        return;
        }

    std::vector<uint64_t> lftOptions;
    std::vector<uint64_t> rhtOptions;

    for (int pass=0; pass<2; pass++)
        {
        const int c = (pass == 0) ? poLeft[r] : poRight[r];
        std::vector<uint64_t>& options = (pass == 0) ? lftOptions : rhtOptions;
        options.clear();

        const uint64_t needed = emitMask & poSubtree[c];
        if (needed == 0)
            options.push_back(0);               // the whole child subtree dies

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
        {
        for (uint64_t h : rhtOptions)
            out.push_back(rBit | l | h);
        }
}

double LikelihoodCalculator::patternWeight(int r, uint64_t aliveMask, size_t column) {

    const size_t n = numStates;
    const size_t lo = (size_t)r + 1 - poSubtreeSize[r];

    for (size_t i=lo; i<=(size_t)r; i++)
        {
        if (((aliveMask >> i) & 1) == 0)
            continue;

        double* const fi = fWork.data() + i * n;

        if (poIsLeaf[i])
            {
            // f(leaf)[alpha] = 1 for the observed character, 0 otherwise
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
                // the child is labelled "-", contributing p(-|X^gamma) = E_c, and
                // everything below it is also "-", contributing p(-|-) = 1
                const double e = tkf91Probs.extinctionProbability[c];
                for (size_t a=0; a<n; a++)
                    fi[a] *= e;
                continue;
                }

            const double* const fc = fWork.data() + (size_t)c * n;
            const double* const tp = poTiMatrix[c];
            const double hc = tkf91Probs.homologousProbability[c];
            const double nc = tkf91Probs.nonHomologousProbability[c];

            // sum_gamma ( H_c p_c(alpha->gamma) + N_c pi(gamma) ) f(c)[gamma].
            // The N part does not depend on alpha, so it is hoisted.
            double piDotF = 0.0;
            for (size_t g=0; g<n; g++)
                piDotF += equilibriumFrequencies[g] * fc[g];
            const double newTerm = nc * piDotF;

            for (size_t a=0; a<n; a++)
                {
                const double* const row = tp + a * n;   // row = from-state, column = to-state
                double dot = 0.0;
                for (size_t g=0; g<n; g++)
                    dot += row[g] * fc[g];
                fi[a] *= (hc * dot + newTerm);
                }
            }
        }

    // the birth itself: p(B^alpha|-) = B_r pi(alpha)
    const double* const fr = fWork.data() + (size_t)r * n;
    double total = 0.0;
    for (size_t a=0; a<n; a++)
        total += equilibriumFrequencies[a] * fr[a];

    return tkf91Probs.birthProbability[r] * total;
}

void LikelihoodCalculator::buildEventTable(uint64_t emitMask, size_t column, std::vector<TKF91Event>& out) {

    out.clear();

    for (size_t r=0; r<numNodes; r++)
        {
        if ((emitMask & ~poSubtree[r]) != 0)
            continue;                           // some emitting leaf lies outside subtree(r)

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
            out.push_back(e);
            }
        }
}

void LikelihoodCalculator::applyEvents(const std::vector<TKF91Event>& events, const StateMap& in, StateMap& out) {

    out.clear();

    for (const auto& entry : in)
        {
        const uint64_t state = entry.first;
        const double a = entry.second;

        for (const TKF91Event& e : events)
            {
            // The birth node is the only node in the event whose liveness has to
            // be tested. This is the legality condition r(e) in S.
            if ((state & e.birthBit) == 0)
                continue;

            const uint64_t next = (state & ~e.clearMask) | e.aliveMask;
            out[next] += a * e.weight;
            }
        }
}

void LikelihoodCalculator::silentClosure(StateMap& a) {

    closureCurrent = a;

    double totalMass = 0.0;
    for (const auto& entry : a)
        totalMass += std::fabs(entry.second);

    for (size_t iter=0; iter<maxClosureIterations; iter++)
        {
        applyEvents(silentEvents, closureCurrent, closureNext);
        if (closureNext.empty())
            return;

        double addedMass = 0.0;
        for (const auto& entry : closureNext)
            {
            a[entry.first] += entry.second;
            addedMass += std::fabs(entry.second);
            }

        if (addedMass <= closureTolerance * (totalMass > 0.0 ? totalMass : 1.0))
            return;
        totalMass += addedMass;

        closureCurrent.swap(closureNext);
        }

    Msg::error("Silent-event series failed to converge in the TKF91 calculator");
}

double LikelihoodCalculator::computeLnLikelihood(void) {

    tree = myTree->getTree(taxonMask);
    if (tree == nullptr)
        Msg::error("Could not find tree when computing the TKF91 likelihood");
    alignment = myAlignment->getAlignment();
    if (alignment == nullptr)
        Msg::error("Could not find alignment when computing the TKF91 likelihood");

    numSegments = alignment->getNumSegments();
    if (numSegments == 0)
        return 0.0;

    initializeTopology();
    setBirthDeathProbabilities();
    setStationaryFrequencies();
    cacheTransitionMatrices();
    buildColumns();

    // silent events do not depend on the column, so their table is built once
    buildEventTable(0, silentColumn, silentEvents);

    // The chain starts with every node live, carrying the immortal link prefactor
    // prod_n (1 - B_n) of Theorem 1.
    const uint64_t fullState = (numNodes == 64) ? ~(uint64_t)0 : (((uint64_t)1 << numNodes) - 1);
    alphaCurrent.clear();
    alphaCurrent[fullState] = 1.0;

    double lnScale = std::log(tkf91Probs.immortalProbability);

    silentClosure(alphaCurrent);

    for (size_t segment=0; segment<numSegments; segment++)
        {
        buildEventTable(columnEmit[segment], segment, columnEvents);

        applyEvents(columnEvents, alphaCurrent, alphaNext);
        if (alphaNext.empty())
            return impossibleLnLikelihood;

        alphaCurrent.swap(alphaNext);
        silentClosure(alphaCurrent);

        double biggest = 0.0;
        for (const auto& entry : alphaCurrent)
            {
            const double m = std::fabs(entry.second);
            if (m > biggest)
                biggest = m;
            }
        if (isPositiveFinite(biggest) == false)
            return impossibleLnLikelihood;

        const double inv = 1.0 / biggest;
        for (auto& entry : alphaCurrent)
            entry.second *= inv;
        lnScale += std::log(biggest);
        }

    // every state is joined to the end state by an edge of probability one
    double total = 0.0;
    for (const auto& entry : alphaCurrent)
        total += entry.second;

    if (isPositiveFinite(total) == false)
        return impossibleLnLikelihood;

    return lnScale + std::log(total);
}
