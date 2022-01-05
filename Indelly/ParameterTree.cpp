#include <cmath>
#include <iomanip>
#include <iostream>
#include <set>
#include "Alignment.hpp"
#include "Model.hpp"
#include "Msg.hpp"
#include "Node.hpp"
#include "NodeSet.hpp"
#include "ParameterTree.hpp"
#include "Probability.hpp"
#include "RandomVariable.hpp"
#include "TransitionProbabilities.hpp"
#include "Tree.hpp"

#undef DEBUG_LOCAL



ParameterTree::ParameterTree(RandomVariable* r, Model* m, std::string treeStr, std::vector<std::string> tNames, std::vector<Alignment*>& wordAlignments, double itl) : Parameter(r, m, "tree") {

    updateChangesRateMatrix = false;
    betaT = itl;
    
    if (treeStr != "")
        std::cout << "   * Setting up the tree parameter from information in json file " << std::endl;
    else
        std::cout << "   * Setting up the tree parameter with randomly chosen tree" << std::endl;
        
    if (treeStr != "")
        fullTree.trees[0] = new Tree(treeStr, tNames, betaT, rv);
    else
        fullTree.trees[0] = new Tree(tNames, betaT, rv);
    fullTree.trees[1] = new Tree(*fullTree.trees[0]);
    
    // initialize the subtrees
    initializeSubtrees(wordAlignments);
    
    if (checkSubtreeCompatibility(0.0001) == false)
        Msg::error("Subtrees are not compatible with the canonical tree");
}

ParameterTree::~ParameterTree(void) {

    delete fullTree.trees[0];
    delete fullTree.trees[1];
    clearSubtrees();
}

void ParameterTree::accept(void) {

    *(fullTree.trees[1]) = *(fullTree.trees[0]);
    for (std::map<RbBitSet,TreePair>::iterator it = subTrees.begin(); it != subTrees.end(); it++)
        *(it->second.trees[1]) = *(it->second.trees[0]);
}

bool ParameterTree::checkSubtreeCompatibility(double tolerance) {

    std::map<TaxonPair,double,CompTaxonPair> fullTreeDistances = fullTree.trees[0]->pairwiseDistances();
    int n = (int)fullTree.trees[0]->getTaxonNames().size();
    if (fullTreeDistances.size() != n * (n-1) / 2)
        Msg::error("Incomplete pairwise distance map for full tree");
    
    for (std::map<RbBitSet,TreePair>::iterator i = subTrees.begin(); i != subTrees.end(); i++)
        {
        std::map<TaxonPair,double,CompTaxonPair> subTreeDistances = i->second.trees[0]->pairwiseDistances();
        n = (int)i->second.trees[0]->getTaxonNames().size();
        if (subTreeDistances.size() != n * (n-1) / 2)
            Msg::error("Incomplete pairwise distance map for sub tree");
        
        for (std::map<TaxonPair,double,CompTaxonPair>::iterator j = subTreeDistances.begin(); j != subTreeDistances.end(); j++)
            {
            double stD = j->second;
            
            std::map<TaxonPair,double,CompTaxonPair>::iterator it = fullTreeDistances.find(j->first);
            if (it != fullTreeDistances.end())
                {
                double ftD = it->second;
                double diff = fabs(stD - ftD);
                if (diff > tolerance)
                    return false;
                }
            else
                {
                Msg::error("Could not find taxon pair");
                }
            }
        }
    return true;
}

void ParameterTree::clearSubtrees(void) {

    for (std::map<RbBitSet,TreePair>::iterator it = subTrees.begin(); it != subTrees.end(); it++)
        {
        delete it->second.trees[0];
        delete it->second.trees[1];
        }
    subTrees.clear();
}

int ParameterTree::countMaskBits(std::vector<bool>& m) {

    int n = 0;
    for (size_t i=0; i<m.size(); i++)
        {
        if (m[i] == true)
            n++;
        }
    return n;
}

void ParameterTree::nniArea(std::vector<Node*>& backbone, Node*& incidentNode) {

    Tree* t = fullTree.trees[0];
    Node* r = t->getRoot();
    Node* p = NULL;
    bool goodNode = true;
    do
        {
        p = t->randomNode(rv);
        goodNode = true;
        if (p->getIsLeaf() == true || p == r)
            goodNode = false;
        else if (p->getAncestor() == r)
            goodNode = false;
        } while(goodNode == false);
    Node* pAnc = p->getAncestor();
    
    if (pAnc == NULL)
        Msg::error("Problem choosing random backbone");
    std::vector<Node*> pDes = p->getDescendantsVector();
    if (pDes.size() != 2)
        Msg::error("Problem choosing random backbone");
    std::vector<Node*> pAncDes = pAnc->getDescendantsVector();
    if (pAncDes.size() != 2)
        Msg::error("Problem choosing random backbone");

    std::vector<Node*> possibleIncidentNodes(2);
    Node* a = pDes[(int)(rv->uniformRv()*2)];
    if (a == pDes[0])
        possibleIncidentNodes[0] = pDes[1];
    else
        possibleIncidentNodes[0] = pDes[0];
    Node* b = NULL;
    if (rv->uniformRv() < 0.5)
        {
        b = pAnc;
        if (pAncDes[0] == p)
            possibleIncidentNodes[1] = pAncDes[1];
        else
            possibleIncidentNodes[1] = pAncDes[0];
        }
    else
        {
        if (pAncDes[0] == p)
            b = pAncDes[1];
        else
            b = pAncDes[0];
        possibleIncidentNodes[1] = pAnc;
        }
        
    // we always assume this order, from tip to root
    backbone.clear();
    backbone.push_back(a);
    backbone.push_back(p);
    backbone.push_back(b);
    incidentNode = possibleIncidentNodes[(int)(rv->uniformRv()*2)];
}

Tree* ParameterTree::getActiveTree(RbBitSet& mask) {

    if (mask.getNumberSetBits() == mask.size())
        return fullTree.trees[0];
        
    std::map<RbBitSet,TreePair>::iterator it = subTrees.find(mask);
    if (it != subTrees.end())
        return it->second.trees[0];
    return NULL;
}

Tree* ParameterTree::getActiveTree(const RbBitSet& mask) {

    if (mask.getNumberSetBits() == mask.size())
        return fullTree.trees[0];
        
    std::map<RbBitSet,TreePair>::iterator it = subTrees.find(mask);
    if (it != subTrees.end())
        return it->second.trees[0];
    return NULL;
}

std::string ParameterTree::getJsonString(void) {

    std::string str = "\"Tree\": \"" + fullTree.trees[0]->getNewick(20) + "\"";
    return str;
}

std::string ParameterTree::getString(void) {

    std::string str = fullTree.trees[0]->getNewick(5);
    return str;
}

void ParameterTree::initializeSubtrees(std::vector<Alignment*>& alns) {
        
    // remove all of the current subtrees
    clearSubtrees();
    
    // loop over all of the alignments, finding those that require subtrees because they only have a subset of the
    // full canonical list of taxa
    for (int i=0; i<alns.size(); i++)
        {
        // get the taxon mask and only proceed if the taxon list is incomplete
        std::vector<bool> m = alns[i]->getTaxonMask();
        if (countMaskBits(m) == m.size())
            continue;
        RbBitSet mask(m);

        std::map<RbBitSet,TreePair>::iterator it = subTrees.find(mask);
        if (it == subTrees.end())
            {
            // build the subtree
            Tree* t0 = new Tree(*fullTree.trees[0], mask);
            Tree* t1 = new Tree(*t0);
            TreePair pair(t0, t1);
            subTrees.insert( std::make_pair(mask,pair) );
            }
        }
        
#   if 0
    int i = 0;
    for (std::map<RbBitSet,TreePair>::iterator it = subTrees.begin(); it != subTrees.end(); it++)
        {
        std::cout << "Subtree: " << i << " (" << it->first << ")" << std::endl;
        it->second.trees[0]->print();
        i++;
        }
#   endif
}

double ParameterTree::lnPriorProbability(void) {

    // Prior from:
    // Rannala, B., T. Zhu, and Z. Yang. 2012. Tail paradox, partial identifiability, and
    //    influential priors in Bayesian branch length inference. Molecular Biology and
    //    Evolution 29(1):325-335.
    // with C = 1.0 and alpha = 1.0.

    // joint prior on tree length and branch lengths from
    Tree* t = fullTree.trees[0];
    double alphaT = 1.0;
    double s = (double)t->getNumTaxa();

    // get a vector with the branch lengths
    double treeLength = t->getTreeLength();
    std::vector<Node*>& dpSeq = fullTree.trees[0]->getDownPassSequence();
    std::vector<double> branchLengths;
    for (int i=0; i<dpSeq.size(); i++)
        {
        Node* p = dpSeq[i];
        if (p->getAncestor() != t->getRoot())
            {
            if (p == t->getRoot())
                {
                std::set<Node*,CompNode>& rootDes = p->getDescendants()->getNodes();
                double len = 0.0;
                for (Node* n : rootDes)
                    len += n->getBranchProportion() * treeLength;
                branchLengths.push_back(len);
                }
            else
                {
                branchLengths.push_back( treeLength * p->getBranchProportion() );
                }
            }
        }
        
    double lnP = 0.0;
    lnP += alphaT * log(betaT) - Probability::Helper::lnGamma(alphaT) - betaT * treeLength;
    lnP += (alphaT - 1.0) * log(treeLength);
    lnP += Probability::Helper::lnGamma(2 * s - 4.0 + 1.0);
    lnP += (-2.0 * s + 4.0) * log(treeLength);
            
    return lnP;
}

void ParameterTree::normalize(std::vector<double>& vec, double minVal) {
    
    // find entries with values that are too small
    int numTooSmall = 0;
    double sum = 0.0;
    for (int i=0; i<vec.size(); i++)
        {
        if (vec[i] < minVal)
            {
            numTooSmall++;
            }
        else
            sum += vec[i];
        }
        
    double factor = (1.0 - numTooSmall * minVal) / sum;
    for (int i=0; i<vec.size(); i++)
        {
        if (vec[i] < minVal)
            vec[i] = minVal;
        else
            vec[i] *= factor;
        }
        
    
#   if 0
    sum = 0.0;
    for (int i=0; i<vec.size(); i++)
        sum += vec[i];
    if ( fabs(1.0 - sum) > 0.000001)
        std::cout << "Problem normalizing vector " << std::fixed << std::setprecision(20) << sum << std::endl;
#   endif
}

void ParameterTree::print(void) {

    fullTree.trees[0]->print();
}

void ParameterTree::printNewick(void) {

    std::cout << fullTree.trees[0]->getNewick(6) << ";" << std::endl;
    for (std::map<RbBitSet,TreePair>::iterator it = subTrees.begin(); it != subTrees.end(); it++)
        std::cout << it->second.trees[0]->getNewick(6) << ";" << std::endl;
}

void ParameterTree::reject(void) {

    *(fullTree.trees[0]) = *(fullTree.trees[1]);
    for (std::map<RbBitSet,TreePair>::iterator it = subTrees.begin(); it != subTrees.end(); it++)
        *(it->second.trees[0]) = *(it->second.trees[1]);
    modelPtr->flipActiveLikelihood();
}

double ParameterTree::update(void) {
    
    // pick a tree parameter to update
    double u = rv->uniformRv();
    if (u <= 0.50)
        return updateNni();
    else if (u > 0.50 && u <= 0.75)
        return updateBrlenProportions();
    else
        return updateTreeLength();
    
    return 0.0;
}

double ParameterTree::updateBrlenProportions(void) {

    lastUpdateType = "branch proportions";
    
    // tuning parameter for move
    double alpha0 = 100.0;
    
    // pick a branch at random (note that the two branches
    // from the root are treated as a single branch)
    Tree* t = getActiveTree();
    std::vector<Node*>& nodes = t->getDownPassSequence();
    Node* nde = NULL;
    do
        {
        nde = nodes[(int)(rv->uniformRv()*nodes.size())];
        if (nde == t->getRoot())
            nde = NULL;
        } while(nde == NULL);
        
    // determine if the branch is one of the two branches
    // incident to the root
    bool isRootBranch = false;
    if (nde->getAncestor() == t->getRoot())
        isRootBranch = true;

    // get original branch proportions
    std::vector<double> oldProportions(2);
    if (isRootBranch == false)
        {
        oldProportions[0] = nde->getBranchProportion();
        }
    else
        {
        std::set<Node*,CompNode>& des = t->getRoot()->getDescendants()->getNodes();
        oldProportions[0] = 0.0;
        for (Node* n : des)
            oldProportions[0] += n->getBranchProportion();
        }
    oldProportions[1] = 1.0 - oldProportions[0];

    // propose new branch proportions
    std::vector<double> alphaForward(2);
    for (int i=0; i<2; i++)
        alphaForward[i] = oldProportions[i] * alpha0;
    std::vector<double> newProportions(2);
    Probability::Dirichlet::rv(rv, alphaForward, newProportions);
    
    // check the proposed values
    normalize(newProportions, 0.000001);
    
    // get reverse move information
    std::vector<double> alphaReverse(2);
    for (int i=0; i<2; i++)
        alphaReverse[i] = newProportions[i] * alpha0;
    
    // update branch proportions on the tree
    if (isRootBranch == false)
        {
        nde->setBranchProportion(newProportions[0]);
        double factor = newProportions[1] / oldProportions[1];
        for (int i=0; i<nodes.size(); i++)
            {
            Node* p = nodes[i];
            if (p != nde && p != t->getRoot())
                {
                double x = p->getBranchProportion();
                p->setBranchProportion( x * factor );
                }
            }
        }
    else
        {
        // update root branches
        std::vector<Node*> rootDes = t->getRoot()->getDescendantsVector();
        if (rootDes.size() != 2)
            Msg::error("Expecting two descendants of the root node");
        double factor0 = newProportions[0] / oldProportions[0];
        double factor1 = newProportions[1] / oldProportions[1];
        for (int i=0; i<nodes.size(); i++)
            {
            Node* p = nodes[i];
            if (p != t->getRoot())
                {
                double x = p->getBranchProportion();
                if (p == rootDes[0] || p == rootDes[1])
                    p->setBranchProportion( x * factor0 );
                else
                    p->setBranchProportion( x * factor1 );
                }
            }
        }
    
    // calculate Hastings ratio, including the Jacobian
    int n = 2 * t->getNumTaxa() - 3;
    double lnH  = Probability::Dirichlet::lnPdf(alphaReverse, oldProportions) - Probability::Dirichlet::lnPdf(alphaForward, newProportions); // Hastings Ratio
    lnH += (n - 2) * log(newProportions[1] / oldProportions[1]); // Jacobian (check this)

    // update the subtrees (before updating transition probabilities)
    updateSubtrees();
    
    // update the transition probabilities
    updateChangesTransitionProbabilities = true;
    TransitionProbabilities& tip = TransitionProbabilities::transitionProbabilties();
    tip.flipActive();
    tip.setNeedsUpdate(true);
    tip.setTransitionProbabilities();

    modelPtr->setUpdateLikelihood();
    modelPtr->flipActiveLikelihood();

#   if 0
    std::vector<double> proposedBrlens;
    for (int i=0; i<nodes.size(); i++)
        {
        Node* p = nodes[i];
        proposedBrlens.push_back(p->getBranchProportion());
        }
    std::cout << "(";
    std::cout << std::fixed << std::setprecision(7);
    for (int i=0; i<proposedBrlens.size(); i++)
        std::cout << proposedBrlens[i] << " ";
    std::cout << ")" << std::endl;
#   endif

    return lnH;
}

double ParameterTree::updateNni(void) {

    lastUpdateType = "local";

    double tuning = log(2.0);
    Tree* t = fullTree.trees[0];
#   if defined(DEBUG_LOCAL)
    t->print("BEFORE");
#   endif
    
    // temporarily, make all of the branch proportions
    // for the tree equal to the branch length
    double treeLen = t->getTreeLength();
    std::vector<Node*> dpSeq = t->getDownPassSequence();
    for (int i=0; i<dpSeq.size(); i++)
        {
        double bp = dpSeq[i]->getBranchProportion();
        dpSeq[i]->setBranchProportion(bp*treeLen);
        }
    
    // randomly choose area of rearrangement
    std::vector<Node*> backbone;
    Node* incidentNode = NULL;
    nniArea(backbone, incidentNode);
    Node* p = backbone[1];
    Node* pAnc = p->getAncestor();
    Node* a = backbone[0];
    Node* b = backbone[2];
    if (b == pAnc && b->getAncestor() == t->getRoot())
        {
        std::vector<Node*> rDes = t->getRoot()->getDescendantsVector();
        double rBrlen = rDes[0]->getBranchProportion() + rDes[1]->getBranchProportion();
        rDes[0]->setBranchProportion(0.0);
        rDes[1]->setBranchProportion(0.0);
        b->setBranchProportion(rBrlen);
        }
        
    // get length of backbone and randomly change it
    double mOld = a->getBranchProportion() + p->getBranchProportion() + b->getBranchProportion();
    double mNew = mOld * exp(tuning*(rv->uniformRv()-0.5));
    double factor = mNew / mOld;
    a->setBranchProportion( a->getBranchProportion() * factor );
    p->setBranchProportion( p->getBranchProportion() * factor );
    b->setBranchProportion( b->getBranchProportion() * factor );

    // rearrange the tree
    if (incidentNode == pAnc)
        {
        // the backbone must be rotated around incident
        double newPt = rv->uniformRv() * mNew;
        if (newPt < a->getBranchProportion())
            {
            p->removeDescendant(a);
            p->addDescendant(b);
            pAnc->removeDescendant(b);
            pAnc->addDescendant(a);
            a->setAncestor(pAnc);
            b->setAncestor(p);
            b->setBranchProportion( b->getBranchProportion() + p->getBranchProportion() );
            p->setBranchProportion(newPt);
            a->setBranchProportion( a->getBranchProportion() - newPt );
            }
        else
            {
            // no topology change, but adjust branch lengths for pAnc descendants
            newPt -= a->getBranchProportion();
            double sum = b->getBranchProportion() + p->getBranchProportion();
            b->setBranchProportion(newPt);
            p->setBranchProportion(sum-newPt);
            }
        }
    else
        {
        // detach and reattach incident as it's above backbone
        Node* incidentAnc = incidentNode->getAncestor();
        Node* incidentSis = incidentNode->getSisterNode();
        Node* incidentAncAnc = incidentAnc->getAncestor();
        incidentAncAnc->removeDescendant(incidentAnc);
        incidentAncAnc->addDescendant(incidentSis);
        incidentSis->setAncestor(incidentAncAnc);
        incidentSis->setBranchProportion( incidentSis->getBranchProportion() + incidentAnc->getBranchProportion() );
        incidentAnc->removeDescendant(incidentSis);
        
        double newPt = rv->uniformRv() * mNew;
        if (newPt < a->getBranchProportion())
            {
            // reattach along branch a
            Node* aAnc = a->getAncestor();
            aAnc->removeDescendant(a);
            aAnc->addDescendant(incidentAnc);
            incidentAnc->addDescendant(a);
            incidentAnc->setAncestor(aAnc);
            a->setAncestor(incidentAnc);
            a->setBranchProportion( a->getBranchProportion() - newPt );
            incidentAnc->setBranchProportion( newPt );
            }
        else
            {
            Node* x = NULL;
            if (a->getAncestor() == p)
                x = p;
            else
                {
                if (pAnc == b)
                    x = pAnc;
                else
                    x = b;
                }
            Node* xAnc = x->getAncestor();
            xAnc->removeDescendant(x);
            xAnc->addDescendant(incidentAnc);
            incidentAnc->addDescendant(x);
            incidentAnc->setAncestor(xAnc);
            x->setAncestor(incidentAnc);
            double xPt = newPt - a->getBranchProportion();
            x->setBranchProportion( x->getBranchProportion() - xPt );
            incidentAnc->setBranchProportion( xPt );
            }
        }
    
    // reinitialize the down pass sequence for the tree, in case the topology has changed
    t->initializeDownPassSequence();
    
    // convert back to branch proportions and a tree length
    dpSeq = t->getDownPassSequence();
    treeLen = 0.0;
    for (int i=0; i<dpSeq.size(); i++)
        treeLen +=  dpSeq[i]->getBranchProportion();
    for (int i=0; i<dpSeq.size(); i++)
        dpSeq[i]->setBranchProportion( dpSeq[i]->getBranchProportion()/treeLen );
    t->setTreeLength(treeLen);
    
    // break the branch at the root evenly in two
    std::vector<Node*> rDes = t->getRoot()->getDescendantsVector();
    double rBrlenSum = rDes[0]->getBranchProportion() + rDes[1]->getBranchProportion();
    rDes[0]->setBranchProportion(rBrlenSum*0.5);
    rDes[1]->setBranchProportion(rBrlenSum*0.5);
    
    // update the subtrees (before updating transition probabilities)
    updateSubtrees();

    // update the transition probabilities
    updateChangesTransitionProbabilities = true;
    TransitionProbabilities& tip = TransitionProbabilities::transitionProbabilties();
    tip.flipActive();
    tip.setNeedsUpdate(true);
    tip.setTransitionProbabilities();

    modelPtr->setUpdateLikelihood();
    modelPtr->flipActiveLikelihood();

#   if defined(DEBUG_LOCAL)
    t->print("AFTER");
    std::cout << "backbone[0]  = " << backbone[0]->getIndex() << std::endl;
    std::cout << "backbone[1]  = " << backbone[1]->getIndex() << std::endl;
    std::cout << "backbone[2]  = " << backbone[2]->getIndex() << std::endl;
    std::cout << "incidentNode = " << incidentNode->getIndex() << std::endl;
#   endif

#   if 0
    if (checkSubtreeCompatibility(0.0001) == false)
        Msg::error("Subtrees are not compatible with the canonical tree");
#   endif

    return 3.0 * log(factor);
}

double ParameterTree::updateSpr(void) {

    return 0.0;
}

double ParameterTree::updateTreeLength(void) {

    lastUpdateType = "tree length";

    // update the tree length
    Tree* t = getActiveTree();
    double tuning = log(2.0);
    double oldL = t->getTreeLength();
    double newL = oldL * exp(tuning * (rv->uniformRv()-0.5));
    t->setTreeLength(newL);

    // update the transition probabilities
    updateChangesTransitionProbabilities = true;
    TransitionProbabilities& tip = TransitionProbabilities::transitionProbabilties();
    tip.flipActive();
    tip.setNeedsUpdate(true);
    tip.setTransitionProbabilities();

    modelPtr->setUpdateLikelihood();
    modelPtr->flipActiveLikelihood();

    return log(newL) - log(oldL);
}

double ParameterTree::updateTopologyFromPrior(void) {

    lastUpdateType = "random tree";
    
    Tree* t = getActiveTree();
    std::vector<std::string> tNames;
    std::vector<std::string>& tNamesToCopy = t->getTaxonNames();
    for (int i=0; i<tNamesToCopy.size(); i++)
        tNames.push_back(tNamesToCopy[i]);
        
    // update the topology and calculate proposal probability
    double lnP1 = lnPriorProbability();
    t->buildRandomTree(tNames, betaT, rv);
    double lnP2 = lnPriorProbability();

    // update the transition probabilities
    updateChangesTransitionProbabilities = true;
    TransitionProbabilities& tip = TransitionProbabilities::transitionProbabilties();
    tip.flipActive();
    tip.setNeedsUpdate(true);
    tip.setTransitionProbabilities();
    modelPtr->setUpdateLikelihood();
    modelPtr->flipActiveLikelihood();

    return lnP1 - lnP2;
}

double ParameterTree::updateBranchlengthsFromPrior(void) {

    lastUpdateType = "random branch lengths";

    Tree* t = getActiveTree();
    double lnP1 = lnPriorProbability();
    std::vector<Node*> dpSeq = t->getDownPassSequence();
    double sum = 0.0;
    for (int i=0; i<dpSeq.size(); i++)
        {
        Node* p = dpSeq[i];
        if (p != t->getRoot())
            {
            p->setBranchProportion( Probability::Exponential::rv(rv, 1.0) );
            sum += p->getBranchProportion();
            }
        else
            p->setBranchProportion(0.0);
        }
    for (int i=0; i<dpSeq.size(); i++)
        {
        Node* p = dpSeq[i];
        p->setBranchProportion( p->getBranchProportion() / sum );
        }

    double treeLength = Probability::Gamma::rv(rv, 1.0, betaT);
    t->setTreeLength(treeLength);

    double lnP2 = lnPriorProbability();

    return lnP1 - lnP2;
}

void ParameterTree::updateSubtrees(void) {

    Tree* t = fullTree.trees[0];
    for (std::map<RbBitSet,TreePair>::iterator it = subTrees.begin(); it != subTrees.end(); it++)
        {
        Tree* st = it->second.trees[0];
        RbBitSet bs = RbBitSet(it->first);
        st->makeSubtree(*t, bs);
        }
}
