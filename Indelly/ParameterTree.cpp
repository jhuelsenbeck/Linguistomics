#include <cmath>
#include <iomanip>
#include <iostream>
#include "Msg.hpp"
#include "Node.hpp"
#include "NodeSet.hpp"
#include "ParameterTree.hpp"
#include "RandomVariable.hpp"
#include "TransitionProbabilities.hpp"
#include "Tree.hpp"



ParameterTree::ParameterTree(RandomVariable* r, Model* m, std::string treeStr, std::vector<std::string> tNames, double itl) : Parameter(r, m, "tree") {

    updateChangesEigens = false;
    
    betaT = itl;
    
    if (treeStr != "")
        std::cout << "   * Setting up the tree parameter from information in json file " << std::endl;
    else
        std::cout << "   * Setting up the tree parameter with randomly chosen tree" << std::endl;
        
    if (treeStr != "")
        trees[0] = new Tree(treeStr, tNames, betaT, rv);
    else
        trees[0] = new Tree(tNames, betaT, rv);
    trees[1] = new Tree(*trees[0]);
    
    //trees[0]->print("trees[0]");
    //trees[1]->print("trees[1]");
    
    /*std::vector<bool> tMask(tNames.size());
    for (int i=0; i<tMask.size(); i++)
        tMask[i] = false;
    tMask[2] = true;
    tMask[6] = true;
    tMask[7] = true;
    Tree t = Tree(*trees[0], tMask);*/
    
}

ParameterTree::~ParameterTree(void) {

    delete trees[0];
    delete trees[1];
}

void ParameterTree::accept(void) {

    *trees[1] = *trees[0];
}

void ParameterTree::nniArea(std::vector<Node*>& backbone, Node*& incidentNode) {

    Tree* t = trees[0];
    Node* p = NULL;
    do
        {
        p = t->randomNode(rv);
        } while(p->getIsLeaf() == true || p == t->getRoot());
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

std::string ParameterTree::getString(void) {

    std::string str = trees[0]->getNewick();
    return str;
}

double ParameterTree::lnPriorProbability(void) {

    // Prior from:
    // Rannala, B., T. Zhu, and Z. Yang. 2012. Tail paradox, partial identifiability, and
    //    influential priors in Bayesian branch length inference. Molecular Biology and
    //    Evolution 29(1):325-335.
    // with C = 1.0 and alpha = 1.0.

    // joint prior on tree length and branch lengths from
    Tree* t = trees[0];
    double alphaT = 1.0;
    double s = (double)t->getNumTaxa();

    // get a vector with the branch lengths
    double treeLength = t->getTreeLength();
    std::vector<Node*>& dpSeq = trees[0]->getDownPassSequence();
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
    lnP += alphaT * log(betaT) - rv->lnGamma(alphaT) - betaT * treeLength;
    lnP += (alphaT - 1.0) * log(treeLength);
    lnP += rv->lnGamma(2 * s - 4.0 + 1.0);
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

    trees[0]->print();
}

void ParameterTree::reject(void) {

    *trees[0] = *trees[1];
}

double ParameterTree::update(void) {

//    updateNni();
    
    // pick a tree parameter to update
    double u = rv->uniformRv();
    
    if (u < 0.75)
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
    rv->dirichletRv(alphaForward, newProportions);
    
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
    double lnH  = rv->lnDirichletPdf(alphaReverse, oldProportions) - rv->lnDirichletPdf(alphaForward, newProportions); // Hastings Ratio
    lnH += (n - 2) * log(newProportions[1] / oldProportions[1]); // Jacobian (check this)

    // update the transition probabilities
    updateChangesTransitionProbabilities = true;
    TransitionProbabilities& tip = TransitionProbabilities::transitionProbabilties();
    tip.flipActive();
    tip.setNeedsUpdate(true);
    tip.setTransitionProbabilities();

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

    double tuning = log(4.0);
    Tree* t = trees[0];
    
    // temporarily, make all of the branch proportion instance variables
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
    
    // calculate the length of the backbone
    double m = 0.0;
    for (int i=0; i<3; i++)
        m += backbone[i]->getBranchProportion();
        
    // choose a new length for the backbone
    double mPrime = m * exp(tuning*(rv->uniformRv()-0.5));
    double factor = mPrime / m;
    backbone[0]->setBranchProportion( backbone[0]->getBranchProportion()*factor );
    backbone[1]->setBranchProportion( backbone[1]->getBranchProportion()*factor );
    backbone[2]->setBranchProportion( backbone[2]->getBranchProportion()*factor );
    
    // slide the incident node a random amount along the backbone
    if (backbone[1]->isDescendant(incidentNode) == true)
        {
        // the incident node is above the backbone
        }
    else
        {
        // the incident node is below the backbone
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
    
    
    t->print();
    std::cout << "backbone[0]  = " << backbone[0]->getIndex() << std::endl;
    std::cout << "backbone[1]  = " << backbone[1]->getIndex() << std::endl;
    std::cout << "backbone[2]  = " << backbone[2]->getIndex() << std::endl;
    std::cout << "incidentNode = " << incidentNode->getIndex() << std::endl;

    exit(1);
    return log(mPrime) - log(m);
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

    return log(newL) - log(oldL);
}
