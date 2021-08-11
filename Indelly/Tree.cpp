#include <iomanip>
#include <iostream>
#include <istream>
#include <sstream>
#include <vector>
#include "Msg.hpp"
#include "Node.hpp"
#include "NodeSet.hpp"
#include "RandomVariable.hpp"
#include "Tree.hpp"



Tree::Tree(Tree& t) {

    clone(t);
}

Tree::Tree(Tree& t, std::vector<bool> taxonMask) {

    // clone the tree
    if (t.getNumTaxa() != taxonMask.size())
        Msg::error("Cannot create pruned tree with ill-formated taxon mask");
    clone(t);
    
    // mark taxa to include
    std::vector<Node*> dpSeq = this->downPassSequence;
    for (int n=0; n<dpSeq.size(); n++)
        {
        Node* p = dpSeq[n];
        if (p->getIsLeaf() == true && taxonMask[p->getIndex()] == true)
            p->setFlag(true);
        else
            p->setFlag(false);
        }
        
    // mark internal nodes to include
    for (int n=0; n<dpSeq.size(); n++)
        {
        Node* p = dpSeq[n];
        if (p->getIsLeaf() == false)
            {
            std::set<Node*,CompNode>& des = p->getDescendants()->getNodes();
            for (Node* nde : des)
                {
                if (nde->getFlag() == true)
                    {
                    p->setFlag(true);
                    break;
                    }
                }
            }
        }
        
    // collect nodes to remove
    std::vector<Node*> nodesDeadToMe;
    for (int n=0; n<dpSeq.size(); n++)
        {
        Node* p = dpSeq[n];
        if (p->getFlag() == false)
            nodesDeadToMe.push_back(p);
        }
        
    // perform reconnections for nodes to be removed
    for (int i=0; i<nodesDeadToMe.size(); i++)
        {
        Node* p = nodesDeadToMe[i];
        Node* pAnc = p->getAncestor();
        if (pAnc != NULL)
            pAnc->removeDescendant(p);
        std::set<Node*,CompNode>& des = p->getDescendants()->getNodes();
        for (Node* nde : des)
            {
            nde->setAncestor(pAnc);
            nde->setBranchProportion( nde->getBranchProportion() + p->getBranchProportion() );
            if (pAnc != NULL)
                pAnc->addDescendant(nde);
            }
        }
        
    // prune nodes from the root
    for (int n=(int)dpSeq.size()-1; n>=0; n--)
        {
        Node* p = dpSeq[n];
        if (p->numDescendants() == 1)
            {
            Node* pAnc = p->getAncestor();
            std::set<Node*,CompNode>& des = p->getDescendants()->getNodes();
            std::set<Node*,CompNode>::iterator it = des.begin();
            if (pAnc == NULL)
                {
                (*it)->setAncestor(NULL);
                (*it)->setBranchProportion(0.0);
                root = (*it);
                }
            else
                {
                pAnc->removeDescendant(p);
                pAnc->addDescendant((*it));
                (*it)->setAncestor(pAnc);
                (*it)->setBranchProportion( (*it)->getBranchProportion() + p->getBranchProportion() );
                }
            p->removeDescendants();
            p->setAncestor(NULL);
            nodesDeadToMe.push_back(p);
            }
        }
    
    initializeDownPassSequence();
        
    t.print("original tree");
    print("pruned tree");
        
}

Tree::Tree(std::vector<std::string> tNames, double betaT, RandomVariable* rv) {

    // read the tree file
    root = NULL;
    numTaxa = 0;
    taxonNames = tNames;

    // start with a simple two-species tree
    Node* n1 = addNode();
    Node* n2 = addNode();
    n1->setName(tNames[0]);
    n2->setName(tNames[1]);
    n1->setIsLeaf(true);
    n2->setIsLeaf(true);
    root = addNode();
    root->addDescendant(n1);
    root->addDescendant(n2);
    n1->setAncestor(root);
    n2->setAncestor(root);
    
    for (int i=2; i<tNames.size(); i++)
        {
        Node* newTip = addNode();
        Node* newInt = addNode();
        newTip->setName(tNames[i]);
        newTip->setIsLeaf(true);
        
        Node* p = nodes[(int)(rv->uniformRv()*nodes.size())];
        Node* pAnc = p->getAncestor();
        
        if (pAnc == NULL)
            {
            p->setAncestor(newInt);
            newTip->setAncestor(newInt);
            newInt->addDescendant(p);
            newInt->addDescendant(newTip);
            newInt->setAncestor(NULL);
            root = newInt;
            }
        else
            {
            p->setAncestor(newInt);
            pAnc->removeDescendant(p);
            pAnc->addDescendant(newInt);
            newInt->addDescendant(p);
            newInt->addDescendant(newTip);
            newInt->setAncestor(pAnc);
            newTip->setAncestor(newInt);
            }
        }

    // initialize down pass sequence
    initializeDownPassSequence();

    // reindex interior nodes
    int intIdx = numTaxa;
    for (int i=0; i<downPassSequence.size(); i++)
        {
        Node* p = downPassSequence[i];
        if (p->getIsLeaf() == false)
            p->setIndex(intIdx++);
        }

    // initialize branch lengths
    // first, initialize all branch lengths, except the root node
    for (int i=0; i<downPassSequence.size(); i++)
        {
        Node* p = downPassSequence[i];
        if (p != root)
            p->setBranchProportion(rv->exponentialRv(1.0));
        else
            p->setBranchProportion(0.0);
        }
    // make certain the two branches incident to the root are the
    // same in length and considered one branch
    std::vector<Node*> rootDes = root->getDescendantsVector();
    if (rootDes.size() != 2)
        Msg::error("Expecting two descendants of the root node");
    double x = rootDes[0]->getBranchProportion() + rootDes[1]->getBranchProportion();
    rootDes[0]->setBranchProportion(x/4.0);
    rootDes[1]->setBranchProportion(x/4.0);

    // rescale so branch proportions sum to 1.0
    double sum = 0.0;
    for (int i=0; i<downPassSequence.size(); i++)
        {
        Node* p = downPassSequence[i];
        sum += p->getBranchProportion();
        }
    for (int i=0; i<downPassSequence.size(); i++)
        {
        Node* p = downPassSequence[i];
        p->setBranchProportion( p->getBranchProportion()/sum );
        }
        
    treeLength = rv->gammaRv(1.0, betaT);
        

}

Tree::Tree(std::string treeStr, std::vector<std::string> tNames, double betaT, RandomVariable* rv) {

    // read the tree file
    Node* p = NULL;
    root = NULL;
    numTaxa = 0;
    taxonNames = tNames;

    std::vector<std::string> tokens = tokenizeTreeString(treeStr);
    bool readingBrlen = false;
    for (int i=0; i<tokens.size(); i++)
        {
        std::string token = tokens[i];
        if (token == "(")
            {
            if (p == NULL)
                {
                Node* newNode = addNode();
                p = newNode;
                root = newNode;
                }
            else
                {
                Node* newNode = addNode();
                p->addDescendant(newNode);
                newNode->setAncestor(p);
                p = newNode;
                }
            readingBrlen = false;
            }
        else if (token == ")" || token == ",")
            {
            if (p->getAncestor() != NULL)
                p = p->getAncestor();
            readingBrlen = false;
            }
        else if (token == ";")
            {
            readingBrlen = false;
            }
        else if (token == ":")
            {
            readingBrlen = true;
            }
        else
            {
            if (readingBrlen == false)
                {
                int n = 0;
                for (std::string s : tNames)
                    {
                    if (s == token)
                        break;
                    n++;
                    }
                
                Node* newNode = addNode();
                p->addDescendant(newNode);
                newNode->setAncestor(p);
                newNode->setName(token);
                newNode->setIsLeaf(true);
                newNode->setIndex(n);
                p = newNode;
                numTaxa++;
                }
            else
                {
                double brlen = std::stod(token);
                p->setBranchProportion(brlen);
                }
            readingBrlen = false;
            }
        }

    // initialize down pass sequence
    initializeDownPassSequence();

    // reindex interior nodes
    int intIdx = numTaxa;
    for (int i=0; i<downPassSequence.size(); i++)
        {
        p = downPassSequence[i];
        if (p->getIsLeaf() == false)
            p->setIndex(intIdx++);
        }

    // check for branch lengths, randomly initializing from prior
    // if they are anot all initialized from the Newick string
    bool branchLengthsPresent = true;
    for (int i=0; i<downPassSequence.size(); i++)
        {
        p = downPassSequence[i];
        if (p != root && p->getBranchProportion() < 0.0001)
            branchLengthsPresent = false;
        }
    if (branchLengthsPresent == false)
        {
        // initialize branch lengths from uniform Dirichlet
        // first, initialize all branch lengths, except the root node
        for (int i=0; i<downPassSequence.size(); i++)
            {
            p = downPassSequence[i];
            if (p != root)
                p->setBranchProportion(rv->exponentialRv(1.0));
            else
                p->setBranchProportion(0.0);
            }
        // make certain the two branches incident to the root are the
        // same in length and considered one branch
        std::vector<Node*> rootDes = root->getDescendantsVector();
        if (rootDes.size() != 2)
            Msg::error("Expecting two descendants of the root node");
        double x = rootDes[0]->getBranchProportion() + rootDes[1]->getBranchProportion();
        rootDes[0]->setBranchProportion(x/4.0);
        rootDes[1]->setBranchProportion(x/4.0);

        // rescale so branch proportions sum to 1.0
        double sum = 0.0;
        for (int i=0; i<downPassSequence.size(); i++)
            {
            p = downPassSequence[i];
            sum += p->getBranchProportion();
            }
        for (int i=0; i<downPassSequence.size(); i++)
            {
            p = downPassSequence[i];
            p->setBranchProportion( p->getBranchProportion()/sum );
            }
            
        treeLength = rv->gammaRv(1.0, betaT);
        }
    else
        {
        // rescale all branch proportions so they sum to one
        treeLength = 0.0;
        for (int i=0; i<nodes.size(); i++)
            treeLength += nodes[i]->getBranchProportion();
        for (int i=0; i<nodes.size(); i++)
            nodes[i]->setBranchProportion( nodes[i]->getBranchProportion()/treeLength );
        }
        
    //print();
}

Tree::~Tree(void) {

    for (int i=0; i<nodes.size(); i++)
        delete nodes[i];
}

Tree& Tree::operator=(Tree& t) {

    if (this != &t)
        clone(t);
    return *this;
}

Node* Tree::addNode(void) {

    Node* newNode = new Node( (int)nodes.size() );
    newNode->setMyTree(this);
    nodes.push_back(newNode);
    return newNode;
}

void Tree::clone(Tree& t) {

    // copy some instance variables
    numTaxa = t.numTaxa;
    taxonNames = t.taxonNames;
    treeLength = t.treeLength;
    
    // make certain we have the saame number of nodes in each tree
    if (nodes.size() != t.nodes.size())
        {
        deleteAllNodes();
        for (int i=0; i<t.nodes.size(); i++)
            addNode();
        }

    // copy the nodes
    root = nodes[t.root->getOffset()];
    for (int i=0; i<t.nodes.size(); i++)
        {
        Node* pLft = nodes[i];
        Node* pRht = t.nodes[i];
        
        pLft->setIndex( pRht->getIndex() );
        pLft->setIsLeaf( pRht->getIsLeaf() );
        pLft->setName( pRht->getName() );
        pLft->setBranchProportion( pRht->getBranchProportion() );
        
        if (pRht->getAncestor() != NULL)
            pLft->setAncestor( nodes[pRht->getAncestor()->getOffset()] );
        else
            pLft->setAncestor(NULL);
            
        pLft->removeDescendants();
        std::set<Node*,CompNode>& rhtNeighbors = pRht->getDescendants()->getNodes();
        for (Node* n : rhtNeighbors)
            pLft->addDescendant( nodes[n->getOffset()] );
        }
            
    // copy the down pass sequence
    downPassSequence.clear();
    for (int i=0; i<t.downPassSequence.size(); i++)
        downPassSequence.push_back( nodes[t.downPassSequence[i]->getOffset()] );
}

void Tree::deleteAllNodes(void) {

    for (int i=0; i<nodes.size(); i++)
        delete nodes[i];
    nodes.clear();
}

void Tree::initializeDownPassSequence(void) {

    downPassSequence.clear();
    passDown(root);
    int idx = numTaxa;
    for (int i=0; i<downPassSequence.size(); i++)
        {
        Node* p = downPassSequence[i];
        if (p->getIsLeaf() == false)
            p->setIndex(idx++);
        }
}

std::vector<int> Tree::getAncestorIndices(void) {

    std::vector<int> ancIndices(nodes.size());
    
    for (int i=0; i<downPassSequence.size(); i++)
        {
        Node* p = downPassSequence[i];
        int idx = p->getIndex();
        if (p == root)
            ancIndices[idx] = -1;
        else
            ancIndices[idx] = p->getAncestor()->getIndex();
        }
    
    return ancIndices;
}

std::vector<double> Tree::getBranchLengthVector(void) {

    std::vector<double> brlenVec(nodes.size());

    for (int i=0; i<downPassSequence.size(); i++)
        {
        Node* p = downPassSequence[i];
        int idx = p->getIndex();
        if (p == root)
            brlenVec[idx] = 0.0;
        else
            brlenVec[idx] = p->getBranchLength();
        }

    return brlenVec;
}

std::string Tree::getNewick(void) {

    std::stringstream ss;
    if (root->getIsLeaf() == true)
        {
        Node* oldRoot = root;
        std::vector<Node*> nbs = root->getDescendantsVector();
        if (nbs.size() > 1)
            Msg::error("Expecting only a single neighbor at the root of the tree");
        Node* newRoot = nbs[0];
        root->setAncestor(newRoot);
        oldRoot->setAncestor(newRoot);
        newRoot->setAncestor(NULL);
        root = newRoot;

        writeTree(root, ss);

        newRoot->setAncestor(oldRoot);
        oldRoot->setAncestor(NULL);
        root = oldRoot;
        }
    else
        {
        writeTree(root, ss);
        }
    std::string newick = ss.str();
    return newick;
}

bool Tree::isBinary(void) {

    for (int i=0; i<downPassSequence.size(); i++)
        {
        Node* p = downPassSequence[i];
        int ndes = (int)p->numDescendants();
        if (p->getIsLeaf() == true)
            {
            if (ndes != 0)
                return false;
            }
        else
            {
            if (ndes != 2)
                return false;
            }
        }
    return true;
}

bool Tree::isTaxonPresent(std::string tn) {

    for (int i=0; i<taxonNames.size(); i++)
        {
        if (tn == taxonNames[i])
            return true;
        }
    return false;
}

void Tree::listNodes(Node* p, size_t indent) {

    if (p != NULL)
        {
        std::set<Node*,CompNode> des = p->getDescendants()->getNodes();
        
        for (int i=0; i<indent; i++)
            std::cout << " ";
//        std::cout << p->getIndex() << " " << p << " ( ";
        std::cout << p->getIndex() << " ( ";
        if (p->getAncestor() != NULL)
            std::cout << "a_" << p->getAncestor()->getIndex() << " ";
        else
        std::cout << "a_NULL ";
        for (Node* n : des)
            {
            std::cout << n->getIndex() << " ";
            }
        std::cout << ") " << std::fixed << std::setprecision(5) << p->getBranchLength() << " ";
    
        if (p->getIsLeaf() == true)
            std::cout << " (" << p->getName() << ")";
            
        if (p == root)
            std::cout << " <- Root";
            
        std::cout << std::endl;

        for (std::set<Node*>::iterator it = des.begin(); it != des.end(); it++)
            {
            listNodes( (*it), indent+3 );
            }
        }
}

void Tree::print(void) {

    listNodes(root, 3);
}

void Tree::print(std::string header) {

    std::cout << header << std::endl;
    print();
}

void Tree::passDown(Node* p) {

    if (p != NULL)
        {
        std::set<Node*,CompNode>& des = p->getDescendants()->getNodes();
        for (std::set<Node*>::iterator it = des.begin(); it != des.end(); it++)
            {
            passDown( (*it) );
            }
        downPassSequence.push_back(p);
        }
}

Node* Tree::randomNode(RandomVariable* rv) {

    return nodes[(int)(rv->uniformRv()*nodes.size())];
}

std::vector<std::string> Tree::tokenizeTreeString(std::string ls) {

    std::vector<std::string> tks;
    
    std::string longToken = "";
    for (int i=0; i<ls.size(); i++)
        {
        char x = ls.at(i);
        if (x ==')' || x == '(' || x == ',' || x == ';' || x == ':')
            {
            if (longToken != "")
                {
                tks.push_back(longToken);
                longToken = "";
                }
            std::string s(1, x);
            tks.push_back( s );
            }
        else
            {
            longToken += x;
            }
        }
#   if 0
    for (std::string s : tks)
        {
        std::cout << s << std::endl;
        }
#   endif
    
    return tks;
}

void Tree::writeTree(Node* p, std::stringstream& ss) {

    if (p != NULL)
        {
        if (p->getIsLeaf() == true)
            {
            ss << p->getIndex()+1;
            ss << ":" << std::fixed << std::setprecision(5) << p->getBranchLength();
            }
        else
            {
            ss << "(";
            }
        std::vector<Node*> myDescendants = p->getDescendantsVector();
        for (int i=0; i<(int)myDescendants.size(); i++)
            {
            writeTree(myDescendants[i], ss);
            if ( (i + 1) != (int)myDescendants.size() )
                ss << ",";
            }
        if (p->getIsLeaf() == false)
            {
            ss << ")";
            if (p != NULL && p != root)
                ss << ":" << std::fixed << std::setprecision(5) << p->getBranchLength();
            }
        }
}


