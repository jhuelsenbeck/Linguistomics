#include <iostream>
#include "Node.hpp"
#include "NodeSet.hpp"
#include "Tree.hpp"



Node::Node(int idx) {

    index      = 0;
    isLeaf     = false;
    name       = "";
    ancestor   = NULL;
    proportion = 0.0;
    offset     = idx;
    myTree     = NULL;
    descendants = new NodeSet;
}

Node::~Node(void) {

    delete descendants;
}

void Node::addDescendant(Node* p) {

    descendants->addNode(p);
}

void Node::clean(void) {

    index      = 0;
    isLeaf     = false;
    name       = "";
    ancestor   = NULL;
    proportion = 0.0;
    myTree     = NULL;
    descendants->deleteAllNodes();
}

double Node::getBranchLength(void) {

    return proportion * myTree->getTreeLength();
}

std::vector<Node*> Node::getDescendantsVector(void) {

    std::set<Node*,CompNode>& des = descendants->getNodes();
    std::vector<Node*> nb;
    for (Node* n : des)
        nb.push_back( n );
    return nb;
}

bool Node::isDescendant(Node* p) {

    std::set<Node*,CompNode>& des = descendants->getNodes();
    std::set<Node*,CompNode>::iterator it = des.find(p);
    if (it != des.end())
        return true;
    return false;
}

size_t Node::numDescendants(void) {

    return descendants->size();
}

void Node::print(void) {

    std::cout << "Node (" << this << ")" << std::endl;
    std::cout << "    Ancestor: " << ancestor->getIndex() << std::endl;
    std::cout << "    Descendants: ";
    std::set<Node*,CompNode>& des = descendants->getNodes();
    for(Node* n : des)
        {
        std::cout << n->getIndex() << " ";
        }
    std::cout << std::endl;
    std::cout << "        index: " << index  << std::endl;
    std::cout << "       isLeaf: " << isLeaf << std::endl;
    std::cout << "         name: \"" << name << "\""  << std::endl;
}

void Node::removeDescendant(Node* p) {

    descendants->deleteNode(p);
}

void Node::removeDescendants(void) {

    descendants->deleteAllNodes();
}
