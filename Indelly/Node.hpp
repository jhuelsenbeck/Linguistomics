#ifndef Node_H
#define Node_H

#include <set>
#include <string>
#include <vector>
class NodeSet;
class Tree;



class Node {

    public:
                            Node(void) = delete;
                            Node(int idx);
                           ~Node(void);
        void                addDescendant(Node* p);
        void                clean(void);
        Node*               getAncestor(void) { return ancestor; }
        double              getBranchLength(void);
        double              getBranchProportion(void) { return proportion; }
        NodeSet*&           getDescendants(void) { return descendants; }
        std::vector<Node*>  getDescendantsVector(void);
        bool                getFlag(void) { return flag; }
        int                 getIndex(void) { return index; }
        bool                getIsLeaf(void) { return isLeaf; }
        std::string         getName(void) { return name; }
        int                 getOffset(void) { return offset; }
        int                 getOffset(void) const { return offset; }
        Node*               getSisterNode(void);
        bool                isDescendant(Node* p);
        size_t              numDescendants(void);
        void                print(void);
        void                removeDescendant(Node* p);
        void                removeDescendants(void);
        void                setAncestor(Node* p) { ancestor = p; }
        void                setBranchProportion(double x) { proportion = x; }
        void                setFlag(bool tf) { flag = tf; }
        void                setIndex(int x) { index = x; }
        void                setIsLeaf(bool tf) { isLeaf = tf; }
        void                setMyTree(Tree* t) { myTree = t; }
        void                setName(std::string s) { name = s; }
        void                setOffset(int x) { offset = x; }

    protected:
        NodeSet*            descendants;
        Node*               ancestor;
        int                 index;
        bool                isLeaf;
        std::string         name;
        double              proportion;
        int                 offset;
        Tree*               myTree;
        bool                flag;
};

struct CompNode {

    bool operator()(const Node* n1, const Node* n2) const {
        
        if (n1->getOffset() > n2->getOffset())
            return true;
        return false;
        }
};

#endif
