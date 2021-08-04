#ifndef Tree_H
#define Tree_H

#include <sstream>
#include <string>
#include <vector>
class Node;
class RandomVariable;



class Tree {

    public:
                                    Tree(std::string fileName, std::vector<std::string> tNames, double betaT, RandomVariable* rv);
                                    Tree(Tree& t);
                                   ~Tree(void);
        Tree&                       operator=(Tree& t);
        std::string                 getNewick(void);
        int                         getNumNodes(void) { return (int)nodes.size(); }
        int                         getNumTaxa(void) { return numTaxa; }
        std::vector<int>            getAncestorIndices(void);
        std::vector<double>         getBranchLengthVector(void);
        std::vector<Node*>&         getDownPassSequence(void) { return downPassSequence; }
        Node*                       getRoot(void) { return root; }
        std::vector<std::string>&   getTaxonNames(void) { return taxonNames; }
        double                      getTreeLength(void) { return treeLength; }
        bool                        isBinary(void);
        bool                        isRoot(Node* p) { return ((p == root) ? true : false); }
        bool                        isTaxonPresent(std::string tn);
        void                        print(void);
        void                        print(std::string header);
        void                        setTreeLength(double x) { treeLength = x; }
                
    private:
                                    Tree(void) {}
        Node*                       addNode(void);
        void                        clone(Tree& t);
        void                        deleteAllNodes(void);
        void                        initializeDownPassSequence(void);
        void                        listNodes(Node* p, size_t indent);
        void                        passDown(Node* p);
        std::vector<std::string>    tokenizeTreeString(std::string ls);
        void                        writeTree(Node* p, std::stringstream& ss);
        std::vector<Node*>          nodes;
        std::vector<Node*>          downPassSequence;
        Node*                       root;
        int                         numTaxa;
        std::vector<std::string>    taxonNames;
        double                      treeLength;
};

#endif
