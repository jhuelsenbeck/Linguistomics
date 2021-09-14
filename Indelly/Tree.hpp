#ifndef Tree_H
#define Tree_H

#include <map>
#include "RbBitSet.h"
#include <sstream>
#include <string>
#include <vector>
class Node;
class RandomVariable;



class Tree {

    public:
                                    Tree(std::vector<std::string> tNames, double betaT, RandomVariable* rv);
                                    Tree(std::string treeStr, std::vector<std::string> tNames, double betaT, RandomVariable* rv);
                                    Tree(Tree& t);
                                    Tree(Tree& t, std::vector<bool> taxonMask, std::vector<std::string> tn);
                                    Tree(RbBitSet& taxonMask, std::map<RbBitSet*,double,CompBitSet>& partitions, std::vector<std::string>& tn);
                                   ~Tree(void);
        Tree&                       operator=(Tree& t);
        void                        buildRandomTree(std::vector<std::string> tNames, double betaT, RandomVariable* rv);
        void                        debugPrint(std::string h);
        std::string                 getNewick(void);
        int                         getNumNodes(void) { return (int)nodes.size(); }
        int                         getNumTaxa(void) { return numTaxa; }
        std::vector<int>            getAncestorIndices(void);
        std::vector<double>         getBranchLengthVector(void);
        std::vector<Node*>&         getDownPassSequence(void) { return downPassSequence; }
        std::map<RbBitSet*,double,CompBitSet>  getTaxonBipartitions(void);
        Node*                       getRoot(void) { return root; }
        std::vector<std::string>&   getTaxonNames(void) { return taxonNames; }
        int                         getTaxonNameIndex(std::string tName);
        double                      getTreeLength(void) { return treeLength; }
        void                        initializeDownPassSequence(void);
        bool                        isBinary(void);
        bool                        isRoot(Node* p) { return ((p == root) ? true : false); }
        bool                        isTaxonPresent(std::string tn);
        void                        print(void);
        void                        print(std::string header);
        Node*                       randomNode(RandomVariable* rv);
        void                        setTreeLength(double x) { treeLength = x; }
                
    private:
                                    Tree(void) {}
        Node*                       addNode(void);
        void                        clone(Tree& t);
        void                        deleteAllNodes(void);
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
