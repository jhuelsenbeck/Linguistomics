#ifndef AlignmentProposal_hpp
#define AlignmentProposal_hpp

#include <map>
#include <set>
#include <vector>
#include "IntVector.hpp"
#include "RbBitSet.h"
class Model;
class ParameterAlignment;
class RandomVariable;
class Tree;



class AlignmentProposal {

    public:
                                            AlignmentProposal(void) = delete;
                                            AlignmentProposal(ParameterAlignment* a, Tree* t, RandomVariable* r, Model* m, double expnt, double gp);
                                           ~AlignmentProposal(void);
        double                              propose(std::vector<std::vector<int> >& newAlignment, double iP);
                
    private:
        void                                cleanTable(std::map<IntVector*,int,CompIntVector>& m);
        int                                 countPaths(std::vector<std::vector<int> >& inputAlignment, int startCol, int endCol);
        void                                drainPool(void);
        void                                freeProfile(int*** x, int n);
        IntVector*                          getVector(void);
        IntVector*                          getVector(IntVector& vec);
        int                                 getNumAllocated(void) { return (int) allocated.size(); }
        int                                 getNumInPool(void) { return (int)pool.size(); }
        void                                initialize(void);
        void                                print(std::string s, std::vector<int>& x);
        void                                print(std::string s, std::vector<std::vector<int> >& x);
        void                                print(std::string s, int*** x, int a, int b, int c);
        void                                returnToPool(IntVector* n);
        Tree*                               tree;
        RandomVariable*                     rv;
        ParameterAlignment*                 alignmentParm;
        Model*                              model;
        double                              gap;
        int                                 gapCode;
        double                              basis;
        int                                 maxlength;
        int                                 maxUnalignDimension;
        enum                                StateLabels { free, possible, edgeUsed, used };
        static int                          bigUnalignableRegion;
        int                                 numTaxa;
        int                                 numNodes;
        RbBitSet                            taxonMask;
        std::vector<IntVector*>             pool;
        std::set<IntVector*>                allocated;
        
        int***                              profile;
        std::vector<int>                    possibles;
        std::vector<int>                    state;
        std::vector<std::vector<int> >      alignment;
        std::vector<std::vector<int> >      profileNumber;
        std::vector<int>                    xProfile;
        std::vector<int>                    yProfile;
        int                                 numStates;
        std::vector<std::vector<double> >   scoring;
        std::vector<std::vector<int> >      tempProfile;
        std::vector<std::vector<double> >   dp;
};
#endif
