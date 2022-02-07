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
        RandomVariable*                     rv;
        ParameterAlignment*                 alignmentParm;
        Model*                              model;
        
        RbBitSet                            taxonMask;
        Tree*                               tree;
        int                                 numTaxa;
        int                                 numNodes;

        double                              gap;
        int                                 gapCode;
        double                              basis;
        int                                 maxlength;
        int                                 maxUnalignDimension;
        enum                                StateLabels { free, possible, edgeUsed, used };
        static int                          bigUnalignableRegion;
        std::vector<IntVector*>             pool;
        std::set<IntVector*>                allocated;
        
        int***                              profile;
        std::vector<int>                    possibles;      // resized once on instantiation, never copied
        std::vector<int>                    state;          // resized to large value on instantiation, resized frequently to smaller values, no realloc should occur
        std::vector<std::vector<int> >      alignment;
        std::vector<std::vector<int> >      profileNumber;
        int*                                xProfile;
        int*                                yProfile;
        int                                 numStates;
        double**                            scoring;
        std::vector<std::vector<int> >      tempProfile;
        int**                               dp;
};
#endif
