#ifndef AlignmentProposal_hpp
#define AlignmentProposal_hpp

#include <map>
#include <vector>
#include "IntVector.hpp"
class ParameterAlignment;
class RandomVariable;
class Tree;



class AlignmentProposal {

    public:
                            AlignmentProposal(void) = delete;
                            AlignmentProposal(ParameterAlignment* a, Tree* t, RandomVariable* r, double expnt, double gp);
        double              propose(std::vector<std::vector<int> >& newAlignment, double iP);
                
    private:
        void                cleanTable(std::map<IntVector*,int,CompIntVector>& m);
        int                 countPaths(std::vector<std::vector<int> >& inputAlignment, int startCol, int endCol);
        void                freeProfile(int*** x, int n);
        void                print(std::string s, std::vector<int>& x);
        void                print(std::string s, std::vector<std::vector<int> >& x);
        void                print(std::string s, int*** x, int a, int b, int c);
        double              gap;
        int                 gapCode;
        double              basis;
        int                 maxlength;
        int                 maxUnalignDimension;
        enum                StateLabels { free, possible, edgeUsed, used };
        Tree*               tree;
        RandomVariable*     rv;
        ParameterAlignment* alignmentParm;
        static int          bigUnalignableRegion;
};
#endif
