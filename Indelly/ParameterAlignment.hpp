#ifndef ParameterAlignment_hpp
#define ParameterAlignment_hpp

#include <map>
#include <string>
#include <vector>
#include "Parameter.hpp"
#include "RbBitSet.h"
class Alignment;
class AlignmentProposal;
class Model;
class SiteLikelihood;



class ParameterAlignment : public Parameter {

    public:
                                        ParameterAlignment(void) = delete;
                                        ParameterAlignment(const ParameterAlignment& pa) = delete;
                                        ParameterAlignment(RandomVariable* r, Model* m, Alignment* a, std::string n, SiteLikelihood* sl, int idx);
                                       ~ParameterAlignment(void);
        void                            accept(void);
        bool                            areAlignmentsIdentical(void);
        std::vector<std::vector<int> >& getAlignment(void) { return alignment[0]; }
        std::vector<std::vector<int> >& getAlignment(int idx) { return alignment[idx]; }
        char                            getCharFromCode(int code);
        int                             getGapCode(void) { return gapCode; }
        std::string                     getHeader(void) { return ""; }
        std::vector<std::vector<int> >  getIndelMatrix(void);
        std::vector<std::vector<int> >  getIndelMatrix(int idx);
        std::vector<std::vector<int> >  getIndelMatrix(std::vector<std::vector<int> >& aln);
        int                             getIndex(void) { return index; }
        bool                            getIsCompletelySampled(void) { return completelySampled; }
        std::string                     getJsonString(void);
        int                             getNumSites(void) { return (int)alignment[0][0].size(); }
        int                             getNumStates(void) { return numStates; }
        int                             getNumTaxa(void) { return numTaxa; }
        int                             getPrintWidth(void) { return printWidth; }
        std::vector<std::vector<int> >  getRawSequenceMatrix(void) { return sequences; }
        SiteLikelihood*                 getSiteProbs(void) { return siteProbs; }
        std::string                     getString(void) { return ""; }
        RbBitSet                        getTaxonMask(void);
        std::string                     getTaxonMaskString(void);
        std::vector<std::string>        getTaxonNames(void) { return taxonNames; }
        double                          lnPriorProbability(void);
        void                            print(void);
        void                            reject(void);
        double                          update(void);
    
    protected:
        AlignmentProposal*              alignmentProposal;
        std::vector<std::vector<int> >  alignment[2];  // numTaxa X numSites
        std::vector<std::vector<int> >  sequences;     // numTaxa X numSites for taxon i (i.e., left-justified)
        bool                            completelySampled;
        double                          tuning;
        double                          exponent;
        double                          gapPenalty;
        int                             numStates;
        int                             gapCode;       // the gap code is simply the maximum number of states possible in any Markov model
        int                             numTaxa;
        int                             printWidth;
        int                             index;
        SiteLikelihood*                 siteProbs;
        std::vector<std::string>        taxonNames;
        std::vector<bool>               taxonMask;
        std::map<int,int>               taxonMapKeyCanonical;
        std::map<int,int>               taxonMapKeyAlignment;
};

#endif
