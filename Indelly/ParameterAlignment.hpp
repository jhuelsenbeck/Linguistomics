#ifndef ParameterAlignment_hpp
#define ParameterAlignment_hpp

#include <string>
#include <vector>
#include "Parameter.hpp"
class Alignment;
class Model;



class ParameterAlignment : public Parameter {

    public:
                                        ParameterAlignment(void) = delete;
                                        ParameterAlignment(const ParameterAlignment& pa) = delete;
                                        ParameterAlignment(RandomVariable* r, Model* m, Alignment* a, std::string n);
                                       ~ParameterAlignment(void);
        double                          update(void);
        void                            accept(void);
        std::vector<std::vector<int> >& getAlignment(void) { return alignment[0]; }
        std::vector<std::vector<int> >& getAlignment(int idx) { return alignment[idx]; }
        char                            getCharFromCode(int code);
        int                             getGapCode(void) { return gapCode; }
        std::string                     getHeader(void) { return ""; }
        std::vector<std::vector<int> >  getIndelMatrix(void);
        std::vector<std::vector<int> >  getIndelMatrix(int idx);
        std::vector<std::vector<int> >  getIndelMatrix(std::vector<std::vector<int> >& aln);
        int                             getNumSites(void) { return (int)alignment[0][0].size(); }
        int                             getNumStates(void) { return numStates; }
        int                             getPrintWidth(void) { return printWidth; }
        std::vector<std::vector<int> >  getRawSequenceMatrix(void) { return sequences; }
        std::string                     getString(void) { return ""; }
        std::vector<std::string>        getTaxonNames(void) { return taxonNames; }
        double                          lnPriorProbability(void);
        void                            print(void);
        void                            reject(void);
    
    protected:
        std::vector<std::vector<int> >  alignment[2];  // numTaxa X numSites
        std::vector<std::vector<int> >  sequences;     // numTaxa X numSites for taxon i (i.e., left-justified)
        std::vector<std::string>        taxonNames;
        double                          tuning;
        double                          exponent;
        double                          gapPenalty;
        int                             numStates;
        int                             gapCode;       // the gap code is simply the maximum number of states possible in any Markov model
        int                             printWidth;
};

#endif
