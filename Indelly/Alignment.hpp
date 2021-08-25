#ifndef Alignment_H
#define Alignment_H

#include <map>
#include <string>
#include <vector>
#include "json.hpp"


class Alignment {

    public:
                                        Alignment(void) = delete;
                                        Alignment(nlohmann::json& j, int ns, std::vector<std::string> canonicalTaxonList);
                                       ~Alignment(void);
        int                             getCharacter(size_t i, size_t j);
        int                             getGapCode(void) { return gapCode; }
        std::string                     getName(void) { return name; }
        int                             getMaximumNumberOfStates(void) { return numStates; }
        int                             getNumTaxa(void) { return numTaxa; }
        int                             getNumChar(void) { return numChar; }
        int                             getNumStates(void) { return numStates; }
        std::vector<int>                getRawSequence(int txnIdx);
        std::vector<std::vector<int> >  getRawSequenceMatrix(void);
        std::vector<std::vector<int> >  getIndelMatrix(void);
        int                             getTaxonIndex(std::string ns);
        std::map<int,int>               getTaxonMap(void) { return taxonMap; }
        std::vector<bool>               getTaxonMask(void) { return taxonMask; }
        std::vector<std::string>        getTaxonNames(void);
        std::string                     getTaxonName(int i);
        bool                            hasAllGapColumn(void);
        bool                            isAlignmentComplete(void);
        bool                            isIndel(size_t i, size_t j);
        void                            listTaxa(void);
        int                             numCompleteTaxa(void);
        void                            print(void);
        void                            printIndels(void);
        void                            setName(std::string s) { name = s; }

    private:
        std::string                     bomLessString(std::string& str);
        bool                            hasBOM(std::string& str);
        bool                            isInteger(const std::string& str);
        std::string                     name;
        int**                           matrix;
        bool**                          indelMatrix;
        int                             numTaxa;
        int                             numChar;
        int                             numStates;
        int                             gapCode;
        std::vector<std::string>        taxonNames;
        std::vector<bool>               taxonMask;
        std::map<int,int>               taxonMap;
};

#endif
