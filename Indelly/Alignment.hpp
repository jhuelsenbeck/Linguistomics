#ifndef Alignment_H
#define Alignment_H

#include <string>
#include <vector>


class Alignment {

    public:
                                        Alignment(void) = delete;
                                        Alignment(std::string fileName, bool dirMode);
                                       ~Alignment(void);
        int                             getCharacter(size_t i, size_t j);
        char                            getCharFromCode(int code);
        std::string                     getDataType(void) { return dataType; }
        int                             getGapCode(void) { return gapCode; }
        std::string                     getName(void) { return name; }
        int                             getMaximumNumberOfStates(void) { return (int)states.length(); }
        int                             getNumTaxa(void) { return numTaxa; }
        int                             getNumSites(void) { return numSites; }
        int                             getNumStates(void) { return numStates; }
        std::vector<int>                getRawSequence(int txnIdx);
        std::vector<std::vector<int> >  getRawSequenceMatrix(void);
        std::vector<std::vector<int> >  getIndelMatrix(void);
        std::string                     getStates(void) { return states; }
        int                             getTaxonIndex(std::string ns);
        std::vector<std::string>        getTaxonNames(void);
        std::string                     getTaxonName(int i);
        bool                            isIndel(size_t i, size_t j);
        void                            listTaxa(void);
        void                            print(void);
        void                            printIndels(void);
        void                            setName(std::string s) { name = s; }
        int                             stateCode(char s);

    private:
        bool                            isInteger(const std::string& str);
        std::vector<std::string>        taxonNames;
        int**                           matrix;
        bool**                          indelMatrix;
        int                             numTaxa;
        int                             numSites;
        int                             numStates;
        std::string                     dataType;
        std::string                     name;
        static std::string              states;
        static int                      gapCode;
};

#endif
