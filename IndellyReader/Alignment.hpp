#ifndef Alignment_hpp
#define Alignment_hpp

#include "json.hpp"
#include "Sequence.hpp"
#include <string>
#include <vector>


class Alignment {

    public:
                                Alignment(void) = delete;
                                Alignment(nlohmann::json j);
        Sequence&               operator[](size_t i) { return matrix[i]; }
        const Sequence&         operator[](size_t i) const { return matrix[i]; };
        bool                    operator==(const Alignment& aln) const;
        int                     getCharCode(int tIdx, int cIdx) { return matrix[tIdx][cIdx]; }
        int                     getNumChar(void) { return numChar; }
        int                     getNumTaxa(void) { return numTaxa; }
        int                     lengthOfLongestName(void);
        void                    print(void);
        void                    print(std::string h);
        size_t                  size(void) { return matrix.size(); }
        size_t                  size(void) const { return matrix.size(); }
    
    private:
        int                     numTaxa;
        int                     numChar;
        std::vector<Sequence>   matrix;
};

struct CompAlignment {

    bool operator()(const Alignment& a1, const Alignment& a2) const {
        
        if (a1.size() > a2.size())
            return true;
        else if ( a1.size() == a2.size())
            {
            if (a1[0].size() > a2[0].size())
                return true;
            else if ( a1[0].size() == a2[0].size() )
                {
                for (int i=0; i<a1.size(); i++)
                    {
                    for (int j=0; j<a1[i].size(); j++)
                        {
                        if (a1[i][j] > a2[i][j])
                            return true;
                        }
                    }
                }
            }
        return false;
        }

    bool operator()(const Alignment* a1, const Alignment* a2) const {

        if (a1->size() > a2->size())
            return true;
        else if ( a1->size() == a2->size())
            {
            if ((*a1)[0].size() > (*a2)[0].size())
                return true;
            else if ( (*a1)[0].size() == (*a2)[0].size() )
                {
                for (int i=0; i<a1->size(); i++)
                    {
                    for (int j=0; j<(*a1)[i].size(); j++)
                        {
                        if ((*a1)[i][j] > (*a2)[i][j])
                            return true;
                        }
                    }
                }
            }
        return false;
        }
};

#endif
