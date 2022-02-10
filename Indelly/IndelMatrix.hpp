#ifndef IndelMatrix_hpp
#define IndelMatrix_hpp

#include "Container.hpp"



class IndelMatrix : public MatrixTemplate<int> {

    public:
                    IndelMatrix(void) = delete;
                    IndelMatrix(int nt, int maxNs);
        int*        getRow(int r) { return &(this->buffer[r * numCols]); }
        int         getNumTaxa(void) { return numTaxa; }
        int         getNumSites(void) { return numSites; }
        void        print(void);
        void        setNumSites(int ns);
    
    private:
        int         numTaxa;
        int         numSites;
        int         maxNumSites;
};

#endif
