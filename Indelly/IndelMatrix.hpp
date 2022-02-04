#ifndef IndelMatrix_hpp
#define IndelMatrix_hpp


class IndelMatrix {

    public:
                    IndelMatrix(void) = delete;
                    IndelMatrix(int nt, int maxNs);
                   ~IndelMatrix(void);
        int*        operator[](int i) { return m[i]; }
        const int*  operator[](int i) const { return m[i]; }
        int*        getRow(int idx) { return m[idx]; }
        int         getNumTaxa(void) { return numTaxa; }
        int         getNumSites(void) { return numSites; }
        void        print(void);
        void        setNumSites(int ns);
    
    private:
        int         numTaxa;
        int         numSites;
        int         maxNumSites;
        int**       m;
};

#endif
