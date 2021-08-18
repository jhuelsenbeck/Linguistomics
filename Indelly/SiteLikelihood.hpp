#ifndef SiteLikelihood_hpp
#define SiteLikelihood_hpp

#include <iostream>


class SiteLikelihood {

    public:
                    SiteLikelihood(void) = delete;
                    SiteLikelihood(int nn, int ns);
                   ~SiteLikelihood(void);
        double**    getProbsH(void) { return probsH; }
        double*     getProbsH(int idx) { return probsH[idx]; }
        double**    getProbsI(void) { return probsI; }
        double*     getProbsI(int idx) { return probsI[idx]; }
        void        print(void);
        void        zeroOutH(int idx) { memcpy( zeroH, probsH[idx], numStates*sizeof(double) ); }
        void        zeroOutI(int idx) { memcpy( zeroI, probsI[idx], (numStates+1)*sizeof(double) ); }
    
    private:
        int         numNodes;
        int         numStates;
        double**    probsH;
        double**    probsI;
        double*     zeroH;
        double*     zeroI;
};

#endif
