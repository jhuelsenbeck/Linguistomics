#include <iostream>
#include "IndelMatrix.hpp"


IndelMatrix::IndelMatrix(int nt, int maxNs) {

    numTaxa = nt;
    maxNumSites = maxNs;
    m = new int*[maxNumSites];
    m[0] = new int[numTaxa*maxNumSites];
    for (int i=1; i<maxNumSites; i++)
        m[i] = m[i-1] + numTaxa;
    for (int i=0; i<maxNumSites; i++)
        for (int j=0; j<numTaxa; j++)
            m[i][j] = 0;
}

IndelMatrix::~IndelMatrix(void) {

    delete [] m[0];
    delete [] m;
}

void IndelMatrix::print(void) {

    std::cout << this << std::endl;
    for (int i=0; i<numSites; i++)
        {
        for (int j=0; j<numTaxa; j++)
            std::cout << m[i][j];
        std::cout << std::endl;
        }
}

void IndelMatrix::setNumSites(int ns) {

    if (ns > maxNumSites)
        {
        std::cout << "Error: invalid attempt to increase the number of sites for an IndelMatrix";
        exit(1);
        }
    numSites = ns;
}
