#ifndef Mcmc_hpp
#define Mcmc_hpp

#include <chrono>
#include <fstream>
#include <string>
class Model;
class RandomVariable;



class Mcmc {

    public:
                            Mcmc(void) = delete;
                            Mcmc(Model* m, RandomVariable* r);
        void                run(void);
    
    private:
        std::vector<double> calculatePowers(int numStones, double alpha, double beta);
        void                closeOutputFiles(void);
        std::string         formattedTime(std::chrono::high_resolution_clock::time_point& t1, std::chrono::high_resolution_clock::time_point& t2);
        void                initialize(void);
        int                 numDigits(double lnX);
        void                openOutputFiles(void);
        int                 phaseLength(std::string phs);
        void                print(int gen, double curLnL, double newLnL, double curLnP, double newLnP, bool accept, std::chrono::high_resolution_clock::time_point& t1, std::chrono::high_resolution_clock::time_point& t2);
        double              safeExponentiation(double lnX);
        void                sample(int gen, double lnL, double lnP);
        void                runPathSampling(void);
        void                runPosterior(void);
        Model*              modelPtr;
        RandomVariable*     rv;
        int                 numMcmcCycles;
        int                 printFrequency;
        int                 sampleFrequency;
        std::ofstream*      algnJsonStrm;
        std::ofstream       parmStrm;
        std::ofstream       treeStrm;
        double*             parmValues;
        int                 numParmValues;
        int                 maxPriorPrint;
        int                 maxLikePrint;
        int                 maxGenPrint;
        int                 preburninLength;
        int                 numTunes;
        int                 tuneLength;
        int                 burninLength;
        int                 sampleLength;
        int                 stoneSampleFrequency;
};

#endif
