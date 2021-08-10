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
        void            run(void);
    
    private:
        void            closeOutputFiles(void);
        int             numDigits(double lnX);
        void            openOutputFiles(void);
        void            print(int gen, double curLnL, double newLnL, double curLnP, double newLnP, bool accept);
        double          safeExponentiation(double lnX);
        void            sample(int gen, double lnL, double lnP);
        std::string     formattedTime(std::chrono::high_resolution_clock::time_point& t1, std::chrono::high_resolution_clock::time_point& t2);
        void            runPathSampling(void);
        void            runPosterior(void);
        Model*          modelPtr;
        RandomVariable* rv;
        int             numMcmcCycles;
        int             printFrequency;
        int             sampleFrequency;
        std::ofstream   algnStrm;
        std::ofstream   parmStrm;
        std::ofstream   treeStrm;
        int             maxPriorPrint;
        int             maxLikePrint;
        int             maxGenPrint;
};

#endif
