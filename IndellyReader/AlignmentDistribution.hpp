#ifndef AlignmentDistribution_hpp
#define AlignmentDistribution_hpp

#include <map>
#include <string>
#include "Alignment.hpp"
class Partition;
class RandomVariable;


class AlignmentDistribution {

    public:
                                                AlignmentDistribution(void) = delete;
                                                AlignmentDistribution(RandomVariable* r, Partition* p);
                                               ~AlignmentDistribution(void);
        void                                    addAlignment(Alignment* aln);
        std::string                             getName(void) { return name; }
        int                                     numSamples(void);
        void                                    print(void);
        void                                    setName(std::string s) { name = s; }
        size_t                                  size(void) { return samples.size(); }
            
    private:
        void                                    print(Alignment* aln);
        RandomVariable*                         rv;
        Partition*                              partition;
        std::string                             name;
        std::map<Alignment*,int>                samples;
};

#endif
