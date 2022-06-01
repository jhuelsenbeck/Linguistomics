#ifndef AlignmentDistribution_hpp
#define AlignmentDistribution_hpp

#include <map>
#include <string>
#include "Alignment.hpp"
#include "json.hpp"
class Partition;
class RandomVariable;


class AlignmentDistribution {

    public:
                                                AlignmentDistribution(void) = delete;
                                                AlignmentDistribution(RandomVariable* r, Partition* p);
                                               ~AlignmentDistribution(void);
        void                                    addAlignment(Alignment* aln);
        int                                     ciSize(void);
        std::string                             getName(void) { return name; }
        Alignment*                              getMapAlignment(void);
        int                                     numSamples(void);
        void                                    print(void);
        void                                    setName(std::string s) { name = s; }
        size_t                                  size(void) { return samples.size(); }
        nlohmann::json                          toJson(double credibleSetSize);
            
    private:
        void                                    print(Alignment* aln);
        RandomVariable*                         rv;
        Partition*                              partition;
        std::string                             name;
        std::map<Alignment*,int>                samples;
};

#endif
