#ifndef AlignmentDistribution_hpp
#define AlignmentDistribution_hpp

#include <map>
#include <string>
#include <iostream>
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
        nlohmann::json                          toJson(int index, double credibleSetSize, std::ostream& findex, std::ostream& fdata);

    private:
        void                                    print(Alignment* aln);
        RandomVariable*                         rv;
        Partition*                              partition;
        std::string                             name;
        std::map<Alignment*,int>                samples;
};

#endif
