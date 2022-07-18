#ifndef McmcSummary_hpp
#define McmcSummary_hpp

#include <map>
#include <string>
#include <vector>
#include "ParameterStatistics.hpp"
#include "UserSettings.hpp"
#include "RbBitSet.h"

class AlignmentDistribution;
class Partition;
class RandomVariable;
class Tree;



class McmcSummary {

    public:
                                                McmcSummary(void) = delete;
                                                McmcSummary(RandomVariable* r);
                                               ~McmcSummary(void);
        std::vector<CredibleInterval>           getCredibleIntervals(void);
        std::vector<double>                     getMeans(void);
        void                                    print(void);
        void                                    printPartitionSet(void);
        void                                    output(UserSettings& settings);
        void                                    readAlnFile(std::string fn, int bi);
        void                                    readTreFile(std::string fn, int bi);
        void                                    readTsvFile(std::string fn, int bi);
        void                                    readConfigFile(std::string fn);
    
    private:
        void                                    addPartion(std::map<RbBitSet,double>& parts);
        std::vector<std::string>                breakString(std::string str);
        std::string                             getCognateName(std::string str);
        bool                                    hasSemicolon(std::string str);
        std::map<int,std::string>               interpretTranslateString(std::vector<std::string> translateTokens);
        std::string                             interpretTreeString(std::string str);
        int                                     parseNumberFromFreqHeader(std::string str);
        void                                    printPartitionFreqs(void);
        RandomVariable*                         rv;
        std::vector<ParameterStatistics*>       stats;
        std::vector<AlignmentDistribution*>     alignments;
        std::map<RbBitSet,ParameterStatistics*> partitions;
        Tree*                                   conTree;
        bool                                    hasPartitions;
        Partition*                              statePartitions;
        std::vector<std::string>                taxa;
};

#endif
