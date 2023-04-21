#ifndef McmcSummary_hpp
#define McmcSummary_hpp

#include <map>
#include <string>
#include <vector>
#include "Container.hpp"
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
        McmcSummary&                            operator+=(const McmcSummary& rhs);
        void                                    calculateRates(DoubleMatrix& m);
        void                                    calculateAverageRates(DoubleMatrix& m);
        std::vector<CredibleInterval>           getCredibleIntervals(void);
        std::vector<double>                     getMeans(void);
        AlignmentDistribution*                  getAlignmentNamed(std::string str) const;
        ParameterStatistics*                    getParameterNamed(std::string str) const;
        void                                    print(void);
        void                                    printPartitionSet(void);

        void                                    output(std::string pathName, std::ofstream& findex);
        void                                    readAlnFile(std::string fn, int bi, Partition* prt);
        void                                    readTreFile(std::string fn, int bi);
        void                                    readTsvFile(std::string fn, int bi);
        void                                    readConfigFile(std::string fn);
        Partition*                              readPartition(std::string fn);
        int                                     readNumStates(std::string fn);
    
    private:
        void                                    writeMatrix(std::ofstream& file, DoubleMatrix& m, std::string name);
        void                                    addPartion(std::map<RbBitSet,double>& parts);
        std::vector<std::string>                breakString(std::string str);
        std::string                             getCognateName(std::string str);
        int                                     getFreqElement(std::string str);
        void                                    getRateElements(std::string str, int& r1, int& r2);
        bool                                    hasSemicolon(std::string str);
        int                                     inferNumberOfRates(void);
        int                                     inferNumberOfStates(void);
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
        int                                     numStates;
};

#endif
