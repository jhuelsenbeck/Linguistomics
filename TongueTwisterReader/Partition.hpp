#ifndef Partition_hpp
#define Partition_hpp

#include <map>
#include <set>
#include "json.hpp"
class Subset;


class Partition {

    public:
                            Partition(void) = delete;
                            Partition(nlohmann::json js);
                           ~Partition(void);
        bool                operator==(const Partition& rhs) const;
        void                addSubset(std::string s, std::vector<int> v);
        int                 getNumElements(void);
        std::set<Subset*>&  getSubsets(void) { return subsets; }
        Subset*             findSubsetIndexed(int x);
        Subset*             findSubsetIndexed(int x) const;
        Subset*             findSubsetLabeled(std::string sName);
        Subset*             findSubsetWithValue(int x);
        int                 indexOfSubsetWithValue(int x);
        bool                isEqualTo(const Partition& rhs) const;
        int                 numSubsets(void) { return (int)subsets.size(); }
        void                print(void);
        int                 maxValue(void);
        nlohmann::json      toJson(void);
        nlohmann::json      toJson(std::map<int,double>& partFreqs, std::ostream& findex);
    
    private:
        std::set<Subset*>   subsets;
};


#endif
