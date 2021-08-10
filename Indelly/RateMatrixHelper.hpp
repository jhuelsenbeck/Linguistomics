#ifndef RateMatrixHelper_hpp
#define RateMatrixHelper_hpp

#include <set>
#include <string>
#include <vector>
#include "json.hpp"

struct GroupPair {

                                    GroupPair(int i, int j);
        int                         group1;
        int                         group2;
};

struct CompGroupPair {

    bool operator()(const GroupPair& p1, const GroupPair& p2) const {
        
        if (p1.group1 < p2.group1)
            return true;
        else if (p1.group1 == p2.group1)
            {
            if (p1.group2 < p2.group2)
                return true;
            }
        return false;
        }
};

class RateMatrixHelper {

    public:
        static RateMatrixHelper&                rateMatrixHelper(void) {
                                                    static RateMatrixHelper rmh;
                                                    return rmh;
                                                }
        std::set<int>&                          getGroup(int idx) { return stateGroupings[idx]; }
        std::map<GroupPair,int, CompGroupPair>& getGroupIndices(void) { return groupIndices; }
        std::vector<std::string>                getLabels(void);
        int**                                   getMap(void) { return m; }
        int                                     getNumRates(void) { return (int)groupIndices.size(); }
        int                                     getNumStates(void) { return numStates; }
        int                                     groupIdForState(int stateIdx);
        void                                    initialize(int ns, nlohmann::json& j);
        void                                    print(void);
        
    private:
                                                RateMatrixHelper(void);
                                               ~RateMatrixHelper(void);
        int**                                   m;
        int                                     numStates;
        bool                                    isInitialized;
        std::vector<std::set<int> >             stateGroupings;
        std::vector<std::string>                stateGroupingsNames;
        int                                     numGroups;
        std::map<GroupPair,int, CompGroupPair>  groupIndices;
};

#endif
