#ifndef ParameterStatistics_hpp
#define ParameterStatistics_hpp

#include "json.hpp"
#include <string>
#include <vector>

struct CredibleInterval {
    double lower;
    double upper;
    double median;

    CredibleInterval(double low, double up, double med) {
        lower = low;
        upper = up;
        median = med;
    }
};




class ParameterStatistics {

    public:
        double              operator[](size_t idx);
        double              operator[](size_t idx) const;
        void                addValue(double x) { values.push_back(x); }
        CredibleInterval    getCredibleInterval(void);
        double              getMean();
        std::string         getName(void) { return name; }
        int                 size(void) { return (int)values.size(); }
        void                setName(std::string s) { name = s; }
        void                sortValues(void);
        nlohmann::json      toJson();
        void                toFile(std::ostream& findex);

    private:
        std::string         name;
        std::vector<double> values;
};

#endif
