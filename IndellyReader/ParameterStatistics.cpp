#include <algorithm>
#include "Msg.hpp"
#include "ParameterStatistics.hpp"



double ParameterStatistics::operator[](size_t idx) {

    if (idx >= values.size())
        Msg::error("Subscript index out-of-range");
    return values[idx];
}

double ParameterStatistics::operator[](size_t idx) const {

    if (idx >= values.size())
        Msg::error("Subscript index out-of-range");
    return values[idx];
}

CredibleInterval ParameterStatistics::getCredibleInterval(void) {

    size_t size = values.size();
    sortValues();
    size_t lp = size / 40;
    size_t up = size - lp - 1;

    CredibleInterval ci(values[lp], values[up], 0.0);
    if (size == 0)
        {
        ci.median = 0.0;
        }
    else
        {
        if (size % 2 == 0)
            ci.median = (values[size / 2 - 1] + values[size / 2]) / 2.0;
        else
            ci.median = values[size / 2];
        }
        
    return ci;
}

double ParameterStatistics::getMean(void) {

    double a = 0.0, s = 0.0;
    for (int i=0; i<values.size(); i++)
        {
        double x = values[i];
        if (i == 0)
            {
            a = x;
            s = 0.0;
            }
        else
            {
            double aOld = a;
            a = aOld + (x - aOld) / (i+1);
            s = s + (x - aOld) * (x - aOld);
            }
        }

    double mean = a;
    double variance = 0.0;
    if (values.size() > 1)
        variance = s / (values.size() - 1);
        
    return mean;
}

void ParameterStatistics::sortValues(void) {

    std::sort(values.begin(), values.end());
}

nlohmann::json ParameterStatistics::toJson(void) {

    double m = getMean();
    CredibleInterval i = getCredibleInterval();
    nlohmann::json j = nlohmann::json::object();
    j["cognate"] = name;
    j["mean"] = m;
    j["lower"] = i.lower;
    j["upper"] = i.upper;
    return j;
}
