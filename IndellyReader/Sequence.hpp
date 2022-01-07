#ifndef Sequence_hpp
#define Sequence_hpp

#include <string>
#include <vector>


class Sequence {

    public:
                            Sequence(void) = delete;
                            Sequence(std::string n, std::vector<int> d);
        int&                operator[](size_t i) { return data[i]; }
        const int&          operator[](size_t i) const { return data[i]; };
        std::string         getName(void) { return name; }
        std::vector<int>    getSequence(void) { return data; }
        size_t              size(void) { return data.size(); }
        size_t              size(void) const { return data.size(); }
        
    private:
        std::string         name;
        std::vector<int>    data;
};

#endif
