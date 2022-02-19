#include <iomanip>
#include <iostream>
#include "Alignment.hpp"
#include "Sequence.hpp"


Alignment::Alignment(nlohmann::json js) {

    for (int i=0; i<js.size(); i++)
        {
        std::string tName = js[i]["Taxon"];
        std::vector<int> tInfo = js[i]["Segments"];
        matrix.push_back( Sequence(tName, tInfo) );
        }
    numTaxa = (int)matrix.size();
    numChar = (int)matrix[0].size();
}

bool Alignment::operator==(const Alignment& aln) const {

    if (matrix.size() != aln.matrix.size())
        return false;
    for (size_t i=0; i<matrix.size(); i++)
        {
        if (matrix[i].size() != aln.matrix[i].size())
            return false;
        for (int j=0; j<matrix[i].size(); j++)
            {
            if (matrix[i][j] != aln.matrix[i][j])
                return false;
            }
        }
    return true;
}

std::map<double,double> Alignment::gapInfo(void) {

    std::map<double,double> info;
    
    for (int j=0; j<numChar; j++)
        {
        int numGaps = 0;
        for (int i=0; i<numTaxa; i++)
            {
            int charCode = matrix[i][j];
            if (charCode == -1)
                numGaps++;
            }
            
        double fracGaps = (double)numGaps / numTaxa;
        double pos = (double)(j+1) / numChar;
        info.insert( std::make_pair(pos,fracGaps) );
        }
    
    return info;
}

int Alignment::lengthOfLongestName(void) {

    int maxLen = 0;
    for (int i=0; i<matrix.size(); i++)
        {
        std::string n = matrix[i].getName();
        if (n.length() > maxLen)
            maxLen = (int)n.length();
        }
    return maxLen;
}

void Alignment::print(std::string h) {

    std::cout << h << std::endl;
    print();
}

void Alignment::print(void) {

    int maxLen = lengthOfLongestName();        
    for (int i=0; i<matrix.size(); i++)
        {
        std::string n = matrix[i].getName();
        std::cout << n << " ";
        for (int j=0; j<maxLen-n.length(); j++)
            std::cout << " ";
        for (int j=0; j<matrix[i].size(); j++)
            {
            int charCode = matrix[i][j];
            if (charCode == -1)
                std::cout << " - ";
            else
                std::cout << std::setw(2) << charCode << " ";
            }
        std::cout << std::endl;
        }
}
