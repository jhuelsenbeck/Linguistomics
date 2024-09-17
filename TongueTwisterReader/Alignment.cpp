#include <iomanip>
#include <iostream>
#include "Alignment.hpp"
#include "Sequence.hpp"


Alignment::Alignment(nlohmann::json js, std::vector <std::string>& taxa) : sorted(taxa.size()) {

    for (int i=0; i<js.size(); i++)
        {
        const std::string& tName = js[i]["Taxon"];
        const std::vector<int>& tInfo = js[i]["Segments"];

        int j = 0;
        for (auto& t: taxa)
            {
            if (t == tName)
                break;
            ++j;
            }


        matrix.push_back( Sequence(j, tName, tInfo) );
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

void Alignment::from_json(const nlohmann::json& /*j*/) {

    //j.at("name").get_to(p.name);
    //j.at("address").get_to(p.address);
    //j.at("age").get_to(p.age);
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

std::vector<double> Alignment::getGapSpectrum(void) {

    std::vector<double> spectrum;
    if (numChar == 1)
        {
        spectrum.resize(1);
        int n = 0;
        for (int i=0; i<numTaxa; i++)
            {
            if (matrix[i][0] == -1)
                n++;
            }
        spectrum[0] = n;
        }
    else if (numChar == 2)
        {
        spectrum.resize(2);
        for (int j=0; j<2; j++)
            {
            int n = 0;
            for (int i=0; i<numTaxa; i++)
                {
                if (matrix[i][j] == -1)
                    n++;
                }
            spectrum[j] = n;
            }
        }
    else
        {
        spectrum.resize(3);
        for (int i=0; i<3; i++)
            spectrum[i] = 0.0;
        for (int j=0; j<numChar; j++)
            {
            int n = 0;
            for (int i=0; i<numTaxa; i++)
                {
                if (matrix[i][j] == -1)
                    n++;
                }
            if (j == 0)
                spectrum[0] = n;
            else if (j == numChar-1)
                spectrum[2] = n;
            else
                spectrum[1] += n;
            }
        spectrum[1] /= (numChar-2);
        }
    
    return spectrum;
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

nlohmann::json Alignment::toJson(std::ostream& nytril) {

    nlohmann::json j = nlohmann::json::array();

    auto size = sorted.size();

    for (int i = 0; i < size; ++i)
        sorted[i] = NULL;

    for (int i=0; i<numTaxa; i++)
        {
        auto& m = matrix[i];
        nlohmann::json sj;
        sj["language"] = m.getName();
        sj["sequence"] = m.getSequence();
        j.push_back(sj);

        sorted[m.getLanguage()] = &m;
        }


    bool com = false;
    for (int i = 0; i < size; i++) {
        auto m = sorted[i];
        if (com)
            nytril << ",";
        if (m) {
            nytril << "[";
            bool comma = false;
            for (auto si : m->getSequence())
            {
                if (comma)
                    nytril << ",";
                nytril << si;
                comma = true;
            }
            nytril << "]";
        }
        else
            nytril << "null";
        com = true;
    }

    return j;
}
