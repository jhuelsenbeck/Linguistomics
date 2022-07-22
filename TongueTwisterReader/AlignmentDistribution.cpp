#include <iomanip>
#include <iostream>
#include <vector>
#include "Alignment.hpp"
#include "AlignmentDistribution.hpp"
#include "Partition.hpp"
#include "RandomVariable.hpp"



bool cmp(std::pair<Alignment*, int>& a, std::pair<Alignment*, int>& b) {

    return a.second > b.second;
}



AlignmentDistribution::AlignmentDistribution(RandomVariable* r, Partition* p) {

    rv = r;
    partition = p;
}

AlignmentDistribution::~AlignmentDistribution(void) {

    if (partition != NULL)
        delete partition;
}

void AlignmentDistribution::addAlignment(Alignment* aln) {

    bool foundAlignmentInMap = false;
    for (std::map<Alignment*,int>::iterator it = samples.begin(); it != samples.end(); it++)
        {
        if (*it->first == *aln)
            {
            it->second++;
            foundAlignmentInMap = true;
            delete aln;
            break;
            }
        }
        
    if (foundAlignmentInMap == false)
        {
        samples.insert( std::make_pair(aln,1) );
        }
}

int AlignmentDistribution::numSamples(void) {

    int n = 0;
    for (auto& it : samples)
        n += it.second;
    return n;
}

int AlignmentDistribution::ciSize(void) {

    std::vector<std::pair<Alignment*, int> > v;
    for (auto& it : samples)
        v.push_back(it);
  
    sort(v.begin(), v.end(), cmp);

    int n = numSamples();
    double cumulativeProb = 0.0;
    int i = 0;
    for (auto& it : v)
        {
        double prob = (double)it.second / n;
        cumulativeProb += prob;
        if (cumulativeProb < 0.95)
            {
            i++;
            }
        else
            {
            double lowerVal = cumulativeProb - prob;
            double u = rv->uniformRv();
            if (u < (0.95-lowerVal)/(cumulativeProb-lowerVal))
                {
                i++;
                }
            break;
           
            }
        }
    return i;
}

Alignment* AlignmentDistribution::getMapAlignment(void) {

    std::vector<std::pair<Alignment*, int> > v;
    for (auto& it : samples)
        v.push_back(it);
    sort(v.begin(), v.end(), cmp);
    return v[0].first;
}

void AlignmentDistribution::print(void) {

    std::vector<std::pair<Alignment*, int> > v;
    for (auto& it : samples)
        v.push_back(it);
  
    sort(v.begin(), v.end(), cmp);
  
    int n = numSamples();
    double cumulativeProb = 0.0;
    int i = 0;
    std::cout << std::endl << "Cognate: " << name << std::endl;
    for (auto& it : v)
        {
        double prob = (double)it.second / n;
        cumulativeProb += prob;
        if (cumulativeProb < 0.95)
            {
            //it.first->print( std::to_string(++i) + " Probability: " + std::to_string(prob) + " (" + std::to_string(cumulativeProb) + ")" );
            std::cout << ++i << " -- Probability: " << prob << " (" << cumulativeProb << ")" << std::endl;
            print(it.first);
            }
        else
            {
            double lowerVal = cumulativeProb - prob;
            double u = rv->uniformRv();
            if (u < (0.95-lowerVal)/(cumulativeProb-lowerVal))
                {
                //it.first->print( std::to_string(++i) + " Probability: " + std::to_string(prob) + " (" + std::to_string(cumulativeProb) + ")" );
                std::cout << ++i << " -- Probability: " << prob << " (" << cumulativeProb << ")" << std::endl;
                print(it.first);
                }
            break;
            }
        }
}

void AlignmentDistribution::print(Alignment* aln) {

    int len = aln->lengthOfLongestName();
    for (int i=0; i<aln->getNumTaxa(); i++)
        {
        std::string tName = (*aln)[i].getName();
        std::cout << tName << " ";
        for (int j=0; j<len-tName.length(); j++)
            std::cout << " ";
        for (int j=0; j<aln->getNumChar(); j++)
            {
            int charCode = aln->getCharCode(i, j);
            int subsetIndex = charCode;
            if (subsetIndex != -1 && partition != NULL)
                subsetIndex = partition->indexOfSubsetWithValue(charCode);
            if (charCode == -1)
                std::cout << " - ";
            else
                std::cout << std::setw(2) << subsetIndex << " ";
            }
        std::cout << std::endl;
        }
}

nlohmann::json AlignmentDistribution::toJson(double credibleSetSize, std::ostream& findex) {
    
    // sort the alignments from highest to lowest posterior probability
    std::vector<std::pair<Alignment*, int> > v(1024);
    for (auto& it : samples) 
        {
        if (it.first->valid())
          v.push_back(it);
        }
    sort(v.begin(), v.end(), cmp);


    findex << "new(\"";
    for (auto& c : name)
        {
        if (c == '-')
          findex << "\", \"";
        else
          findex << c;
        }


    findex << "\", [\n";

    auto j = nlohmann::json::object();
    j["cognate"] = name;
    auto alnVec = nlohmann::json::array();                 
    int n = numSamples();
    double cumulativeProb = 0.0;
    int i = 0;
    int count = 0;
    for (auto& it : v)
        {
        ++i;
        double prob = (double)it.second / n;
        cumulativeProb += prob;
        
        findex << "new(" << prob << ",[";


        auto jaln = nlohmann::json::object();
        jaln["index"] = i;
        jaln["prob"] = prob;
        jaln["cumprob"] = cumulativeProb;

        bool done = false;
        bool add  = true;
                
        if (cumulativeProb >= credibleSetSize) 
            {
            double lowerVal = cumulativeProb - prob;
            double u = rv->uniformRv();
            add  = u < (credibleSetSize-lowerVal)/(cumulativeProb-lowerVal);
            done = true;
            }

        if (add) 
            {
            jaln["aln"] = it.first->toJson(findex);
            alnVec.push_back(jaln);
            }

        ++count;
        findex << "]),\n";
        
        if (done)
            break;
        }


    j["aln_set"] = alnVec;

    findex << "]),\n";
    return j;
}

