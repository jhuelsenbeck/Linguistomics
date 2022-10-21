#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include "Alignment.hpp"
#include "AlignmentDistribution.hpp"
#include "Container.hpp"
#include "json.hpp"
#include "McmcSummary.hpp"
#include "Msg.hpp"
#include "Partition.hpp"
#include "RandomVariable.hpp"
#include "Subset.hpp"
#include "Tree.hpp"



McmcSummary::McmcSummary(RandomVariable* r) : conTree(NULL), rv(r) {

    statePartitions = NULL;
}

McmcSummary::~McmcSummary(void) {

    for (int i=0; i<stats.size(); i++)
        delete stats[i];
    if (statePartitions != NULL)
        delete statePartitions;
}

McmcSummary& McmcSummary::operator+=(const McmcSummary& rhs) {

    // combine the ParameterStatistics
    for (int i=0; i<stats.size(); i++)
        {
        std::string thisName = stats[i]->getName();
        ParameterStatistics* rhsStats = rhs.getParameterNamed(thisName);
        if (rhsStats != NULL)
            *(stats[i]) += *rhsStats;
        }
        
    // combine the alignments
    for (int i=0; i<alignments.size(); i++)
        {
        std::string thisName = alignments[i]->getName();
        AlignmentDistribution* rhsAln = rhs.getAlignmentNamed(thisName);
        if (rhsAln != NULL)
            *(alignments[i]) += *rhsAln;
        }
        
    // combine the tree partition information
    for (std::map<RbBitSet,ParameterStatistics*>::const_iterator it = rhs.partitions.begin(); it != rhs.partitions.end(); it++)
        {
        std::map<RbBitSet,ParameterStatistics*>::iterator partitionsIt = this->partitions.find(it->first);
        if (partitionsIt == this->partitions.end())
            {
            ParameterStatistics* newStats = new ParameterStatistics(*(it->second));
            partitions.insert( std::make_pair(it->first,newStats) );
            }
        else
            {
            *(partitionsIt->second) += *(it->second);
            }
        }

    return *this;
}

void McmcSummary::addPartion(std::map<RbBitSet,double>& parts) {

    for (std::map<RbBitSet,double>::iterator it = parts.begin(); it != parts.end(); it++)
        {
        std::map<RbBitSet,ParameterStatistics*>::iterator partitionsIt = partitions.find(it->first);
        if (partitionsIt == partitions.end())
            {
            ParameterStatistics* newStats = new ParameterStatistics;
            newStats->addValue(it->second);
            newStats->setName(it->first.bitString());
            partitions.insert( std::make_pair(it->first,newStats) );
            }
        else
            {
            partitionsIt->second->addValue(it->second);
            }
        }
}

std::vector<std::string> McmcSummary::breakString(std::string str) {

    std::vector<std::string> broken;
    std::string word = "";
    for (int i=0; i<str.length(); i++)
        {
        char c = str[i];
        if (c == ',' || c == ';')
            {
            if (word != "")
                broken.push_back(word);
            broken.push_back( std::string(1,c) );
            word = "";
            }
        else
            {
            word +=  c;
            }
        }
    if (word != "")
        broken.push_back(word);
    
    return broken;
}

void McmcSummary::calculateRates(DoubleMatrix& m) {

    int numStates = 0;
    
    // 1. get state frequencies (and the number of states). We look
    // up the apropriate ParameterStatistics objects for this
    std::set<int> observedStates;
    for (int i=0; i<stats.size(); i++)
        {
        ParameterStatistics* s = stats[i];
        if (s->getName()[0] == 'F')
            {
            int st = getFreqElement(s->getName());
            observedStates.insert(st);
            }
        }
    for (int i=0; i<observedStates.size(); i++) // make certain all of the states from 0 to (numStates-1) are found
        {
        std::set<int>::iterator it = observedStates.find(i);
        if (it == observedStates.end())
            Msg::error("Could not find state " + std::to_string(i) + " in frequencies summary");
        }
    numStates = (int)observedStates.size();
    std::vector<double> f(numStates);
    for (int i=0; i<stats.size(); i++)
        {
        ParameterStatistics* s = stats[i];
        if (s->getName()[0] == 'F')
            {
            int st = getFreqElement(s->getName());
            f[st] = s->getMean();
            }
        }
    if (statePartitions != NULL) // check consistencey with a state partition object, if we have one
        {
        int tempNumStates = statePartitions->getNumElements();
        if (tempNumStates != numStates)
            Msg::error("Mismatch in the number of states from the frequencies and the state partitions");
        }

    // allocate a vector to hold the exchangeability parameters
    DoubleMatrix r(numStates,numStates);
    for (int i=0; i<numStates; i++)
        for (int j=0; j<numStates; j++)
            r(i,j) = -1.0;
    for (int i=0; i<stats.size(); i++)
        {
        ParameterStatistics* s = stats[i];
        if (s->getName()[0] == 'R')
            {
            int r1, r2;
            getRateElements(s->getName(), r1, r2);
            if (statePartitions != NULL)
                {
                Subset* s1 = statePartitions->findSubsetIndexed(r1);
                Subset* s2 = statePartitions->findSubsetIndexed(r2);
                for (int x1 : s1->getValues())
                    {
                    for (int x2 : s2->getValues())
                        {
                        double ave = s->getMean();
                        r(x1,x2) = ave;
                        r(x2,x1) = ave;
                        }
                    }
                
                }
            else
                {
                r(r1,r2) = s->getMean();
                }
            }
        }
    for (int i=0; i<numStates; i++)
        {
        for (int j=0; j<numStates; j++)
            {
            if (r(i,j) < 0.0 && i != j)
                Msg::error("Could not find all rate parameters");
            }
        }


    // fill in the rates (Q)
    //DoubleMatrix m(numStates,numStates);
    double averageRate = 0.0;
    for (int i=0; i<numStates; i++)
        {
        double sum = 0.0;
        for (int j=0; j<numStates; j++)
            {
            if (i != j)
                {
                m(i,j) = r(i,j) * f[j];
                sum += m(i,j);
                }
            }
        m(i,i) = -sum;
        averageRate += f[i] * sum;
        }
        
    // rescale so average is one (Q)
    double scaleFactor = 1.0 / averageRate;
    for (int i=0; i<numStates; i++)
        for (int j=0; j<numStates; j++)
            m(i,j) *= scaleFactor;
                        
#   if 1
    double sum = 0.0;
    std::cout << "Frequencies:" << std::endl;
    for (int i=0; i<f.size(); i++)
        {
        std::cout << "F[" << i << "] = " << f[i] << std::endl;
        sum += f[i];
        }
    std::cout << "Frequencies sum = " << sum << std::endl;
    std::cout << "Rates:" << std::endl;
    sum = 0.0;
    for (int i=0; i<numStates; i++)
        {
        for (int j=0; j<numStates; j++)
            {
            std::cout << r(i,j) << " ";
            sum += r(i,j);
            }
        std::cout << std::endl;
        }
    std::cout << "Rates sum = " << sum << std::endl;
    std::cout << "Average Rates:" << std::endl;
    sum = 0.0;
    for (int i=0; i<numStates; i++)
        {
        for (int j=0; j<numStates; j++)
            {
            if (i != j)
                {
                std::cout << m(i,j) << " ";
                sum += m(i,j);
                }
            }
        std::cout << std::endl;
        }
    std::cout << "Average rates  um = " << sum << std::endl;
#   endif

}

void McmcSummary::calculateAverageRates(DoubleMatrix& m) {

    int numStates = 0;
    
    // 1. get state frequencies (and the number of states). We look
    // up the apropriate ParameterStatistics objects for this
    std::set<int> observedStates;
    for (int i=0; i<stats.size(); i++)
        {
        ParameterStatistics* s = stats[i];
        if (s->getName()[0] == 'F')
            {
            int st = getFreqElement(s->getName());
            observedStates.insert(st);
            }
        }
    for (int i=0; i<observedStates.size(); i++) // make certain all of the states from 0 to (numStates-1) are found
        {
        std::set<int>::iterator it = observedStates.find(i);
        if (it == observedStates.end())
            Msg::error("Could not find state " + std::to_string(i) + " in frequencies summary");
        }
    numStates = (int)observedStates.size();
    std::vector<double> f(numStates);
    for (int i=0; i<stats.size(); i++)
        {
        ParameterStatistics* s = stats[i];
        if (s->getName()[0] == 'F')
            {
            int st = getFreqElement(s->getName());
            f[st] = s->getMean();
            }
        }
    if (statePartitions != NULL) // check consistencey with a state partition object, if we have one
        {
        int tempNumStates = statePartitions->getNumElements();
        if (tempNumStates != numStates)
            Msg::error("Mismatch in the number of states from the frequencies and the state partitions");
        }

    // allocate a vector to hold the exchangeability parameters
    DoubleMatrix r(numStates,numStates);
    for (int i=0; i<numStates; i++)
        for (int j=0; j<numStates; j++)
            r(i,j) = -1.0;
    for (int i=0; i<stats.size(); i++)
        {
        ParameterStatistics* s = stats[i];
        if (s->getName()[0] == 'R')
            {
            int r1, r2;
            getRateElements(s->getName(), r1, r2);
            if (statePartitions != NULL)
                {
                Subset* s1 = statePartitions->findSubsetIndexed(r1);
                Subset* s2 = statePartitions->findSubsetIndexed(r2);
                for (int x1 : s1->getValues())
                    {
                    for (int x2 : s2->getValues())
                        {
                        double ave = s->getMean();
                        r(x1,x2) = ave;
                        r(x2,x1) = ave;
                        }
                    }
                
                }
            else
                {
                r(r1,r2) = s->getMean();
                }
            }
        }
    for (int i=0; i<numStates; i++)
        {
        for (int j=0; j<numStates; j++)
            {
            if (r(i,j) < 0.0 && i != j)
                Msg::error("Could not find all rate parameters");
            }
        }


    // fill in the rates (Q)
    //DoubleMatrix m(numStates,numStates);
    double averageRate = 0.0;
    for (int i=0; i<numStates; i++)
        {
        double sum = 0.0;
        for (int j=0; j<numStates; j++)
            {
            if (i != j)
                {
                m(i,j) = r(i,j) * f[j];
                sum += m(i,j);
                }
            }
        m(i,i) = -sum;
        averageRate += f[i] * sum;
        }
        
    // rescale so average is one (Q)
    double scaleFactor = 1.0 / averageRate;
    for (int i=0; i<numStates; i++)
        for (int j=0; j<numStates; j++)
            m(i,j) *= scaleFactor;
            
    // now get average rate from i to j (R)
    for (int i=0; i<numStates; i++)
        for (int j=0; j<numStates; j++)
            m(i,j) *= f[i];
            
#   if 1
    double sum = 0.0;
    std::cout << "Frequencies:" << std::endl;
    for (int i=0; i<f.size(); i++)
        {
        std::cout << "F[" << i << "] = " << f[i] << std::endl;
        sum += f[i];
        }
    std::cout << "Frequencies sum = " << sum << std::endl;
    std::cout << "Rates:" << std::endl;
    sum = 0.0;
    for (int i=0; i<numStates; i++)
        {
        for (int j=0; j<numStates; j++)
            {
            std::cout << r(i,j) << " ";
            sum += r(i,j);
            }
        std::cout << std::endl;
        }
    std::cout << "Rates sum = " << sum << std::endl;
    std::cout << "Average Rates:" << std::endl;
    sum = 0.0;
    for (int i=0; i<numStates; i++)
        {
        for (int j=0; j<numStates; j++)
            {
            if (i != j)
                {
                std::cout << m(i,j) << " ";
                sum += m(i,j);
                }
            }
        std::cout << std::endl;
        }
    std::cout << "Average rates  um = " << sum << std::endl;
#   endif
}

int McmcSummary::getFreqElement(std::string str) {

    int n = 0;
    std::string tempStr = "";
    bool readingBetweenSquareBrackets = false;
    for (int i=0; i<str.length(); i++)
        {
        char c = str[i];
        if (c == '[')
            readingBetweenSquareBrackets = true;
        else if (c == ']')
            readingBetweenSquareBrackets = false;
        else
            {
            if (readingBetweenSquareBrackets == true)
                tempStr += std::string(1, c);
            }
        }
    if (tempStr != "")
        n = std::stoi(tempStr);
    return n;
}

void McmcSummary::getRateElements(std::string str, int& r1, int& r2) {

    r1 = 0;
    r2 = 0;

    std::string tempStr = "";
    bool readingFirstNum = false, readingSecondNum = false;
    for (int i=0; i<str.length(); i++)
        {
        char c = str[i];
        if (c == '[')
            {
            readingFirstNum = true;
            readingSecondNum = false;
            }
        else if (c == ']')
            {
            readingFirstNum = false;
            readingSecondNum = false;
            if (tempStr != "")
                r2 = std::stoi(tempStr);
            tempStr = "";
            }
        else if (c == '-')
            {
            readingFirstNum = false;
            readingSecondNum = true;
            if (tempStr != "")
                r1 = std::stoi(tempStr);
            tempStr = "";
            }
        else
            {
            if (readingFirstNum == true || readingSecondNum == true)
                tempStr += std::string(1, c);
            }
        }
}

std::vector<CredibleInterval> McmcSummary::getCredibleIntervals(void) {

    std::vector<CredibleInterval> cis;
    for (int i=0; i<stats.size(); i++)
        cis.push_back( stats[i]->getCredibleInterval() );
    return cis;
}

std::vector<double> McmcSummary::getMeans(void) {

    std::vector<double> means;
    for (int i=0; i<stats.size(); i++)
        means.push_back( stats[i]->getMean() );
    return means;
}

std::string McmcSummary::getCognateName(std::string str) {

    std::vector<std::string> segs;
    std::string word = "";
    for (int i=0; i<str.length(); i++)
        {
        char c = str[i];
        if (c == '/' || c == '\\' || c == '.')
            {
            if (word != "")
                {
                segs.push_back(word);
                word = "";
                }
            }
        else
            {
            word += c;
            }
        }
    if (word != "")
        segs.push_back(word);
        
//    for (int i=0; i<segs.size(); i++)
//        {
//        std::cout << i << " -- " << segs[i] << std::endl;
//        }
    return segs[segs.size()-2];
}

AlignmentDistribution* McmcSummary::getAlignmentNamed(std::string str) const {

    for (AlignmentDistribution* a : alignments)
        {
        if (a->getName() == str)
            return a;
        }
    return NULL;
}

ParameterStatistics* McmcSummary::getParameterNamed(std::string str) const {

    for (ParameterStatistics* s : stats)
        {
        if (s->getName() == str)
            return s;
        }
    return NULL;
}

bool McmcSummary::hasSemicolon(std::string str) {

    std::size_t found = str.find(';');
    if (found != std::string::npos)
        return true;
    return false;
}

int McmcSummary::inferNumberOfRates(void) {

    int numRateClasses = 0;
    for (int i=0; i<stats.size(); i++)
        {
        ParameterStatistics* s = stats[i];
        if (s->getName()[0] == 'R')
            numRateClasses++;
        }
    if (statePartitions != NULL) // check consistencey with a state partition object, if we have one
        {
        int numSubsets = statePartitions->numSubsets();
        int tempNumRates = numSubsets * (numSubsets-1) / 2;
        std::set<Subset*>& subsets = statePartitions->getSubsets();
        for (Subset* s : subsets)
            {
            if (s->getNumElements() > 1)
                tempNumRates++;
            }
        if (tempNumRates != numRateClasses)
            Msg::error("Mismatch in the number of rate classes");
        }
    return numRateClasses;
}

int McmcSummary::inferNumberOfStates(void) {

    // 1. get state frequencies (and the number of states). We look
    // up the apropriate ParameterStatistics objects for this
    std::set<int> observedStates;
    for (int i=0; i<stats.size(); i++)
        {
        ParameterStatistics* s = stats[i];
        if (s->getName()[0] == 'F')
            {
            int st = getFreqElement(s->getName());
            observedStates.insert(st);
            }
        }
    for (int i=0; i<observedStates.size(); i++) // make certain all of the states from 0 to (numStates-1) are found
        {
        std::set<int>::iterator it = observedStates.find(i);
        if (it == observedStates.end())
            Msg::error("Could not find state " + std::to_string(i) + " in frequencies summary");
        }
    int numStates = (int)observedStates.size();
    if (statePartitions != NULL) // check consistencey with a state partition object, if we have one
        {
        int tempNumStates = statePartitions->getNumElements();
        if (tempNumStates != numStates)
            Msg::error("Mismatch in the number of states from the frequencies and the state partitions");
        }
    return numStates;
}

std::map<int,std::string> McmcSummary::interpretTranslateString(std::vector<std::string> translateTokens) {

    std::map<int,std::string> translateMap;
    bool readingKey = true;
    std::string key = "";
    std::string val = "";
    for (int i=0; i<translateTokens.size(); i++)
        {
        std::string token = translateTokens[i];
        if (token == "," || token == ";")
            {
            readingKey = true;
            int intKey = atoi(key.c_str());
            translateMap.insert( std::make_pair(intKey,val) );
            }
        else
            {
            if (readingKey == true)
                {
                key = token;
                readingKey = false;
                }
            else
                val = token;
            }
        //std::cout << i << " " << token << std::endl;
        }
//    for (std::map<int,std::string>::iterator it = translateMap.begin(); it != translateMap.end(); it++)
//        std::cout << it->first << " -> " << it->second << std::endl;
        
    return translateMap;
}

std::string McmcSummary::interpretTreeString(std::string str) {

    std::string ns = "";
    bool reading = false;
    for (int i=0; i<str.length(); i++)
        {
        char c = str[i];
        if (c == '(')
            reading = true;
            
        if (reading == true)
            {
            ns += c;
            }
            
        if (c == ';')
            reading = false;
        }
    
    return ns;
}

void McmcSummary::print(void) {

    // print summary of tsv file
    std::vector<CredibleInterval> cis = getCredibleIntervals();
    std::vector<double> means = getMeans();
    std::cout << std::fixed << std::setprecision(4);
    for (int i=0; i<stats.size(); i++)
        {
        std::cout << stats[i]->getName() << ": ";
        std::cout << cis[i].median << " (";
        std::cout << cis[i].lower << ", " << cis[i].upper << ")";
        std::cout << std::endl;
        }
        
    if (statePartitions != NULL)
        printPartitionFreqs();
        
    // print summary of aln files
    for (int i=0; i<alignments.size(); i++)
        {
        alignments[i]->print();
        }
    
    // print summary of tre file
    int n = 0;
    for (auto it = partitions.begin(); it != partitions.end(); it++)
        {
        if (it->second->size() > n)
            n = it->second->size();
        }
    std::cout << std::fixed << std::setprecision(6);
    for (auto it = partitions.begin(); it != partitions.end(); it++)
        {
        auto ci = it->second->getCredibleInterval();
        std::cout << it->second->getName() << " ";
        std::cout << (double)it->second->size() / n << " ";
        std::cout << it->second->getMean() << " (" << ci.lower << ", " << ci.upper << ")";
        std::cout << std::endl;
        }
    conTree->print();
    std::cout << conTree->getNewick(4) << std::endl;
    
    std::vector<int> dist(1000,0);
    int maxSize = 0;
    std::cout << "raw = c(";
    for (int i=0; i<alignments.size(); i++)
        {
        int ciSize = alignments[i]->ciSize();
        if (ciSize > maxSize)
            maxSize = ciSize;
        dist[ciSize]++;

        std::cout << ciSize;
        if (i + 1 != alignments.size())
            std::cout << ",";
        }
    std::cout << ")" << std::endl;

    for (int i=1; i<maxSize; i++)
        std::cout << i << " -- " << dist[i] << std::endl;
    std::cout << "v = c(";
    for (int i=1; i<maxSize; i++)
        {
        std::cout << i;
        if (i + 1 != maxSize)
            std::cout << ",";
        }
    std::cout << ")" << std::endl;
    
    std::cout << "x = c(";
    for (int i=1; i<maxSize; i++)
        {
        std::cout << dist[i];
        if (i + 1 != maxSize)
            std::cout << ",";
        }
    std::cout << ")" << std::endl;
    
    std::cout << "x <- c(";
    for (int i=0; i<alignments.size(); i++)
        {
        Alignment* map = alignments[i]->getMapAlignment();
        std::map<double,double> gaps = map->gapInfo();
        for (std::map<double,double>::iterator it = gaps.begin(); it != gaps.end(); it++)
            std::cout << it->first << ",";
        }
    std::cout << std::endl;
    std::cout << "y <- c(";
    for (int i=0; i<alignments.size(); i++)
        {
        Alignment* map = alignments[i]->getMapAlignment();
        std::map<double,double> gaps = map->gapInfo();
        for (std::map<double,double>::iterator it = gaps.begin(); it != gaps.end(); it++)
            std::cout << it->second << ",";
        }
    std::cout << std::endl;

}

void McmcSummary::printPartitionSet(void) {

    if (statePartitions != NULL)
        statePartitions->print();
}

void McmcSummary::printPartitionFreqs(void) {

    if (statePartitions != NULL)
        {
        std::vector<double> subsetFreqs(statePartitions->numSubsets(), 0.0);
        for (int i=0; i<stats.size(); i++)
            {
            std::string statName = stats[i]->getName();
             if (statName[0] == 'F')
                {
                statName.erase(remove(statName.begin(), statName.end(), 'F'), statName.end());
                statName.erase(remove(statName.begin(), statName.end(), '['), statName.end());
                statName.erase(remove(statName.begin(), statName.end(), ']'), statName.end());
                int idx = atoi(statName.c_str());
                int pIdx = statePartitions->indexOfSubsetWithValue(idx);
                subsetFreqs[pIdx-1] += stats[i]->getMean();
                //std::cout << statName << " " << idx << " " << pIdx << std::endl;
                }
            }
            
        double maxFreq = 0.0;
        for (int i=0; i<subsetFreqs.size(); i++)
            {
            if (subsetFreqs[i] > maxFreq)
                maxFreq = subsetFreqs[i];
            }
        double maxR = sqrt(maxFreq / 3.14);

        for (int i=0; i<subsetFreqs.size(); i++)
            {
            double r = (sqrt(subsetFreqs[i] / 3.14) / maxR) * 100.0;
            std::cout << "Subset " << i+1 << " freqeuncy = " << subsetFreqs[i] << " " << r << std::endl;
            }
            
            
            
            
        
        double maxVal = 0.0;
        for (int i=0; i<stats.size(); i++)
            {
            std::string statName = stats[i]->getName();
             if (statName[0] == 'R')
                {
                if (stats[i]->getMean() > maxVal)
                    maxVal = stats[i]->getMean();
                }
            }
        std::cout << "maxVal = " << maxVal << std::endl;
        maxR = sqrt(maxVal / 3.14);
        for (int i=0; i<stats.size(); i++)
            {
            std::string statName = stats[i]->getName();
             if (statName[0] == 'R')
                {
                double r = (sqrt(stats[i]->getMean() / 3.14) / maxR) * 36.0;
                std::cout << stats[i]->getName() << " " << stats[i]->getMean() << " " << r << std::endl;
                }
            }
           
            

        }
}

void McmcSummary::output(std::string pathName, std::ofstream& findex) {

    auto stree = conTree->getNewick(4);
    auto tree  = new std::ofstream(pathName + "/consensus.tre", std::ios::out);
    *tree << stree << ";";
    tree->close();
    delete tree;

    double cutoff = 0.95;

    findex << "// This file is auto-generated by TongueTwisterReader. Do not edit\n\n";
    findex << "AlignmentCutoff = " << cutoff << ";\n\n";

    auto j = nlohmann::json::object();
    
    // output state partition information
    if (statePartitions != NULL)
        j["state_part"] = statePartitions->toJson();

    // output the consensus tree
    j["consensus"]["tree"] = stree;
    j["consensus"]["n_taxa"] = conTree->numTaxa;
    
    // output information on mean and credible interval for all real-valued parameters
    auto jStats = nlohmann::json::array();
    for (int i=0; i<stats.size(); i++)
        {
        nlohmann::json cogStats = stats[i]->toJson();
        jStats.push_back(cogStats);
        }
    j["stats"] = jStats;


    int  matrix = 6;
    auto pcount = 10;

    findex << "StatClass[][] TransitionStats = [\n";
    for (int mi = 0; mi < pcount; ++mi) {
        if (mi)
            findex << ",\n";
        findex << "  [";
        for (int mj = 0; mj < pcount; ++mj) {
            if (mj)
                findex << ",";
            if (mj < mi)
                findex << "null";
            else
                stats[matrix++]->toFile(findex);
        }
        findex << "]";
    }
    findex << "\n];\n\n";
    
    int numRateClasses = inferNumberOfRates();
    std::cout << "numRateClasses = " << numRateClasses << std::endl;
    
    // output average rates of change
    // No credible intervals on this information.
    int numStates = inferNumberOfStates();
    DoubleMatrix m(numStates,numStates);
    calculateAverageRates(m);
    findex << "AverageRates = [\n";
    for (int mi = 0; mi < numStates; ++mi)
        {
        if (mi)
            findex << ",\n";
        findex << "  [";
        for (int mj = 0; mj < numStates; ++mj)
            {
            double r = m(mi,mj);
            if (mj)
                findex << ",";
            findex << r;
            }
        findex << "]";
        }
    findex << "\n];\n\n";
    
    DoubleMatrix q(numStates,numStates);
    calculateRates(q);
    findex << "QRates = [\n";
    for (int mi = 0; mi < numStates; ++mi)
        {
        if (mi)
            findex << ",\n";
        findex << "  [";
        for (int mj = 0; mj < numStates; ++mj)
            {
            double r = q(mi,mj);
            if (mj)
                findex << ",";
            findex << r;
            }
        findex << "]";
        }
    findex << "\n];\n\n";


    findex << "AlignIndexClass[] Alignments = [\n";

        
    // output the credible set for all alignments
    auto sampledAlignments = nlohmann::json::array();
    for (int i=0; i<alignments.size(); i++)
        {
        auto& a = alignments[i];
        if (a->size() > 0)
            {
            auto cogAlns = a->toJson(cutoff, findex);
            sampledAlignments.push_back(cogAlns);
            }
        }
    j["cog_alns"] = sampledAlignments;

    findex << "];\n\n";

    // output the pooled partition subset frequencies (if a subset is available)
    if (statePartitions != NULL)
        {
        std::map<int,double> subsetFreqs;
        
        for (int i=0; i<stats.size(); i++)
            {
            std::string n = stats[i]->getName();
            if (n[0] == 'F' && n[1] == '[')
                {
                int x = parseNumberFromFreqHeader(n);
                double val = stats[i]->getMean();
                
                //std::cout << n << " " << x << std::endl;
                Subset* s = statePartitions->findSubsetWithValue(x);
                int subsetIndex = s->getIndex();
                if (s == NULL)
                    Msg::error("Could not find subset with value " + std::to_string(x));
                
                std::map<int,double>::iterator it = subsetFreqs.find(subsetIndex);
                if (it == subsetFreqs.end())
                    {
                    subsetFreqs.insert( std::make_pair(subsetIndex, val) );
                    }
                else
                    {
                    it->second += val;
                    }
                
                }
            }
            double sum = 0.0;
            for (std::map<int,double>::iterator it = subsetFreqs.begin(); it != subsetFreqs.end(); it++)
                {
                std::cout << it->first << " " << it->second << std::endl;
                sum += it->second;
                }
            std::cout << "sum = " << sum << std::endl;
        j["part_freqs"] = statePartitions->toJson(subsetFreqs, findex);
        }
        
    // output json representation to a file
    auto file = new std::ofstream(pathName + "/summary.json", std::ios::out);
    *file << j;
    std::cout << "JSON file written to: " << pathName << "summary.json" << std::endl;
    file->close();
    delete file;
}
    
int McmcSummary::parseNumberFromFreqHeader(std::string str) {

    bool parse = false;
    std::string numStr = "";
    for (int i=0; i<str.length(); i++)
        {
        char c = str[i];
        if (c == '[')
            parse = true;
        else if (c == ']')
            parse = false;
        if (c != '[' && c != ']' && parse == true)
            numStr += std::string(1, c);
        }
    int num = atoi(numStr.c_str());
    return num;
}

void McmcSummary::readAlnFile(std::string fn, int bi) {

    std::ifstream ifs(fn);
    nlohmann::json j;
    try
        {
        j = nlohmann::json::parse(ifs);
        }
    catch (nlohmann::json::parse_error& ex)
        {
        Msg::warning("Error parsing JSON file" + fn + " at byte " + std::to_string(ex.byte));
        return;
        }
        
    std::string cognateName = getCognateName(fn);

    auto it = j.find("PartitionSets");
    hasPartitions = true;
    Partition* part = NULL;
    if (it == j.end())
        {
        hasPartitions = false;
        //Msg::error("Could not find partition set in the JSON file");
        }
    if (hasPartitions == true)
        {
        nlohmann::json jsPart = j["PartitionSets"];
        part = new Partition(jsPart);
        }

    AlignmentDistribution* dist = new AlignmentDistribution(rv, part);
    alignments.push_back(dist);
    dist->setName(cognateName);
    

    it = j.find("Samples");
    if (it == j.end())
        Msg::error("Could not find sampled alignments in the JSON file");
    nlohmann::json js = j["Samples"];
    for (int i=0; i<js.size(); i++)
        {
        if (i > bi)
            {
            Alignment* aln = new Alignment( js[i]["Data"], taxa );
            dist->addAlignment(aln);
            }
        }
        
    //dist->print();
}

void McmcSummary::readConfigFile(std::string fn) {

    std::ifstream ifs(fn);
    nlohmann::json j;
    try
        {
        j = nlohmann::json::parse(ifs);
        }
    catch (nlohmann::json::parse_error& ex)
        {
        Msg::warning("Error parsing JSON file" + fn + " at byte " + std::to_string(ex.byte));
        return;
        }

    auto jtaxa = j["Taxa"];
    for (int i = 0; i < jtaxa.size(); i++)
        taxa.push_back(jtaxa[i].get<std::string>());



    auto it = j.find("PartitionSets");
    if (it == j.end())
        {
        Msg::warning("Could not find partition set in the JSON configuration file");
        statePartitions = NULL;
        }
    else
        {
        nlohmann::json jsPart = j["PartitionSets"];
        statePartitions = new Partition(jsPart);
        }
    
    if (statePartitions != NULL)
        statePartitions->print();
}

void McmcSummary::readTreFile(std::string fn, int bi) {

	// open the file
	std::ifstream fstrm(fn.c_str());
	if (!fstrm)
        Msg::error("Cannot open file \"" + fn + "\"");

	std::string lineString = "";
	int line = 0;
    bool readingTranslateTable = false, readingTree = false;
    std::vector<std::string> translateTokens;
    std::string treeString = "";
    std::map<int,std::string> translateMap;
    int treeCount = 0;
	while( getline(fstrm, lineString).good() )
		{
        //std::cout << line << " -- " << lineString << std::endl;
		std::istringstream linestream(lineString);
		int ch;
		std::string word = "";
		int wordNum = 0;
		std::string cmdString = "";
		do
			{
			word = "";
			linestream >> word;
            if (word == "")
                continue;
            if (word == "translate")
                {
                readingTranslateTable = true;
                }
            else if (word == "tree")
                {
                readingTree = true;
                }
            else
                {
                if (readingTranslateTable == true)
                    {
                    std::vector<std::string> brokenWord = breakString(word);
                    if (hasSemicolon(word) == true)
                        {
                        for (int i=0; i<brokenWord.size(); i++)
                            translateTokens.push_back(brokenWord[i]);
                        readingTranslateTable = false;
                        translateMap = interpretTranslateString(translateTokens);
                        }
                    else
                        {
                        for (int i=0; i<brokenWord.size(); i++)
                            translateTokens.push_back(brokenWord[i]);
                        }
                    }
                else if (readingTree == true)
                    {
                    if (hasSemicolon(word) == true)
                        {
                        treeString += word;
                        readingTree = false;
                        treeCount++;
                        if (treeCount > bi)
                            {
                            std::string newickStr = interpretTreeString(treeString);
                            Tree t(newickStr, translateMap);
                            std::map<RbBitSet,double> parts = t.getPartitions();
                            addPartion(parts);
                            }
                        treeString = "";
                        }
                    else
                        treeString += word;
                    }
                }

			wordNum++;
            } while ( (ch=linestream.get()) != EOF );
						
		line++;
		}

	// close the file
	fstrm.close();
 
    conTree = new Tree(partitions, translateMap);
}

void McmcSummary::readTsvFile(std::string fn, int bi) {

	// open the file
	std::ifstream fstrm(fn.c_str());
	if (!fstrm)
        Msg::error("Cannot open file \"" + fn + "\"");

	std::string lineString = "";
	int line = 0;
	while( getline(fstrm, lineString).good() )
		{
        //std::cout << line << " -- " << lineString << std::endl;
		std::istringstream linestream(lineString);
		int ch;
		std::string word = "";
		int wordNum = 0;
		std::string cmdString = "";
		do
			{
			word = "";
			linestream >> word;
            if (word == "")
                continue;
			if (line == 0)
				{
                ParameterStatistics* s = new ParameterStatistics;
                s->setName(word);
                stats.push_back(s);
				}
			else
				{
                if (line > bi)
                    {
                    double x = atof(word.c_str());
                    ParameterStatistics* s = stats[wordNum];
                    s->addValue(x);
                    }
				}
			wordNum++;
            } while ( (ch=linestream.get()) != EOF );
						
		line++;
		}

	
	// close the file
	fstrm.close();
}
