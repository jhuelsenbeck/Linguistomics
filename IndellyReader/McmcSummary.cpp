#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include "Alignment.hpp"
#include "AlignmentDistribution.hpp"
#include "json.hpp"
#include "McmcSummary.hpp"
#include "Msg.hpp"
#include "Partition.hpp"
#include "RandomVariable.hpp"
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

bool McmcSummary::hasSemicolon(std::string str) {

    std::size_t found = str.find(';');
    if (found != std::string::npos)
        return true;
    return false;
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

void McmcSummary::output(UserSettings& settings) {

    nlohmann::json j = nlohmann::json::object();

    // output the consensus tree
    j["consensus"]["tree"] = conTree->getNewick(4);
    j["consensus"]["n_taxa"] = conTree->numTaxa;
    
    // output information on mean and credible interval for all real-valued parameters
    auto jStats = nlohmann::json::array();
    for (int i=0; i<stats.size(); i++)
        {
        nlohmann::json cogStats = stats[i]->toJson();
        jStats.push_back(cogStats);
        }
    j["stats"] = jStats;
        
    // output the credible set for all alignments
    auto sampledAlignments = nlohmann::json::array();
    for (int i=0; i<alignments.size(); i++)
        {
        nlohmann::json cogAlns = alignments[i]->toJson(0.95);
        sampledAlignments.push_back(cogAlns);
        }
    j["cog_alns"] = sampledAlignments;

    // output json representation to a file
    auto& file = *new std::ofstream(settings.getPath() + "/summary.json" , std::ios::out);
    file << j;
    std::cout << "JSON file written to: ";
    std::cout << settings.getPath() << "summary.json" << std::endl;
}
    
void McmcSummary::readAlnFile(std::string fn, int /*bi*/) {

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
        Alignment* aln = new Alignment( js[i]["Data"] );
        dist->addAlignment(aln);
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

void McmcSummary::readTreFile(std::string fn, int /*bi*/) {

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
                        std::string newickStr = interpretTreeString(treeString);
                        //std::cout << "Tree: \"" << newickStr << "\"" << std::endl;
                        Tree t(newickStr, translateMap);
                        //t.print();
                        std::map<RbBitSet,double> parts = t.getPartitions();
                        addPartion(parts);
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
