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



McmcSummary::McmcSummary(RandomVariable* r):
    conTree(NULL),
    rv(r)
{
}

McmcSummary::~McmcSummary(void) {

    for (int i=0; i<stats.size(); i++)
        delete stats[i];
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

    // print summary of csv file
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
}

void McmcSummary::output(UserSettings& settings) {
    auto& file = *new std::ofstream(settings.getPath() + "/consensus.tre" , std::ios::out);
    file << "#NEXUS\n\nbegin taxa;" << std::endl;

    file << "dimensions=" << conTree->numTaxa << ";" << std::endl;
    file << "end;\nbegin trees;\ntree TREE1 = ";
    file << conTree->getNewick(4) << ";"  << std::endl;
    file << "end;" << std::endl;
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
        Msg::error("Error parsing JSON file at byte " + std::to_string(ex.byte));
        }
        
    std::string cognateName = getCognateName(fn);
        
    auto it = j.find("PartitionSets");
    if (it == j.end())
        Msg::error("Could not find partition set in the JSON file");
    nlohmann::json jsPart = j["PartitionSets"];
    Partition* part = new Partition(jsPart);

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
