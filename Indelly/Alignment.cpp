#include <iostream>
#include <iomanip>
#include <istream>
#include <sstream>
#include <fstream>
#include "Alignment.hpp"
#include "Msg.hpp"



Alignment::Alignment(nlohmann::json& j, int ns, std::vector<std::string> canonicalTaxonList) {
    
    //std::cout << j.dump() << std::endl;
    
    matrix = NULL;
    indelMatrix = NULL;
    
    // set the state information
    numStates = ns;
    gapCode = numStates;
    if (gapCode <= 1)
        Msg::error("There must be at least two states in the model");

    // read the word name
    auto it = j.find("Name");
    if (it == j.end())
        Msg::error("Could not find word name in the JSON object");
    name = j["Name"];
    
    // read the data
    it = j.find("Data");
    if (it == j.end())
        Msg::error("Could not find word data in the JSON object");
    numTaxa = (int)j["Data"].size();
    if (numTaxa <= 0)
        Msg::error("Must have at least one taxon in the word");
    numChar = 0;
    for (int i=0; i<numTaxa; i++)
        {
        nlohmann::json jw = j["Data"][i];
        
        // check that there is taxon name information
        it = jw.find("Taxon");
        if (it == jw.end())
            Msg::error("Could not find taxon name in the JSON object");
        std::string tName = jw["Taxon"];
        taxonNames.push_back(tName);
        
        // check that there is segment information
        it = jw.find("Segments");
        if (it == jw.end())
            Msg::error("Could not find segment data in the JSON object");
        std::vector<int> segInfo = jw["Segments"];
        
        // check that the number of word segments is the same for each taxon
        if (i == 0)
            {
            numChar = (int)segInfo.size();
            if (numChar <= 0)
                Msg::error("Must have at least one segment in the word");

            if (matrix == NULL)
                {
                matrix = new int*[numTaxa];
                matrix[0] = new int[numTaxa * numChar];
                for (size_t i=1; i<numTaxa; i++)
                    matrix[i] = matrix[i-1] + numChar;
                for (size_t i=0; i<numTaxa; i++)
                    for (size_t j=0; j<numChar; j++)
                        matrix[i][j] = 0;
                }
            if (indelMatrix == NULL)
                {
                indelMatrix = new bool*[numTaxa];
                indelMatrix[0] = new bool[numTaxa * numChar];
                for (size_t i=1; i<numTaxa; i++)
                    indelMatrix[i] = indelMatrix[i-1] + numChar;
                for (size_t i=0; i<numTaxa; i++)
                    for (size_t j=0; j<numChar; j++)
                        indelMatrix[i][j] = true;
                }
            }
        else
            {
            if (segInfo.size() != numChar)
                Msg::error("Inconsistent segment lengths for word " + name);
            }
            
        //read the segment information
        for (int j=0; j<segInfo.size(); j++)
            {
            if (segInfo[j] == -1)
                matrix[i][j] = gapCode;
            else
                matrix[i][j] = segInfo[j];
            }
        }
        
    // fill in the indel matrix
    for (int i=0; i<numTaxa; i++)
        {
        for (int j=0; j<numChar; j++)
            {
            if (isIndel(i, j) == true)
                indelMatrix[i][j] = false;
            }
        }
            
    std::cout << "   * Word alignment \"" << name << "\" has " << numTaxa << " taxa and " << numChar << " syllables" << std::endl;
}


Alignment::~Alignment(void) {

	delete [] matrix[0];
	delete [] matrix;
    delete [] indelMatrix[0];
    delete [] indelMatrix;
}

std::string Alignment::bomLessString(std::string& str) {

    std::stringstream is(str);
    std::string rStr = "";
    char c;
    while (is.get(c))
        {
        if (c != '\xEF' && c != '\xBB' && c != '\xBF')
            rStr += c;
        }
    return rStr;
}

int Alignment::getCharacter(size_t i, size_t j) {

    return matrix[i][j];
}

std::vector<std::string> Alignment::getTaxonNames(void) {

    return taxonNames;
}

std::vector<std::vector<int> > Alignment::getIndelMatrix(void) {

    // note that this returns a numSites X numTaxa matrix (i.e.,
    // the rows contain the information for the i-th site while
    // the columns contain the information the j-th taxon)
    std::vector<std::vector<int> > m;
    m.resize(numChar);
    for (int i=0; i<numChar; i++)
        m[i].resize(numTaxa);
    for (int i=0; i<numChar; i++)
        {
        for (int j=0; j<numTaxa; j++)
            {
            if (isIndel(j, i) == true)
                m[i][j] = 0;
            else
                m[i][j] = 1;
            }
        }
    return m;
}

std::vector<int> Alignment::getRawSequence(int txnIdx) {

    std::vector<int> seq;
    for (int j=0; j<numChar; j++)
        {
        if (isIndel(txnIdx, j) == false)
            {
            seq.push_back( matrix[txnIdx][j] );
            }
        }
    return seq;
}

std::vector<std::vector<int> > Alignment::getRawSequenceMatrix(void) {

    std::vector<std::vector<int> > seqMat;
    
    seqMat.resize(numTaxa);
    for (int i=0; i<numTaxa; i++)
        {
        seqMat[i] = getRawSequence(i);
        }
        
    return seqMat;
}

bool Alignment::hasBOM(std::string& str) {

    std::stringstream is(str);

    // read the first byte
    char const c0 = is.get();
    if (c0 != '\xEF')
        {
        is.putback(c0);
        return false;
        }

    // read the second byte
    char const c1 = is.get();
    if (c1 != '\xBB')
        {
        is.putback(c1);
        is.putback(c0);
        return false;
        }

    // peek the third byte
    char const c2 = is.peek();
    if (c2 != '\xBF')
        {
        is.putback(c1);
        is.putback(c0);
        return false;
        }

    return true;
}

bool Alignment::isIndel(size_t i, size_t j) {

    if (matrix[i][j] == gapCode)
        return true;
    return false;
}

bool Alignment::isInteger(const std::string& str) {

    return str.find_first_not_of("0123456789") == std::string::npos;
}

void Alignment::listTaxa(void) {

	int i = 1;
	for (std::vector<std::string>::iterator p=taxonNames.begin(); p != taxonNames.end(); p++)
		std::cout << std::setw(4) << i++ << " -- " << (*p) << '\n';
}

std::string Alignment::getTaxonName(int i) {

	return taxonNames[i];
}

int Alignment::getTaxonIndex(std::string ns) {

	int taxonIndex = -1;
	int i = 0;
	for (std::vector<std::string>::iterator p=taxonNames.begin(); p != taxonNames.end(); p++)
		{
		if ( (*p) == ns )
			{
			taxonIndex = i;
			break;
			}
		i++;
		}
	return taxonIndex;
}

int Alignment::numCompleteTaxa(void) {

    int n = 0;
    for (int i=0; i<numTaxa; i++)
        {
        bool hasNongap = false;
        for (int j=0; j<numChar; j++)
            {
            if (isIndel(i,j) == false)
                hasNongap = true;
            }
        if (hasNongap == true)
            n++;
        }
    
    return n;
}

void Alignment::print(void) {

	int** x = matrix;
		
    std::cout << "Name: " << name << std::endl;
    std::cout << "Taxa: ( ";
    for (int i=0; i<taxonNames.size(); i++)
        std::cout << taxonNames[i] << " ";
    std::cout << ")" << std::endl;
    std::cout << "Name: " << name << std::endl;
	std::cout << "        ";
	for (size_t i=0; i<numTaxa; i++)
		std::cout << std::setw(3) << i;
	std::cout << '\n';
	std::cout << "------------------------";
	for (size_t i=0; i<numTaxa; i++)
		std::cout << "---";
	std::cout << '\n';
	for (size_t j=0; j<numChar; j++)
		{
		std::cout << std::setw(4) << j+1 << " -- ";
		for (size_t i=0; i<numTaxa; i++)
            {
            if (x[i][j] == gapCode)
                std::cout << std::setw(3) << "-";
            else
                std::cout << std::setw(3) << x[i][j];
            }
		std::cout << '\n';
		}
    std::cout << std::endl;
}

void Alignment::printIndels(void) {

    bool** x = indelMatrix;
        
    std::cout << "        ";
    for (size_t i=0; i<numTaxa; i++)
        std::cout << std::setw(3) << i;
    std::cout << '\n';
    std::cout << "------------------------";
    for (size_t i=0; i<numTaxa; i++)
        std::cout << "---";
    std::cout << '\n';
    for (size_t j=0; j<numChar; j++)
        {
        std::cout << std::setw(4) << j+1 << " -- ";
        for (size_t i=0; i<numTaxa; i++)
            {
            std::cout << std::setw(3) << x[i][j];
            }
        std::cout << '\n';
        }
    std::cout << std::endl;
}

