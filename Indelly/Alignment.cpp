#include <iostream>
#include <iomanip>
#include <istream>
#include <sstream>
#include <fstream>
#include "Alignment.hpp"
#include "Msg.hpp"

std::string Alignment::states = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789!$%&`+,./<";
int Alignment::gapCode = (int)states.size();



Alignment::Alignment(std::string fileName, bool dirMode) {

	// open the file
	std::ifstream seqStream(fileName.c_str());
	if (seqStream.is_open() == true)
        {
        if (dirMode == false)
            std::cout << "   * Reading data file \"" << fileName << "\"" << std::endl;
        }
    else
		{
		std::cerr << "Cannot open file \"" + fileName + "\"" << std::endl;
		exit(1);
		}

	std::string linestring = "";
	int line = 0;
	int taxonNum = 0;
	matrix = NULL;
    indelMatrix = NULL;
	numTaxa = numSites = 0;
    dataType = "";
    while ( getline (seqStream, linestring) )
		{
        if (hasBOM(linestring) == true)
            linestring = bomLessString(linestring);
		std::istringstream linestream(linestring);
        //std::cout << line << " -- \"" << linestring << "\"" << std::endl;
		int ch;
		std::string word = "";
		int wordNum = 0;
		int siteNum = 0;
		std::string cmdString = "";
		do
			{
			word = "";
			linestream >> word;
			wordNum++;
            //std::cout << "word:" << wordNum << "\"" << word << "\"" << std::endl;
			if (line == 0)
				{
				// read the number of taxa/chars from the first line
				if (wordNum == 1)
					numTaxa = atoi(word.c_str());
				else if (wordNum == 2)
					numSites = atoi(word.c_str());
                else
                    {
                    if (isInteger(word) == true)
                        {
                        dataType = "syllable";
                        numStates = atoi(word.c_str());
                        }
                    else
                        {
                        std::locale loc;
                        for(auto elem : word)
                           dataType += std::tolower(elem,loc);
                        if (dataType == "dna")
                            Msg::error("Only linguistic data can be input");
                        else if (dataType == "aa")
                            Msg::error("Only linguistic data can be input");
                        else
                            Msg::error("Unknown data type " + word);
                        }
                    }
				if (numTaxa > 0 && numSites > 0 && dataType != "" && matrix == NULL)
					{	
					matrix = new int*[numTaxa];
					matrix[0] = new int[numTaxa * numSites];
					for (size_t i=1; i<numTaxa; i++)
						matrix[i] = matrix[i-1] + numSites;
					for (size_t i=0; i<numTaxa; i++)
						for (size_t j=0; j<numSites; j++)
							matrix[i][j] = 0;
 
                    indelMatrix = new bool*[numTaxa];
                    indelMatrix[0] = new bool[numTaxa * numSites];
                    for (size_t i=1; i<numTaxa; i++)
                        indelMatrix[i] = indelMatrix[i-1] + numSites;
                    for (size_t i=0; i<numTaxa; i++)
                        for (size_t j=0; j<numSites; j++)
                            indelMatrix[i][j] = true;
                            
                    std::string lastPathComponent = fileName;
                    const size_t lastSlashIdx = lastPathComponent.find_last_of("/");
                    if (std::string::npos != lastSlashIdx)
                        lastPathComponent.erase(0, lastSlashIdx + 1);
                    const size_t periodIdx = lastPathComponent.rfind('.');
                    if (std::string::npos != periodIdx)
                        lastPathComponent.erase(periodIdx);

                    if (dataType == "syllable")
                        std::cout << "   * Word alignment \"" << lastPathComponent << "\" has " << numTaxa << " taxa and " << numSites << " syllables" << std::endl;
                    else
                        Msg::error("Unknown data type");
					}
				}
			else
				{
				if (wordNum == 1)
					{
                    taxonNames.push_back(word);
                    taxonNum++;
					}
				else
					{
                    for (int i=0; i<word.length(); i++)
                        {
                        char site = word.at(i);
                        matrix[taxonNum-1][siteNum++] = stateCode(site);
                        }
					}
				}
			} while ( (ch=linestream.get()) != EOF );
			
		line++;
		}	
	
	// close the file
	seqStream.close();
 
    // fill in the indel matrix
    for (int i=0; i<numTaxa; i++)
        for (int j=0; j<numSites; j++)
            {
            if (isIndel(i, j) == true)
                indelMatrix[i][j] = false;
            }
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

char Alignment::getCharFromCode(int code) {

    if (code == states.size())
        return '-';
    return states[code];
}

std::vector<std::vector<int> > Alignment::getIndelMatrix(void) {

    // note that this returns a numSites X numTaxa matrix (i.e.,
    // the rows contain the information for the i-th site while
    // the columns contain the information the j-th taxon)
    std::vector<std::vector<int> > m;
    m.resize(numSites);
    for (int i=0; i<numSites; i++)
        m[i].resize(numTaxa);
    for (int i=0; i<numSites; i++)
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
    for (int j=0; j<numSites; j++)
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

void Alignment::print(void) {

	int** x = matrix;
		
	std::cout << "        ";
	for (size_t i=0; i<numTaxa; i++)
		std::cout << std::setw(3) << i;
	std::cout << '\n';
	std::cout << "------------------------";
	for (size_t i=0; i<numTaxa; i++)
		std::cout << "---";
	std::cout << '\n';
	for (size_t j=0; j<numSites; j++)
		{
		std::cout << std::setw(4) << j+1 << " -- ";
		for (size_t i=0; i<numTaxa; i++)
            std::cout << std::setw(3) << getCharFromCode( x[i][j] );
		std::cout << '\n';
		}
    std::cout << std::endl;
}

void Alignment::printCode(void) {
    
    std::cout << name << std::endl;
    std::cout << "numTaxa = " << numTaxa << ";" << std::endl;
    std::cout << "numChar = " << numSites << ";" << std::endl;
    std::cout << "aln = new char[numTaxa][numChar];" << std::endl;
    for (int i=0; i<numTaxa; i++)
        {
        for (int j=0; j<numSites; j++)
            {
            int charCode = getCharacter(i, j);
            std::cout << "aln[" << i << "][" << j << "] = '" << getCharFromCode(charCode) << "';" << std::endl;
            }
        }
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
    for (size_t j=0; j<numSites; j++)
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

int Alignment::stateCode(char s) {

    if (s == '?')
        Msg::error("Missing data (?) is not allowed in the data matrix");
    if (s == '-')
        return (int)states.size();
        
    for (int i=0; i<numStates; i++)
        {
        char c = states[i];
        if (c == s)
            return i;
        }
        
    Msg::error("Unidentified state " + std::to_string(s));
    return -1;
}
