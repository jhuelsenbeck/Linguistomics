#include <iomanip>
#include <iostream>
#include <string>
#include "Alignment.hpp"
#include "AlignmentProposal.hpp"
#include "Model.hpp"
#include "ParameterAlignment.hpp"
#include "SiteLikelihood.hpp"
#include "TransitionProbabilities.hpp"



ParameterAlignment::ParameterAlignment(RandomVariable* r, Model* m, Alignment* a, std::string n, SiteLikelihood* sl) : Parameter(r, m, n) {
    
    updateChangesRateMatrix = false;
    
    siteProbs = sl;

    std::string name = a->getName();
    const size_t periodIdx = name.rfind('.');
    if (std::string::npos != periodIdx)
        name.erase(periodIdx);
    std::cout << "   * Setting up alignment parameter for word " << name << std::endl;

    // initialize the alignment
    numTaxa = a->getNumTaxa();
    int numSites = a->getNumChar();
    taxonMask = a->getTaxonMask();
    taxonMapKeyCanonical = a->getTaxonMap();
    for (std::map<int,int>::iterator it = taxonMapKeyCanonical.begin(); it != taxonMapKeyCanonical.end(); it++)
        taxonMapKeyAlignment.insert( std::make_pair(it->second,it->first) );
        
    completelySampled = true;
    for (int i=0; i<taxonMask.size(); i++)
        {
        if (taxonMask[i] == false)
            completelySampled = false;
        }
    
    printWidth = numSites * 1.8;
    gapCode = a->getGapCode();
    numStates = a->getNumStates();
    taxonNames = a->getTaxonNames();
    gapCode = a->getGapCode();
    alignment[0].resize(numTaxa);
    for (int i=0; i<numTaxa; i++)
        alignment[0][i].resize(numSites);
    for (int i=0; i<numTaxa; i++)
        {
        for (int j=0; j<numSites; j++)
            alignment[0][i][j] = a->getCharacter(i, j);
        }
    alignment[1] = alignment[0];
    
#   if 0
    std::cout << "alignment from ParameterAlignment " << name << std::endl;
    for (int i=0; i< alignment[0].size(); i++)
        {
        for (int j=0; j<alignment[0][i].size(); j++)
            std::cout << std::setw(2) << alignment[0][i][j] << " ";
        std::cout << std::endl;
        }
#   endif
    
    // initialize the raw sequences
    sequences = a->getRawSequenceMatrix();
 
    tuning = 0.1;
    exponent = 1.5;
    gapPenalty = -10.0;
 
    //print();
}

ParameterAlignment::~ParameterAlignment(void) {

    delete siteProbs;
}

void ParameterAlignment::accept(void) {

    alignment[1] = alignment[0];
}

std::vector<std::vector<int> > ParameterAlignment::getIndelMatrix(void) {

    return getIndelMatrix(0);
}

std::vector<std::vector<int> > ParameterAlignment::getIndelMatrix(int idx) {

    // note that this returns a numSites X numTaxa matrix (i.e.,
    // the rows contain the information for the i-th site while
    // the columns contain the information the j-th taxon)
    size_t nt = alignment[idx].size();
    size_t ns = alignment[idx][0].size();
    std::vector<std::vector<int> > m;
    m.resize(ns);
    for (int i=0; i<ns; i++)
        m[i].resize(nt);
    for (int i=0; i<ns; i++)
        {
        for (int j=0; j<nt; j++)
            {
            if (alignment[idx][j][i] == gapCode)
                m[i][j] = 0;
            else
                m[i][j] = 1;
            }
        }
    return m;
}

std::vector<std::vector<int> > ParameterAlignment::getIndelMatrix(std::vector<std::vector<int> >& aln) {

    size_t nt = aln.size();
    size_t ns = aln[0].size();
    std::vector<std::vector<int> > m;
    m.resize(ns);
    for (int i=0; i<ns; i++)
        m[i].resize(nt);
    for (int i=0; i<ns; i++)
        {
        for (int j=0; j<nt; j++)
            {
            if (aln[j][i] == gapCode)
                m[i][j] = 0;
            else
                m[i][j] = 1;
            }
        }
    return m;
}

std::string ParameterAlignment::getJsonString(void) {

    int longestName = 0;
    for (int i=0; i<taxonNames.size(); i++)
        {
        if (taxonNames[i].length() > longestName)
            longestName = (int)taxonNames[i].length();
        }
        
    std::string jsonStr = "";
    jsonStr += "{\"Name\": \"" + parmName + "\", \"Data\": [\n";
    for (int i=0; i<alignment[0].size(); i++)
        {
        jsonStr += "{\"Taxon\": \"" + taxonNames[i] + "\", ";
        for (int j=0; j<longestName-taxonNames[i].length(); j++)
            jsonStr += " ";
        jsonStr += "\"Segments\": [";
        for (int j=0; j<alignment[0][i].size(); j++)
            {
            int x = alignment[0][i][j];
            if (x == numStates)
                jsonStr += "-1";
            else
                {
                if (x < 10)
                    jsonStr += " ";
                jsonStr += std::to_string(x);
                }
            if (j + 1 != alignment[0][i].size())
                jsonStr += ",";
            }
        jsonStr += "]}";
        if (i + 1 != alignment[0].size())
            jsonStr += ",\n";
        }
    jsonStr += "]}\n";
    
    return jsonStr;
}

RbBitSet ParameterAlignment::getTaxonMask(void) {

    RbBitSet bs(taxonMask);
    return bs;
}

std::string ParameterAlignment::getTaxonMaskString(void) {

    std::string str = "";
    for (int i=0; i<taxonMask.size(); i++)
        {
        if (taxonMask[i] == true)
            str += "1";
        else
            str += "0";
        }
    return str;
}

double ParameterAlignment::lnPriorProbability(void) {

    return 0.0;
}

void ParameterAlignment::print(void) {

    int longestName = 0;
    for (int i=0; i<taxonNames.size(); i++)
        {
        if (taxonNames[i].length() > longestName)
            longestName = (int)taxonNames[i].length();
        }
                
    for (int k=0; k<2; k++)
        {
        std::cout << "alignment[" << k << "]" << std::endl;
        for (int i=0; i<alignment[k].size(); i++)
            {
            std::cout << taxonNames[i] << " ";
            for (int j=0; j<longestName-taxonNames[i].length(); j++)
                std::cout << " ";
            std::map<int,int>::iterator it = taxonMapKeyAlignment.find(i);
            std::cout << std::setw(3) << i << " [ " << it->first << "->" << it->second << "] -- ";
            for (int j=0; j<alignment[k][i].size(); j++)
                std::cout << std::setw(2) << alignment[k][i][j] << " ";
            std::cout << std::endl;
            }
        }

    std::cout << "sequences" << std::endl;
    for (int i=0; i<sequences.size(); i++)
        {
        std::cout << std::setw(3) << i << " -- ";
        for (int j=0; j<sequences[i].size(); j++)
            std::cout << std::setw(2) << sequences[i][j] << " ";
        std::cout << std::endl;
        }
}

void ParameterAlignment::reject(void) {

    alignment[0] = alignment[1];
}

double ParameterAlignment::update(void) {

    lastUpdateType = "alignment for " + parmName;

    // update the alignment
    AlignmentProposal alignmentProposal(this, modelPtr->getTree(), rv, 1.5, -5.0);
    std::vector<std::vector<int> > newAlignment;
    double lnProposalRatio = alignmentProposal.propose(newAlignment, 0.5);
    alignment[0] = newAlignment;

    // set flags indicating the transition probabilities are not affected
    TransitionProbabilities& tip = TransitionProbabilities::transitionProbabilties();
    tip.setNeedsUpdate(false);
    updateChangesTransitionProbabilities = false;

    return lnProposalRatio;
}
