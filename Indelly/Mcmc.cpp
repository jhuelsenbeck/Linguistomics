#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include "IntVector.hpp"
#include "Mcmc.hpp"
#include "Model.hpp"
#include "Msg.hpp"
#include "Parameter.hpp"
#include "ParameterAlignment.hpp"
#include "RandomVariable.hpp"
#include "Tree.hpp"
#include "UpdateInfo.hpp"
#include "UserSettings.hpp"




Mcmc::Mcmc(Model* m, RandomVariable* r) {

    modelPtr = m;
    rv = r;
    maxLikePrint = 0;
    maxPriorPrint = 0;
    
    UserSettings& settings = UserSettings::userSettings();
    numMcmcCycles = settings.getNumMcmcCycles();
    printFrequency = settings.getPrintFrequency();
    sampleFrequency = settings.getSampleFrequency();

    maxGenPrint = numDigits(numMcmcCycles);
}

void Mcmc::closeOutputFiles(void) {

    parmStrm.close();
    treeStrm.close();

    std::vector<ParameterAlignment*> alns = modelPtr->getAlignments();
    for (int i=0; i<alns.size(); i++)
        algnJsonStrm[i].close();
    delete [] algnJsonStrm;
}

std::string Mcmc::formattedTime(std::chrono::high_resolution_clock::time_point& t1, std::chrono::high_resolution_clock::time_point& t2) {

    std::chrono::duration<double> durationSecs  = std::chrono::duration_cast<std::chrono::seconds>(t2 - t1);
    int s = durationSecs.count();
    int m = s / 60;
    int h = s / 3600;
        
    std::string tStr = "";
    if (h > 0)
        {
        tStr += std::to_string(h) + "h:";
        m -= h * 60;
        s -= h * 60 * 60;
        }
    if (m > 0 || (m == 0 && h > 0))
        {
        tStr += std::to_string(m) + "m:";
        s -= m * 60;
        }
    tStr += std::to_string(s) + "s";
    
    return tStr;
}

int Mcmc::numDigits(double x) {

    if (x < 0.0)
        x = -x;
    return log(x) / log(10.0) + 1;
}

void Mcmc::openOutputFiles(void) {

    // open files for samples
    UserSettings& settings = UserSettings::userSettings();
    std::string outPath = settings.getOutFile();
    std::string parmFileName = outPath + ".tsv";
    std::string treeFileName = outPath + ".tre";

    parmStrm.open( parmFileName.c_str(), std::ios::out );
    if (!parmStrm)
        Msg::error("Cannot open file \"" + parmFileName + "\"");
    treeStrm.open( treeFileName.c_str(), std::ios::out );
    if (!treeStrm)
        Msg::error("Cannot open file \"" + treeFileName + "\"");
  
    std::vector<ParameterAlignment*> alns = modelPtr->getAlignments();
    algnJsonStrm = new std::ofstream[alns.size()];
    for (int i=0; i<alns.size(); i++)
        {
        std::string fn = outPath + "." + alns[i]->getName() + ".aln";
        algnJsonStrm[i].open( fn.c_str(), std::ios::out );
        if (!algnJsonStrm[i])
            Msg::error("Cannot open file \"" + fn + "\"");
        }
}

void Mcmc::print(int gen, double curLnL, double newLnL, double curLnP, double newLnP, bool accept) {

    if (numDigits(curLnL) > maxLikePrint)
        maxLikePrint = numDigits(curLnL);
    if (numDigits(newLnL) > maxLikePrint)
        maxLikePrint = numDigits(newLnL);
    if (numDigits(curLnP) > maxPriorPrint)
        maxPriorPrint = numDigits(curLnP);
    if (numDigits(newLnP) > maxPriorPrint)
        maxPriorPrint = numDigits(newLnP);
        
    std::cout << std::fixed << std::setprecision(3);
    std::cout << "   * ";
    std::cout << std::setw(maxGenPrint) << gen << " --   ";
    std::cout << std::setw(maxLikePrint + 5) << curLnL << " -> " << std::setw(maxLikePrint + 5) << newLnL << "   ";
    std::cout << std::setw(maxPriorPrint + 5) << curLnP << " -> " << std::setw(maxPriorPrint + 5) << newLnP << "   ";
    if (accept == true)
        std::cout << "Accepted update of ";
    else
        std::cout << "Rejected update of ";
    std::cout << modelPtr->getUpdatedParameterName();
    std::cout << " parameter";
    std::string updateType = modelPtr->getLastUpdate();
    if (updateType != "")
        {
        std::cout << " (";
        std::cout << updateType;
        std::cout << " update)";
        }
    std::cout << std::endl;
}

void Mcmc::run(void) {

    UserSettings& settings = UserSettings::userSettings();
    if (settings.getCalculateMarginalLikelihood() == false)
        runPosterior();
    else
        runPathSampling();
}

void Mcmc::runPathSampling(void) {

    std::cout << "   Markov chain Monte Carlo" << std::endl;
    
    // initialize the chain
    openOutputFiles();
    double curLnL = modelPtr->lnLikelihood();
    double curLnP = modelPtr->lnPriorProbability();
    UpdateInfo& updateInfo = UpdateInfo::updateInfo();

    // Metropolis-Hastings algorithm
    std::chrono::high_resolution_clock::time_point start = std::chrono::high_resolution_clock::now();
    for (int n=1; n<=numMcmcCycles; n++)
        {
        // propose a new value for the chain
        double lnProposalRatio = modelPtr->update();
        
        // calculate the likelihood and prior ratios (natural log scale)
        double newLnL = modelPtr->lnLikelihood();
        double lnLikelihoodRatio = newLnL - curLnL;
        double newLnP = modelPtr->lnPriorProbability();
        double lnPriorRatio = newLnP - curLnP;
        
        // accept or reject the state
        bool accept = false;
        if ( log(rv->uniformRv()) < lnLikelihoodRatio + lnPriorRatio + lnProposalRatio )
            accept = true;
        
        // print to the screen to give the user some clue where the chain is
        if (n == 1 || n == numMcmcCycles || n % printFrequency == 0)
            print(n, curLnL, newLnL, curLnP, newLnP, accept);
        
        // update the state of the chain
        if (accept == false)
            {
            modelPtr->reject();
            updateInfo.reject(modelPtr->getLastUpdate());
            }
        else
            {
            modelPtr->accept();
            updateInfo.accept(modelPtr->getLastUpdate());
            curLnL = newLnL;
            curLnP = newLnP;
            }
        
        // sample the chain, printing the current state to files
        if (n == 1 || n == numMcmcCycles || n % sampleFrequency == 0)
            sample(n, curLnL, curLnP);
        }
    std::chrono::high_resolution_clock::time_point stop = std::chrono::high_resolution_clock::now();
        
    // print summary information
    std::cout << std::endl;
    std::cout << "   Markov chain Monte Carlo run summary" << std::endl;
    std::cout << "   * Run time = " << formattedTime(start, stop) << std::endl;
    updateInfo.print();
    std::cout << std::endl;
    
    // clean up
    closeOutputFiles();
}

void Mcmc::runPosterior(void) {

    std::cout << "   Markov chain Monte Carlo" << std::endl;
    
    // initialize the chain
    openOutputFiles();
    double curLnL = modelPtr->lnLikelihood();
    double curLnP = modelPtr->lnPriorProbability();
    UpdateInfo& updateInfo = UpdateInfo::updateInfo();

    // Metropolis-Hastings algorithm
    std::chrono::high_resolution_clock::time_point start = std::chrono::high_resolution_clock::now();
    for (int n=1; n<=numMcmcCycles; n++)
        {
        // propose a new value for the chain
        double lnProposalRatio = modelPtr->update();
        
        // calculate the likelihood and prior ratios (natural log scale)
        double newLnL = modelPtr->lnLikelihood();
        double lnLikelihoodRatio = newLnL - curLnL;
        double newLnP = modelPtr->lnPriorProbability();
        double lnPriorRatio = newLnP - curLnP;
        
        // accept or reject the state
        bool accept = false;
        if ( log(rv->uniformRv()) < lnLikelihoodRatio + lnPriorRatio + lnProposalRatio )
            accept = true;
        
        // print to the screen to give the user some clue where the chain is
        if (n == 1 || n == numMcmcCycles || n % printFrequency == 0)
            print(n, curLnL, newLnL, curLnP, newLnP, accept);
        
        // update the state of the chain
        if (accept == false)
            {
            modelPtr->reject();
            updateInfo.reject(modelPtr->getLastUpdate());
            }
        else
            {
            modelPtr->accept();
            updateInfo.accept(modelPtr->getLastUpdate());
            curLnL = newLnL;
            curLnP = newLnP;
            }
        
        // sample the chain, printing the current state to files
        if (n == 1 || n == numMcmcCycles || n % sampleFrequency == 0)
            sample(n, curLnL, curLnP);
        }
    std::chrono::high_resolution_clock::time_point stop = std::chrono::high_resolution_clock::now();
        
    // print summary information
    std::cout << std::endl;
    std::cout << "   Markov chain Monte Carlo run summary" << std::endl;
    std::cout << "   * Run time = " << formattedTime(start, stop) << std::endl;
    updateInfo.print();
    std::cout << std::endl;
    
    // clean up
    closeOutputFiles();
}

double Mcmc::safeExponentiation(double lnX) {

    if (lnX > 0.0)
        return 1.0;
    else if (lnX < -300.0)
        return 0.0;
    else
        return exp(lnX);
}

void Mcmc::sample(int gen, double lnL, double lnP) {

    Tree* t = modelPtr->getTree();
    double tl = t->getTreeLength();
    std::string ts = t->getNewick();
    
    if (gen == 1)
        {
        std::string parmStr = "";
        parmStr += "Gen\t";
        parmStr += "lnL\t";
        parmStr += "lnP\t";
        parmStr += "TL\t";
        parmStr += modelPtr->getParameterHeader();
        parmStrm << parmStr << std::endl;

        std::vector<std::string>& tn = t->getTaxonNames();
        treeStrm << "begin trees;" << std::endl;
        treeStrm << "   translate" << std::endl;
        for (int i=0; i<tn.size(); i++)
            {
            treeStrm << "      " << i+1 << " " << tn[i];
            if (i+1 != tn.size())
                treeStrm << ",";
            else
                treeStrm << ";";
            treeStrm << std::endl;
            }

        std::vector<ParameterAlignment*> alns = modelPtr->getAlignments();
        for (int i=0; i<alns.size(); i++)
            {
            std::string stateSetsStr = modelPtr->getStateSetsJsonString();
            if (stateSetsStr != "")
                {
                algnJsonStrm[i] << "{" << stateSetsStr;
                algnJsonStrm[i] << ", \"Samples\": [\n";
                }
            else
                {
                algnJsonStrm[i] << "{\"Samples\": [\n";
                }
            algnJsonStrm[i] << alns[i]->getJsonString() << "," << std::endl;
            }
        }
        
    // output to parameter file
    std::string parmStr = "";
    parmStr += std::to_string(gen) + '\t';
    parmStr += std::to_string(lnL) + '\t';
    parmStr += std::to_string(lnP) + '\t';
    parmStr += std::to_string(tl)  + '\t';
    parmStr += modelPtr->getParameterString();
    parmStrm << parmStr << std::endl;
    
    // output to tree file
    treeStrm << "   tree t_" << gen << " = " << ts << ";";
    treeStrm << std::endl;
    
    // output to alignment file
    std::vector<ParameterAlignment*> alns = modelPtr->getAlignments();
    for (int i=0; i<alns.size(); i++)
        {
        algnJsonStrm[i] << alns[i]->getJsonString();
        if (gen == numMcmcCycles)
            algnJsonStrm[i] << "]" << std::endl;
        else
            algnJsonStrm[i] << "," << std::endl;
        }
    
    if (gen == numMcmcCycles)
        {
        treeStrm << "end;" << std::endl;
        }
}
