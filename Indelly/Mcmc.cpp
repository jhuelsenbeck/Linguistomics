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
#include "Probability.hpp"
#include "RandomVariable.hpp"
#include "SteppingStones.h"
#include "Tree.hpp"
#include "UpdateInfo.hpp"
#include "UserSettings.hpp"




Mcmc::Mcmc(Model* m, RandomVariable* r) {

    modelPtr      = m;
    rv            = r;
    maxLikePrint  = 0;
    maxPriorPrint = 0;
    
    UserSettings& settings = UserSettings::userSettings();
    numMcmcCycles        = settings.getNumMcmcCycles();
    printFrequency       = settings.getPrintFrequency();
    sampleFrequency      = settings.getSampleFrequency();
    preburninLength      = settings.getPreburninLength();
    numTunes             = settings.getNumTunes();
    tuneLength           = settings.getTuneLength();
    burninLength         = settings.getBurninLength();
    sampleLength         = settings.getSampleLength();
    stoneSampleFrequency = settings.getSampleToStoneFrequency();
    maxGenPrint          = numDigits(settings.getNumMcmcCycles());
}

std::vector<double> Mcmc::calculatePowers(int numStones, double alpha, double beta) {

    int ns = numStones - 1;
    std::vector<double> pwrs;
    double intervalProb = (double)1.0 / ns;
    pwrs.push_back(1.0);
    for (int i=ns-1; i>0; i--)
        pwrs.push_back( Probability::Beta::quantile(alpha, beta, i * intervalProb) );
    pwrs.push_back(0.0);
#   if 0
    for (int i=0; i<pwrs.size(); i++)
        std::cout << i+1 << " -- " << std::fixed << std::setprecision(15) << pwrs[i] << std::endl;
#   endif
    return pwrs;
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
    int s = (int)durationSecs.count();
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

void Mcmc::initialize(void) {

    openOutputFiles();
    numParmValues = modelPtr->getNumParameterValues();
    parmValues = new double[numParmValues];
    for (int i=0; i<numParmValues; i++)
        parmValues[i] = 0.0;
}

int Mcmc::numDigits(double x) {

    if (x < 0.0)
        x = -x;
    return (int)(log(x) / log(10.0)) + 1;
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
    auto outsize = alns.size();
    algnJsonStrm = new std::ofstream[outsize];
    for (int i=0; i< outsize; i++)
        {
        std::string fn = outPath + "." + alns[i]->getName() + ".aln";
        algnJsonStrm[i].open( fn.c_str(), std::ios::out );
        if (!algnJsonStrm[i])
            Msg::error("Cannot open file \"" + fn + "\"");
        }
}

int Mcmc::phaseLength(std::string phs) {

    if (phs == "preburn")
        return preburninLength;
    else if (phs == "tune")
        return numTunes * tuneLength;
    else if (phs == "burn")
        return burninLength;
    else
        return sampleLength;
}

void Mcmc::print(int gen, double curLnL, double newLnL, double curLnP, double newLnP, bool accept, std::chrono::high_resolution_clock::time_point& t1, std::chrono::high_resolution_clock::time_point& t2) {

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

    // estimate time remaining based on time to this point
    std::chrono::duration<double> durationSecs  = std::chrono::duration_cast<std::chrono::seconds>(t1 - t2);
    double timePerCycle = (double)durationSecs.count() / gen;
    if (timePerCycle == 0)
        timePerCycle = 1.0 / printFrequency;
    int s = (int)((numMcmcCycles - gen) * timePerCycle);
    int m = s / 60;
    int h = s / 3600;

    if (h > 0)
        {
        std::cout << h;
        std::cout << ":";
        m -= h * 60;
        s -= h * 60 * 60;
        }
    if (m > 0 || (m == 0 && h > 0))
        {
        if (m < 10)
            std::cout << "0";
        std::cout << m;
        std::cout << ":";
        s -= m * 60;
        }
    if (s < 10) 
        std::cout << "0";
    std::cout << s;
    std::cout << " remaining  ";

    //std::cout << std::setw(maxPriorPrint + 5) << curLnP << " -> " << std::setw(maxPriorPrint + 5) << newLnP << "   ";
    if (accept == true)
        std::cout << "Accepted ";
    else
        std::cout << "Rejected ";
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

    // initialize the chain
    initialize();
    
    // get powers and initialize vector with phases
    int numStones = 127;
    double alpha = 0.3;
    double beta = 1.0;
    std::vector<double> powers = calculatePowers(numStones, alpha, beta);
    SteppingStones samples(powers);
    std::vector<std::string> mcmcPhases = { "preburn", "tune", "burn", "sample" };
    
    modelPtr->setUpdateLikelihood();
    double curLnL = modelPtr->lnLikelihood();
    double curLnP = modelPtr->lnPriorProbability();
    UpdateInfo& updateInfo = UpdateInfo::updateInfo();

    // Metropolis-Hastings algorithm
    auto start = std::chrono::high_resolution_clock::now();

    int iteration = 0;
//    int numIterations = (preburninLength + (numTunes*tuneLength) + burninLength + sampleLength) * (int)powers.size();

    for (int powIdx=0; powIdx<powers.size(); powIdx++)
        {
        // set the power
        double power = powers[powIdx];
        for (std::string phase : mcmcPhases)
            {
            int phsLength = phaseLength(phase);

            for (int n=1; n<=phsLength; n++)
                {
                iteration++;
                
                // propose a new value for the chain
                double lnProposalRatio = modelPtr->update();
                
                // calculate the likelihood and prior ratios (natural log scale)
                double newLnL = modelPtr->lnLikelihood();
                double lnLikelihoodRatio = (newLnL - curLnL) * power;
                double newLnP = modelPtr->lnPriorProbability();
                double lnPriorRatio = newLnP - curLnP;
                
                // accept or reject the state
                bool accept = false;
                if ( log(rv->uniformRv()) < lnLikelihoodRatio + lnPriorRatio + lnProposalRatio )
                    accept = true;
                
                // print to the screen to give the user some clue where the chain is
                if (n == 1 || n == numMcmcCycles || n % printFrequency == 0)
                    {
                    auto timePt = std::chrono::high_resolution_clock::now();
                    print(n, curLnL, newLnL, curLnP, newLnP, accept, timePt, start);
                    }
                
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

                 if ( (n % stoneSampleFrequency == 0  || n == sampleLength) && phase == "sample" )
                    samples.addSample(powIdx, curLnL);

                if ( n % tuneLength == 0 && phase == "tune" )
                    {
                    }
                }
            
            }
            
        }
    auto stop = std::chrono::high_resolution_clock::now();
        
    // print summary information
    std::cout << std::endl;
    std::cout << "   Markov chain Monte Carlo run summary" << std::endl;
    std::cout << "   * Run time = " << formattedTime(start, stop) << std::endl;
    updateInfo.print();
    std::cout << std::endl;

    double marginalLnL = samples.marginalLikelihood();
    std::cout << "   * Marginal likelihood = " << marginalLnL << std::endl;

    // clean up
    closeOutputFiles();
    delete [] parmValues;
}

void Mcmc::runPosterior(void) {

    std::cout << "   Markov chain Monte Carlo" << std::endl;
    
    // initialize the chain
    initialize();
    
    modelPtr->setUpdateLikelihood();
    double curLnL = modelPtr->lnLikelihood();
    double curLnP = modelPtr->lnPriorProbability();
    UpdateInfo& updateInfo = UpdateInfo::updateInfo();

    // Metropolis-Hastings algorithm
    auto start = std::chrono::high_resolution_clock::now();
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
            {
            auto timePt = std::chrono::high_resolution_clock::now();
            print(n, curLnL, newLnL, curLnP, newLnP, accept, timePt, start);
            }
        
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
    auto stop = std::chrono::high_resolution_clock::now();
        
    // print summary information
    std::cout << std::endl;
    std::cout << "   Markov chain Monte Carlo run summary" << std::endl;
    std::cout << "   * Run time = " << formattedTime(start, stop) << std::endl;
    updateInfo.print();
    std::cout << std::endl;
    
    // clean up
    closeOutputFiles();
    delete [] parmValues;
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
    std::string ts = t->getNewick(6);
    
    if (gen == 1)
        {
        parmStrm << "Gen\t";
        parmStrm << "lnL\t";
        parmStrm << "lnP\t";
        parmStrm << "TL\t";
        parmStrm << modelPtr->getParameterHeader();
        parmStrm << std::endl;

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
#   if 0
    std::string parmStr = "";
    parmStr += std::to_string(gen) + '\t';
    parmStr += std::to_string(lnL) + '\t';
    parmStr += std::to_string(lnP) + '\t';
    parmStr += std::to_string(tl)  + '\t';
    parmStr += modelPtr->getParameterString();
    parmStrm << parmStr << std::endl;
#   else
    std::cout << std::fixed << std::setprecision(8);
    modelPtr->fillParameterValues(parmValues, numParmValues);
    parmStrm << gen << '\t';
    parmStrm << lnL << '\t';
    parmStrm << lnP << '\t';
    parmStrm << tl << '\t';
    for (int i=0; i<numParmValues; i++)
        parmStrm << parmValues[i] << '\t';
    parmStrm << std::endl;
#   endif
    
    // output to tree file
#   if 0
    treeStrm << "   tree t_" << gen << " = " << ts << ";";
    treeStrm << std::endl;
#   else
    treeStrm << "   tree t_" << gen << " = ";
    t->newickStream(treeStrm, 6);
    treeStrm  << ts << ";" << std::endl;
#   endif
    
    // output to alignment file
    std::vector<ParameterAlignment*> alns = modelPtr->getAlignments();
    for (int i=0; i<alns.size(); i++)
        {
        //algnJsonStrm[i] << alns[i]->getJsonString();
        alns[i]->jsonStream(algnJsonStrm[i]);
        if (gen == numMcmcCycles)
            algnJsonStrm[i] << "]\n}" << std::endl;
        else
            algnJsonStrm[i] << "," << std::endl;
        }
    
    if (gen == numMcmcCycles)
        treeStrm << "end;" << std::endl;
}
