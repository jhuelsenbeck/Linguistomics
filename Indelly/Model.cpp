#include <future>
#include <fstream>
#include <iostream>
#include <list>
#include <mutex>
#include "Alignment.hpp"
#include "EigenSystem.hpp"
#include "LikelihoodTkf.hpp"
#include "Model.hpp"
#include "Msg.hpp"
#include "ParameterAlignment.hpp"
#include "ParameterEquilibirumFrequencies.hpp"
#include "ParameterExchangabilityRates.hpp"
#include "ParameterIndelRates.hpp"
#include "ParameterTree.hpp"
#include "RandomVariable.hpp"
#include "RateMatrixHelper.hpp"
#include "ThreadPool.hpp"
#include "TransitionProbabilities.hpp"
#include "Tree.hpp"
#include "UserSettings.hpp"




Model::Model(RandomVariable* r) {

    std::cout << "   Model" << std::endl;
    rv = r;

    // read the json file
    nlohmann::json j = parseJsonFile();

    // initialize alignments
    std::vector<Alignment*> wordAlignments = initializeAlignments(j);
    
    // initialize parameters of model
    initializeParameters(wordAlignments, j);
    
    // initialize transition probabilities
    initializeTransitionProbabilities(wordAlignments[0]->getNumStates(), j);
        
    std::cout << std::endl;
}

Model::~Model(void) {

    delete [] threadLnL;
    for (int i=0; i<parameters.size(); i++)
        delete parameters[i];
}

void Model::accept(void) {

    parameters[updatedParameterIdx]->accept();
}

ParameterAlignment* Model::getAlignment(int idx) {

    int n = 0;
    for (int i=0; i<parameters.size(); i++)
        {
        ParameterAlignment* p = dynamic_cast<ParameterAlignment*>(parameters[i]);
        if (p != NULL)
            {
            if (n == idx)
                return p;
            n++;
            }
        }
    return NULL;
}

std::vector<ParameterAlignment*> Model::getAlignments(void) {

    std::vector<ParameterAlignment*> alns;
    for (int i=0; i<parameters.size(); i++)
        {
        ParameterAlignment* p = dynamic_cast<ParameterAlignment*>(parameters[i]);
        if (p != NULL)
            alns.push_back(p);
        }
    return alns;
}

double Model::getDeletionRate(void) {

    for (int i=0; i<parameters.size(); i++)
        {
        ParameterIndelRates* p = dynamic_cast<ParameterIndelRates*>(parameters[i]);
        if (p != NULL)
            {
            return p->getDeletionRate();
            }
        }
    return 0.0;
}

std::vector<double>& Model::getEquilibriumFrequencies(void) {

    ParameterEquilibirumFrequencies* p = NULL;
    for (int i=0; i<parameters.size(); i++)
        {
        p = dynamic_cast<ParameterEquilibirumFrequencies*>(parameters[i]);
        if (p != NULL)
            break;
        }
    return p->getValue();
}

std::vector<double>& Model::getExchangabilityRates(void) {

    ParameterExchangabilityRates* p = NULL;
    for (int i=0; i<parameters.size(); i++)
        {
        p = dynamic_cast<ParameterExchangabilityRates*>(parameters[i]);
        if (p != NULL)
            break;
        }
    return p->getValue();
}

double Model::getInsertionRate(void) {

    for (int i=0; i<parameters.size(); i++)
        {
        ParameterIndelRates* p = dynamic_cast<ParameterIndelRates*>(parameters[i]);
        if (p != NULL)
            {
            return p->getInsertionRate();
            }
        }
    return 0.0;
}

int Model::getNumAlignments(void) {

    int n = 0;
    for (int i=0; i<parameters.size(); i++)
        {
        ParameterAlignment* p = dynamic_cast<ParameterAlignment*>(parameters[i]);
        if (p != NULL)
            n++;
        }
    return n;
}

std::string Model::getLastUpdate(void) {

    std::string str = parameters[updatedParameterIdx]->getLastUpdate();
    return str;
}

std::string Model::getParameterHeader(void) {

    std::string str = "";
    for (int i=0; i<parameters.size(); i++)
        {
        ParameterAlignment* pA = dynamic_cast<ParameterAlignment*>(parameters[i]);
        ParameterTree* pT = dynamic_cast<ParameterTree*>(parameters[i]);
        if (pA == NULL && pT == NULL)
            str += parameters[i]->getHeader();
        }
    return str;
}

std::string Model::getParameterString(void) {

    std::string str = "";
    for (int i=0; i<parameters.size(); i++)
        {
        ParameterAlignment* pA = dynamic_cast<ParameterAlignment*>(parameters[i]);
        ParameterTree* pT = dynamic_cast<ParameterTree*>(parameters[i]);
        if (pA == NULL && pT == NULL)
            str += parameters[i]->getString();
        }
    return str;
}

Tree* Model::getTree(void) {

    Tree* t = NULL;
    for (int i=0; i<parameters.size(); i++)
        {
        if (parameters[i]->getName() == "tree")
            {
            ParameterTree* pt = dynamic_cast<ParameterTree*>(parameters[i]);
            return pt->getActiveTree();
            }
        }
    return t;
}

std::string Model::getUpdatedParameterName(void) {

    return parameters[updatedParameterIdx]->getName();
}

std::vector<Alignment*> Model::initializeAlignments(nlohmann::json& j) {

    // get the list of taxa
    auto it = j.find("Taxa");
    if (it == j.end())
        Msg::error("Could not find list of taxa in the JSON file");
    size_t numTaxa = j["Taxa"].size();
    std::vector<std::string> canonicalTaxonList = j["Taxa"];
    std::cout << "   * Found " << canonicalTaxonList.size() << " taxa in JSON file" << std::endl;
    
    // check that there are words in the json object
    it = j.find("Words");
    if (it == j.end())
        Msg::error("Could not find word information in the JSON file");
    size_t numWords = j["Words"].size();
    std::cout << "   * Found " << numWords << " words in JSON file" << std::endl;
    
    // check that the json object has the number of states
    it = j.find("NumberOfStates");
    if (it == j.end())
        Msg::error("Could not the number of states in the JSON file");
    int numStates = j["NumberOfStates"];
    if (numStates <= 1)
        Msg::error("There must be at least two states in the model");
    std::cout << "   * The substitution model has " << numStates << " states" << std::endl;
    
    // read each word
    std::vector<Alignment*> words;
    std::vector<std::string> rejectedWords;
    std::vector<std::string> taxonNames;
    for (size_t i=0; i<numWords; i++)
        {
        // get the json object for the word
        nlohmann::json jw = j["Words"][i];

        // instantiate the word alignment from the json object
        Alignment* aln = new Alignment(jw, numStates);
                
        // check if we have problems, deleting the alignment if we do
        std::vector<std::string> alnTaxonNames = aln->getTaxonNames();
        if (canonicalTaxonList.size() != aln->numCompleteTaxa())
            {
            rejectedWords.push_back(aln->getName());
            delete aln;
            continue;
            }
            
        if (taxonNames.size() == 0)
            taxonNames = aln->getTaxonNames();
            
        if (taxonNames.size() != alnTaxonNames.size())
            Msg::error("List of taxa is inconsistent across words");
        for (int j=0; j<taxonNames.size(); j++)
            {
            if (taxonNames[j] != alnTaxonNames[j])
                Msg::error("List of taxa is inconsistent across words");
            }
            
        // make certain that each taxon in the alignment is in the canonical list of taxa
        for (int j=0; j<alnTaxonNames.size(); j++)
            {
            if (std::find(canonicalTaxonList.begin(), canonicalTaxonList.end(), alnTaxonNames[j]) == canonicalTaxonList.end())
                Msg::error("Could not find taxon " + alnTaxonNames[j] + " in canonical taxon list for word " + aln->getName());
            }
            
        words.push_back(aln);
        }

    if (words.size() < 1)
        Msg::error("No word alignments were read!");
        
    std::cout << "   * Number of words                       = " << words.size() << std::endl;
    std::cout << "   * Number of taxa in each word alignment = " << taxonNames.size() << std::endl;
    if (rejectedWords.size() > 0)
        {
        std::cout << "   * These words were not included because at least one " << std::endl;
        std::cout << "     taxon had no word segments: ";
        int cnt = 28;
        for (int i=0; i<rejectedWords.size(); i++)
            {
            std::cout << rejectedWords[i];
            if (i+1 != rejectedWords.size())
                std::cout << ", ";
            cnt += rejectedWords[i].length();
            if (cnt > 40)
                {
                std::cout << std::endl << "     ";
                cnt = 0;
                }
            }
        std::cout << std::endl;
        }
    
#   if 0
    for (int i=0; i<words.size(); i++)
        words[i]->print();
#   endif

    return words;
}

void Model::initializeParameters(std::vector<Alignment*>& wordAlignments, nlohmann::json& j) {

    UserSettings& settings = UserSettings::userSettings();

    // attempt to find the initial tree in the json object
    std::string treeStr = "";
    auto it = j.find("Tree");
    if (it != j.end())
        treeStr = j["Tree"];
        
    // find the number of states in the json object
    it = j.find("NumberOfStates");
    if (it == j.end())
        Msg::error("Could not find the number of states in the substitution model");
    int numStates = j["NumberOfStates"];

    // initialize a few important parameters
    substitutionModel = settings.getSubstitutionModel();
    std::vector<std::string> taxonNames = wordAlignments[0]->getTaxonNames();
    
    // if the model is a custom one, make certain there is a partition of states to read
    if (substitutionModel == custom)
        {
        it = j.find("PartitionSets");
        if (it == j.end())
            Msg::error("Could not finda  partition of the substitution model states");
        nlohmann::json jsonPart = j["PartitionSets"];
            
        RateMatrixHelper& helper = RateMatrixHelper::rateMatrixHelper();
        helper.initialize(numStates, jsonPart);
        helper.print();
        }
    
    // set up the tree parameter
    Parameter* pTree = new ParameterTree(rv, this, treeStr, taxonNames, settings.getInverseTreeLength());
    pTree->setProposalProbability(10.0);
    parameters.push_back(pTree);
    
    // set up the indel parameter
    Parameter* pIndel = new ParameterIndelRates(rv, this, "indel", 7.0, 100.0, 100.0);
    pIndel->setProposalProbability(1.0);
    parameters.push_back(pIndel);

    // set up the alignment parameter(s)
    for (int i=0; i<wordAlignments.size(); i++)
        {
        std::string alnName = wordAlignments[i]->getName();
        Parameter* pAlign = new ParameterAlignment(rv, this, wordAlignments[i], alnName);
        pAlign->setProposalProbability(1.0);
        parameters.push_back(pAlign);
        wordParameterAlignments.push_back( dynamic_cast<ParameterAlignment*>(pAlign) );
        }
        
    // check for consistency between alignment(s) and tree
    std::vector<std::string> treeTaxonNames = getTree()->getTaxonNames();
    if (treeTaxonNames.size() != taxonNames.size())
        Msg::error("Mismatch between the size of the starting tree and the alignments");
    for (int i=0; i<treeTaxonNames.size(); i++)
        {
        if (treeTaxonNames[i] != taxonNames[i])
            Msg::error("Taxon names in tree do not match names in the alignments");
        }
    if (getTree()->isBinary() == false)
        Msg::error("The initial tree is not binary");
    
    if (substitutionModel == gtr)
        {
        // set up exchangability parameters
        Parameter* pExchange = new ParameterExchangabilityRates(rv, this, "exchangability", numStates);
        pExchange->setProposalProbability(50.0);
        parameters.push_back(pExchange);
        
        // set up equilibrium frequencies
        Parameter* pEquil = new ParameterEquilibirumFrequencies(rv, this, "stationary", numStates);
        pEquil->setProposalProbability(50.0);
        parameters.push_back(pEquil);
        }
    else if (substitutionModel == custom)
        {
        // set up exchangability parameters
        RateMatrixHelper& helper = RateMatrixHelper::rateMatrixHelper();
        std::vector<std::string> labels = helper.getLabels();
        
        Parameter* pExchange = new ParameterExchangabilityRates(rv, this, "exchangability", numStates, labels);
        pExchange->setProposalProbability(5.0);
        parameters.push_back(pExchange);
        
        // set up equilibrium frequencies
        Parameter* pEquil = new ParameterEquilibirumFrequencies(rv, this, "stationary", numStates);
        pEquil->setProposalProbability(5.0);
        parameters.push_back(pEquil);
        }
    
    // delete the alignment objects, leaving behind only the alignment parameters
    for (int i=0; i<wordAlignments.size(); i++)
        delete wordAlignments[i];
    wordAlignments.clear();
    
    // set proposal probabilities
    double sum = 0.0;
    for (int i=0; i<parameters.size(); i++)
        sum += parameters[i]->getProposalProbability();
    for (int i=0; i<parameters.size(); i++)
        parameters[i]->setProposalProbability( parameters[i]->getProposalProbability() / sum );
        
    // make certain to allocate memory to hold thread log likelihoods
    threadLnL = new double[wordParameterAlignments.size()];
    for (int i=0; i<wordParameterAlignments.size(); i++)
        threadLnL[i] = 0.0;
        
#   if 0
    for (int i=0; i<parameters.size(); i++)
        parameters[i]->print();
#   endif
}

void Model::initializeTransitionProbabilities(int numStates, nlohmann::json& j) {

    UserSettings& settings = UserSettings::userSettings();

    // initialize the eigen system object and calculate the first set of eigens
    std::cout << "   * Initializing likelihood-calculation machinery" << std::endl;
    if (substitutionModel == gtr || substitutionModel == custom)
        {
        EigenSystem& eigs = EigenSystem::eigenSystem();
        eigs.initialize(numStates);
        eigs.updateRateMatrix(getExchangabilityRates(), getEquilibriumFrequencies());
        }
    
    // initialize the transition probabilities
    TransitionProbabilities& tProbs = TransitionProbabilities::transitionProbabilties();
    tProbs.initialize( this, getTree()->getNumNodes(), numStates, settings.getSubstitutionModel() );
    tProbs.setNeedsUpdate(true);
    tProbs.setTransitionProbabilities();
}

double Model::lnLikelihood(void) {

#   if 1
    // set up thread pool
    ThreadPool& workers = ThreadPool::threadPoolInstance();
    Tree* tree = getTree();

    for (int i=0; i<wordParameterAlignments.size(); i++)
        {
        workers.post([=] {
            wordLnLike(i, wordParameterAlignments[i], tree);
            });
        }

    workers.wait();
    workers.reset();
    
    double lnL = 0.0;
    for (int i=0; i<wordParameterAlignments.size(); i++)
        lnL += threadLnL[i];
        
    return lnL;
    
#   else

    // calculate likelihood under TKF91 model
    Tree* tree = getTree();
    double lnL = 0.0;
    for (int i=0; i<wordParameterAlignments.size(); i++)
        {
        LikelihoodTkf likelihood(wordParameterAlignments[i], tree, this, substitutionModel);
        double lnP = likelihood.tkfLike();
        lnL += lnP;
        }
    return lnL;
    
#   endif
}

double Model::lnPriorProbability(void) {

    double lnP = 0.0;
    for (int i=0; i<parameters.size(); i++)
        lnP += parameters[i]->lnPriorProbability();
    return lnP;
}

nlohmann::json Model::parseJsonFile(void) {

    // parse the file
    UserSettings& settings = UserSettings::userSettings();
    std::string fn = settings.getDataFile();
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
        
    // print the json information to a file
    std::string configPath = settings.getOutFile();
    configPath += ".config";
    std::ofstream configStrm;
    configStrm.open( configPath.c_str(), std::ios::out );
    if (!configStrm)
        Msg::error("Cannot open file \"" + configPath + "\"");
    configStrm << j.dump() << std::endl;
    configStrm.close();
    
    // return the json object
    return j;
}

void Model::reject(void) {

    Parameter* parm = parameters[updatedParameterIdx];
    parm->reject();
    if (parm->getUpdateChangesEigens() == true)
        {
        // flip eigen index to original state
        EigenSystem& eigs = EigenSystem::eigenSystem();
        eigs.flipActiveValues();
        }
    if (parm->getUpdateChangesTransitionProbabilities() == true)
        {
        TransitionProbabilities& tip = TransitionProbabilities::transitionProbabilties();
        tip.flipActive();
        }
}

double Model::update(void) {

    double u = rv->uniformRv();
    double sum = 0.0;
    for (int i=0; i<parameters.size(); i++)
        {
        sum += parameters[i]->getProposalProbability();
        if (u < sum)
            {
            updatedParameterIdx = i;
            return parameters[i]->update();
            }
        }
    Msg::error("Failed to pick a parameter to update");
    return 0.0;
}

void Model::wordLnLike(int i, ParameterAlignment* aln, Tree* t) {

    LikelihoodTkf likelihood(aln, t, this);
    threadLnL[i] = likelihood.tkfLike();
}
