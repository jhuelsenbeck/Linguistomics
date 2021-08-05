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
#include "ThreadPool.hpp"
#include "TransitionProbabilities.hpp"
#include "Tree.hpp"
#include "UserSettings.hpp"




Model::Model(RandomVariable* r) {

    std::cout << "   Model" << std::endl;

    rv = r;
    UserSettings& settings = UserSettings::userSettings();

    // read the json file
    nlohmann::json j = parseJsonFile(settings.getDataFile());

    // initialize alignments
    std::vector<Alignment*> wordAlignments = initializeAlignments(j);
    std::vector<std::string> taxonNames = wordAlignments[0]->getTaxonNames();
    int numStates = wordAlignments[0]->getNumStates();
    
    // initialize parameters of model
    initializeParameters(&settings, wordAlignments, j);
    
    // initialize the eigen system object and calculate the first set of eigens
    if (substitutionModel == "GTR")
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

    // check that there are words in the json object
    auto it = j.find("Words");
    if (it == j.end())
        Msg::error("Could not find word information in the JSON file");
    size_t numWords = j["Words"].size();
    std::cout << "   * Found " << numWords << " words in JSON file" << std::endl;
    
    // check that there is a list of valid states
    it = j.find("States");
    if (it == j.end())
        Msg::error("Could not find state information in the JSON file");
    std::string validStates = j["States"];
    if (validStates.length() <= 1)
        Msg::error("There must be at least two states in the model");
    
    // read each word
    std::vector<Alignment*> words;
    std::vector<std::string> taxonNames;
    std::string tempStates = "";
    for (size_t i=0; i<numWords; i++)
        {
        // instantiate the word alignment from the json object
        nlohmann::json jw = j["Words"][i];
        Alignment* aln = new Alignment(jw, validStates);
        if (i == 0)
            {
            taxonNames = aln->getTaxonNames();
            tempStates = aln->getStates();
            }
        
        // check if we have problems
        if (taxonNames.size() != aln->numCompleteTaxa())
            {
            Msg::warning("   * Word " + aln->getName() + " will is taxa with no word segments and will not be included");
            delete aln;
            continue;
            }
        std::vector<std::string> alnTaxonNames = aln->getTaxonNames();
        if (taxonNames.size() != alnTaxonNames.size())
            Msg::error("List of taxa is inconsistent across words");
        for (int j=0; j<taxonNames.size(); j++)
            {
            if (taxonNames[j] != alnTaxonNames[j])
                Msg::error("List of taxa is inconsistent across words");
            }
        std::string alnStates = aln->getStates();
        if (tempStates.size() != alnStates.size())
            Msg::error("The set of states is inconsistent acros word alignments");
        for (int j=0; j<tempStates.length(); j++)
            {
            if (tempStates.at(j) != alnStates.at(j))
                Msg::error("The set of states is inconsistent acros word alignments");
            }
            
        words.push_back(aln);
        }

    if (words.size() < 1)
        Msg::error("No word alignments were read!");
        
    std::cout << "   * Number of words                       = " << words.size() << std::endl;
    std::cout << "   * Number of taxa in each word alignment = " << taxonNames.size() << std::endl;
    
#   if 0
    for (int i=0; i<words.size(); i++)
        words[i]->print();
#   endif

    return words;
}

void Model::initializeParameters(UserSettings* settings, std::vector<Alignment*>& wordAlignments, nlohmann::json& j) {

    // attempt to find the initial tree in the json object
    std::string treeStr = "";
    auto it = j.find("Tree");
    if (it != j.end())
        treeStr = j["Tree"];

    // initialize a few important parameters
    substitutionModel = settings->getSubstitutionModel();
    int numStates = wordAlignments[0]->getNumStates();
    std::vector<std::string> taxonNames = wordAlignments[0]->getTaxonNames();
    
    // set up the tree parameter
    Parameter* pTree = new ParameterTree(rv, this, treeStr, taxonNames, settings->getInverseTreeLength());
    pTree->setProposalProbability(10.0);
    parameters.push_back(pTree);
    
    // set up the indel parameter
    Parameter* pIndel = new ParameterIndelRates(rv, this, "indel", 7.0, 100.0, 100.0);
    pIndel->setProposalProbability(1.0);
    parameters.push_back(pIndel);

    // set up the alignment parameter
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
    
    if (substitutionModel == "GTR")
        {
        // set up exchangability parameters
        Parameter* pExchange = new ParameterExchangabilityRates(rv, this, "exchangability", numStates, wordAlignments[0]->getStates());
        pExchange->setProposalProbability(50.0);
        parameters.push_back(pExchange);
        
        // set up equilibrium frequencies
        Parameter* pEquil = new ParameterEquilibirumFrequencies(rv, this, "stationary", numStates, wordAlignments[0]->getStates());
        pEquil->setProposalProbability(50.0);
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
        {
        //std::cout << i << " " << threadLnL[i] << std::endl;
        lnL += threadLnL[i];
        }
        
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

nlohmann::json Model::parseJsonFile(std::string fn) {

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

    LikelihoodTkf likelihood(aln, t, this, substitutionModel);
    threadLnL[i] = likelihood.tkfLike();
}
