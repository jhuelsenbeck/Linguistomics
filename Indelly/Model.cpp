#include <future>
#include <fstream>
#include <iostream>
#include <list>
#include <mutex>
#include "Alignment.hpp"
#include "EigenSystem.hpp"
#include "LikelihoodCalculator.hpp"
#include "LikelihoodTkf.hpp"
#include "Model.hpp"
#include "Msg.hpp"
#include "ParameterAlignment.hpp"
#include "ParameterEquilibirumFrequencies.hpp"
#include "ParameterExchangabilityRates.hpp"
#include "ParameterIndelRates.hpp"
#include "ParameterIndelGammaShape.hpp"
#include "ParameterTree.hpp"
#include "Partition.hpp"
#include "RandomVariable.hpp"
#include "RateMatrix.hpp"
#include "RateMatrixHelper.hpp"
#include "SiteLikelihood.hpp"
#include "threads.hpp"
#include "TransitionProbabilities.hpp"
#include "Tree.hpp"
#include "UserSettings.hpp"

class WordLnLikeTask : public ThreadTask {
    public:
        WordLnLikeTask(LikelihoodCalculator* calculator, double* threadLnL, double* wordLnL) {
            Calculator = calculator;
            ThreadLnL  = threadLnL;
            WordLnL    = wordLnL;
        }

        virtual void Run() {
            double lnL = Calculator->lnLikelihood();
            *ThreadLnL = lnL;
            *WordLnL   = lnL;
        }

    private:
        LikelihoodCalculator* Calculator;
        double*               ThreadLnL;
        double*               WordLnL;
};

Model::Model(RandomVariable* r, ThreadPool* p) {

    std::cout << "   Model" << std::endl;
    rv = r;
    threadPool = p;
    partitionInfo = NULL;

    // read the json file
    nlohmann::json j = parseJsonFile();

    // initialize alignments
    std::vector<Alignment*> wordAlignments = initializeAlignments(j);
    
    // initialize parameters of model
    initializeParameters(wordAlignments, j);
    
    // initialize transition probabilities
    initializeTransitionProbabilities(wordAlignments);

    // delete the alignments
    for (int i=0; i<wordAlignments.size(); i++)
        delete wordAlignments[i];

    std::cout << std::endl;
}

Model::~Model(void) {

    delete [] threadLnL;
    for (int i=0; i<parameters.size(); i++)
        delete parameters[i];
    for (int i=0; i<wordLikelihoodCalculators.size(); i++)
        delete wordLikelihoodCalculators[i];
    if (partitionInfo != NULL)
        delete partitionInfo;
}

void Model::accept(void) {

    parameters[updatedParameterIdx]->accept();
}

void Model::flipActiveLikelihood(void) {

    for (int i=0; i<activeLikelihood.size(); i++)
        {
        if (activeLikelihood[i] == 0)
            activeLikelihood[i] = 1;
        else
            activeLikelihood[i] = 0;
        }
}

void Model::flipActiveLikelihood(int idx) {

    if (activeLikelihood[idx] == 0)
        activeLikelihood[idx] = 1;
    else
        activeLikelihood[idx] = 0;
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

std::vector<double>& Model::getIndelGammaRates(void) {

    ParameterIndelGammaShape* p = NULL;
    for (int i=0; i<parameters.size(); i++)
        {
        p = dynamic_cast<ParameterIndelGammaShape*>(parameters[i]);
        if (p != NULL)
            break;
        }
    return p->getRates();
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

ParameterTree* Model::getParameterTree(void) {

    for (int i=0; i<parameters.size(); i++)
        {
        ParameterTree* pT = dynamic_cast<ParameterTree*>(parameters[i]);
        if (pT != NULL)
            return pT;
        }
    return NULL;
}

std::string Model::getStateSetsJsonString(void) {

    if (stateSets.size() == 0)
        return "";

    std::string str = "\"PartitionSets\": [";
    
    for (std::map<std::string,std::set<int> >::iterator it = stateSets.begin(); it != stateSets.end(); it++)
        {
        if (it != stateSets.begin())
            str += ",";
        str += "{\"Name\": \"" + it->first + "\", \"Set\": [";
        std::set<int>& ss = it->second;
        for (std::set<int>::iterator sit = ss.begin(); sit != ss.end(); sit++)
            {
            if (sit != ss.begin())
                str += ",";
            str += std::to_string(*sit);
            }
        str += "]}";
        }
    
    str += "]";
    
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

Tree* Model::getTree(RbBitSet& mask) {

    Tree* t = NULL;
    for (int i=0; i<parameters.size(); i++)
        {
        if (parameters[i]->getName() == "tree")
            {
            ParameterTree* pt = dynamic_cast<ParameterTree*>(parameters[i]);
            return pt->getActiveTree(mask);
            }
        }
    return t;
}

Tree* Model::getTree(const RbBitSet& mask) {

    Tree* t = NULL;
    for (int i=0; i<parameters.size(); i++)
        {
        if (parameters[i]->getName() == "tree")
            {
            ParameterTree* pt = dynamic_cast<ParameterTree*>(parameters[i]);
            return pt->getActiveTree(mask);
            }
        }
    return t;
}

std::string Model::getUpdatedParameterName(void) {

    return parameters[updatedParameterIdx]->getName();
}

std::vector<Alignment*> Model::initializeAlignments(nlohmann::json& j) {

    UserSettings& settings = UserSettings::userSettings();
    
    // get the list of taxa (the so-called canonical list of taxa)
    auto it = j.find("Taxa");
    if (it == j.end())
        Msg::error("Could not find list of taxa in the JSON file");
    std::vector<std::string> tempList = j["Taxa"];
    canonicalTaxonList = tempList;
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
    for (size_t i=0; i<numWords; i++)
        {
        // get the json object for the word
        nlohmann::json jw = j["Words"][i];

        // instantiate the word alignment from the json object
        Alignment* aln = new Alignment(jw, numStates, canonicalTaxonList);
                
        // check if we have problems, deleting the alignment if we do
        if (settings.getUseOnlyCompleteWords() == true)
            {
            // we only use the completely-sampled words
            std::vector<std::string> alnTaxonNames = aln->getTaxonNames();
            if (canonicalTaxonList.size() != aln->numCompleteTaxa())
                {
                rejectedWords.push_back(aln->getName());
                delete aln;
                continue;
                }
            else if (aln->hasAllGapColumn() == true)
                {
                Msg::error("Alignment " + aln->getName() + " has a column with all gaps");
                rejectedWords.push_back(aln->getName());
                delete aln;
                continue;
                }
            }
        else
            {
            // we use all of the cognates, except those with fewer than two languages sampled
            if (aln->numCompleteTaxa() < 3)
                {
                rejectedWords.push_back(aln->getName());
                delete aln;
                continue;
                }
            }
                                    
        words.push_back(aln);
        }

    if (words.size() < 1)
        Msg::error("No word alignments were read!");
        
    std::cout << "   * Number of words                       = " << words.size() << std::endl;
    std::cout << "   * Number of taxa in each word alignment = " << canonicalTaxonList.size() << std::endl;
    if (rejectedWords.size() > 0)
        {
        std::cout << "   * These words were not included because at least one " << std::endl;
        std::cout << "     taxon had no word segments: ";
        size_t cnt = 28;
        for (int i=0; i<rejectedWords.size(); i++)
            {
            std::cout << rejectedWords[i];
            if (i+(size_t)1 != rejectedWords.size())
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
    
    // if the model is a custom one, make certain there is a partition of states to read
    if (substitutionModel == custom)
        {
        it = j.find("PartitionSets");
        if (it == j.end())
            Msg::error("Could not finda  partition of the substitution model states");
        nlohmann::json jsonPart = j["PartitionSets"];
        
        partitionInfo = new Partition(jsonPart);
        //partitionInfo->print();
            
        RateMatrixHelper& helper = RateMatrixHelper::rateMatrixHelper();
        helper.initialize(numStates, partitionInfo);
        helper.print();
        //helper.printMap();
        
        initializeStateSets(jsonPart);
        }
    
    // set up the tree parameter
    Parameter* pTree = new ParameterTree(rv, this, treeStr, canonicalTaxonList, wordAlignments, settings.getInverseTreeLength());
    pTree->setProposalProbability(10.0);
    parameters.push_back(pTree);
    int numNodes = ((ParameterTree*)pTree)->getActiveTree()->getNumNodes();

    // set up the indel parameter
    Parameter* pIndel = new ParameterIndelRates(rv, this, "indel", 7.0, 100.0, 100.0);
    pIndel->setProposalProbability(1.0);
    parameters.push_back(pIndel);
    
    // set up the indel rate variation parameter
    if (settings.getNumIndelCategories() > 1)
        {
        Parameter* pIndelGamma = new ParameterIndelGammaShape(rv, this, "indel gamma", 2.0, settings.getNumIndelCategories());
        pIndelGamma->setProposalProbability(1.0);
        parameters.push_back(pIndelGamma);
        }

    // set up the alignment parameter(s)
    for (int i=0; i<wordAlignments.size(); i++)
        {
        std::string alnName = wordAlignments[i]->getName();
        Parameter* pAlign = new ParameterAlignment(rv, this, wordAlignments[i], alnName, new SiteLikelihood(numNodes,numStates), i);
        pAlign->setProposalProbability(1.0);
        parameters.push_back(pAlign);
        wordParameterAlignments.push_back( dynamic_cast<ParameterAlignment*>(pAlign) );
        wordLikelihoodCalculators.push_back( new LikelihoodCalculator(dynamic_cast<ParameterAlignment*>(pAlign), this) );
        }
        
    // check for consistency between alignment(s) and tree
    std::vector<std::string> treeTaxonNames = getTree()->getTaxonNames();
    if (treeTaxonNames.size() != canonicalTaxonList.size())
        Msg::error("Mismatch between the size of the starting tree and the alignments");
    for (int i=0; i<treeTaxonNames.size(); i++)
        {
        if (treeTaxonNames[i] != canonicalTaxonList[i])
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
        
    // set up vectors for the likelihood calculations
    updateLikelihood.resize(wordParameterAlignments.size());
    activeLikelihood.resize(wordParameterAlignments.size());
    wordLnLikelihoods[0].resize(wordParameterAlignments.size());
    wordLnLikelihoods[1].resize(wordParameterAlignments.size());
    for (int i=0; i<wordParameterAlignments.size(); i++)
        {
        updateLikelihood[i] = false;
        activeLikelihood[i] = 0;
        wordLnLikelihoods[0][i] = 0.0;
        wordLnLikelihoods[1][i] = 0.0;
        }

#   if 0
    for (int i=0; i<parameters.size(); i++)
        parameters[i]->print();
#   endif
}

void Model::initializeStateSets(nlohmann::json& json) {

    // get the state groupings from the json object
    int numGroups = (int)json.size();
    for (int i=0; i<numGroups; i++)
        {
        std::string groupName = json[i]["Name"];
        std::vector<int> group = json[i]["Set"];
        std::set<int> groupSet;
        for (int j=0; j<group.size(); j++)
            groupSet.insert(group[j]);
        stateSets.insert( std::make_pair(groupName,groupSet) );
        }
}

void Model::initializeTransitionProbabilities(std::vector<Alignment*>& wordAlignments) {

    std::cout << "   * Initializing likelihood-calculation machinery" << std::endl;
    
    int numStates = wordAlignments[0]->getNumStates();
    UserSettings& settings = UserSettings::userSettings();
    
    // set up the rate matrix
    RateMatrix& rmat = RateMatrix::rateMatrix();
    rmat.initialize(numStates, settings.getUseEigenSystem());
    if (settings.getUseEigenSystem() == true && settings.getSubstitutionModel() != jc69)
        {
        EigenSystem& eigs = EigenSystem::eigenSystem();
        eigs.initialize(numStates);
        rmat.updateRateMatrix(getExchangabilityRates(), getEquilibriumFrequencies());
        }
    else if (settings.getUseEigenSystem() == false && settings.getSubstitutionModel() != jc69)
        rmat.updateRateMatrix(getExchangabilityRates(), getEquilibriumFrequencies());
    
    // initialize the transition probabilities
    TransitionProbabilities& tProbs = TransitionProbabilities::transitionProbabilties();
    tProbs.initialize( this, threadPool, wordAlignments, getTree()->getNumNodes(), numStates, settings.getSubstitutionModel() );
    tProbs.setNeedsUpdate(true);
    tProbs.setTransitionProbabilities();
}

double Model::lnLikelihood(void) {

#   if 1

    // set up thread pool for calculating the likelihood under the TKF91 model
    for (int i=0; i<wordParameterAlignments.size(); i++)
        {
        if (updateLikelihood[i] == true)
            {
            threadPool->PushTask(new WordLnLikeTask(wordLikelihoodCalculators[i], &threadLnL[i], &wordLnLikelihoods[ activeLikelihood[i] ][i]));
            }
        else
            {
            threadLnL[i] = wordLnLikelihoods[ activeLikelihood[i] ][i];
            }
        updateLikelihood[i] = false;
        }
    
    threadPool->Wait();
    
#   else

    // calculate likelihood under TKF91 model
    for (int i=0; i<wordParameterAlignments.size(); i++)
        {
        //double lnP = wordLikelihoodCalculators[i]->lnLikelihood();
        if (updateLikelihood[i] == true)
            {
            threadLnL[i] = wordLikelihoodCalculators[i]->lnLikelihood();
            wordLnLikelihoods[ activeLikelihood[i] ][i] = threadLnL[i];
            }
        else
            {
            threadLnL[i] = wordLnLikelihoods[ activeLikelihood[i] ][i];
            }
        updateLikelihood[i] = false;
        }
                
#   endif

    double lnL = 0.0;
    for (int i=0; i<wordParameterAlignments.size(); i++)
        lnL += threadLnL[i];
    return lnL;
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
    if (parm->getUpdateChangesRateMatrix() == true)
        {
        // flip rate matrix index to original state
        RateMatrix& rmat = RateMatrix::rateMatrix();
        rmat.flipActiveValues();
        }
    if (parm->getUpdateChangesTransitionProbabilities() == true)
        {
        TransitionProbabilities& tip = TransitionProbabilities::transitionProbabilties();
        tip.flipActive();
        }
}

void Model::setUpdateLikelihood(void) {

    for (int i=0; i<updateLikelihood.size(); i++)
        updateLikelihood[i] = true;
}

void Model::setUpdateLikelihood(int idx) {

    updateLikelihood[idx] = true;
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

void Model::wordLnLike(int i) {

    threadLnL[i] = wordLikelihoodCalculators[i]->lnLikelihood();
    wordLnLikelihoods[ activeLikelihood[i] ][i] = threadLnL[i];
}

