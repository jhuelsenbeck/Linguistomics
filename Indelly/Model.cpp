#include <future>
#include <iostream>
#include <list>
#include <mutex>
#include "Alignment.hpp"
#include "EigenSystem.hpp"
#include "IoManager.hpp"
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

    UserSettings& settings = UserSettings::userSettings();

    rv = r;
    
    // read the data file
    std::vector<Alignment*> wordAlignments;
    std::vector<std::string> taxonNames = readWords(settings.getDataFile(), wordAlignments);
    threadLnL = new double[wordAlignments.size()];
    int numStates = wordAlignments[0]->getNumStates();
    if (numStates > wordAlignments[0]->getMaximumNumberOfStates())
        Msg::error("The maximum number of states is " + std::to_string(wordAlignments[0]->getMaximumNumberOfStates()));
        
    std::cout << "   Model" << std::endl;
    // set up the tree parameter
    Parameter* pTree = new ParameterTree(rv, this, settings.getTreeFile(), taxonNames, settings.getInverseTreeLength());
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
    if (taxonNames.size() != getTree()->getTaxonNames().size())
        Msg::error("Mismatch in the number of taxa read in from the files and in the tree");
    for (int i=0; i<taxonNames.size(); i++)
        {
        if (getTree()->isTaxonPresent(taxonNames[i]) == false)
            Msg::error("Mismatch in the taxon composition between the input tree and alignment files");
        }
    
    // set up exchangability parameters
    Parameter* pExchange = new ParameterExchangabilityRates(rv, this, "exchangability", numStates, wordAlignments[0]->getStates());
    pExchange->setProposalProbability(50.0);
    parameters.push_back(pExchange);
    
    // set up equilibrium frequencies
    Parameter* pEquil = new ParameterEquilibirumFrequencies(rv, this, "stationary", numStates, wordAlignments[0]->getStates());
    pEquil->setProposalProbability(50.0);
    parameters.push_back(pEquil);
    
    // delete the alignment objects, leaving behind only the alignment parameters
    for (int i=0; i<wordAlignments.size(); i++)
        delete wordAlignments[i];
        
    std::cout << std::endl;
    
    // initialize the eigen system object and calculate the first set of eigens
    EigenSystem& eigs = EigenSystem::eigenSystem();
    eigs.initialize(numStates);
    eigs.updateRateMatrix(getExchangabilityRates(), getEquilibriumFrequencies());
    
    // initialize the transition probabilities
    TransitionProbabilities& tProbs = TransitionProbabilities::transitionProbabilties();
    tProbs.initialize( this, getTree()->getNumNodes(), numStates );
    tProbs.setNeedsUpdate(true);
    tProbs.setTransitionProbabilities();

    // set proposal probabilities
    double sum = 0.0;
    for (int i=0; i<parameters.size(); i++)
        sum += parameters[i]->getProposalProbability();
    for (int i=0; i<parameters.size(); i++)
        parameters[i]->setProposalProbability( parameters[i]->getProposalProbability() / sum );
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
        LikelihoodTkf likelihood(wordParameterAlignments[i], tree, this);
        lnL += likelihood.tkfLike();
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

std::vector<std::string> Model::readWords(std::string fp, std::vector<Alignment*>& wordAlignments) {

    std::cout << "   Data" << std::endl;

    // file manager
    IoManager fileMngr;
    fileMngr.setFilePath(fp);
    
    std::cout << "   * Reading word alignments in directory \"" << fp << "\"" << std::endl;
    
	// check first by reading all of the alignments, calculating the number of nucleotides for each
	int numAlignmentsToRead = 0;
    std::vector<std::string> taxonNames;
	for (int fn=0; fn<fileMngr.getNumFilesInDirectory(); fn++)
		{
		std::string fileName = fileMngr.getFileNumber(fn);
		Alignment* alignmentPtr = new Alignment(fp+ "/" + fileName, true);
		if (fn == 0)
			{
            taxonNames = alignmentPtr->getTaxonNames();
			}
		else
			{
			std::vector<std::string> namesFromAln = alignmentPtr->getTaxonNames();
			bool isSame = true;
			if ( namesFromAln.size() != taxonNames.size() )
				isSame = false;
			else
				{
				for (int i=0; i<taxonNames.size(); i++)
					if ( namesFromAln[i] != taxonNames[i] )
						isSame = false;
				}
			if (isSame == false)
				{
                for (int i=0; i<taxonNames.size(); i++)
                    std::cout << i << " -- " << taxonNames[i] << " " << namesFromAln[i] << std::endl;;
                Msg::error("Missmatch in the taxon names of the alignments: " + fileName);
				}
			}
		numAlignmentsToRead++;

        std::string lastPathComponent = fileName;
        const size_t lastSlashIdx = lastPathComponent.find_last_of("/");
        if (std::string::npos != lastSlashIdx)
            lastPathComponent.erase(0, lastSlashIdx + 1);
        const size_t periodIdx = lastPathComponent.rfind('.');
        if (std::string::npos != periodIdx)
            lastPathComponent.erase(periodIdx);

        alignmentPtr->setName(lastPathComponent);
        wordAlignments.push_back(alignmentPtr);
		}
		
    std::cout << "   * Number of words                       = " << numAlignmentsToRead << std::endl;
    std::cout << "   * Number of taxa in each word alignment = " << taxonNames.size() << std::endl;
    std::cout << std::endl;
			
#   if 0
    for (int i=0; i<wordAlignments.size(); i++)
        wordAlignments[i]->print();
#   endif

    return taxonNames;
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
