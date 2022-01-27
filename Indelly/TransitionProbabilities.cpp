#include <cmath>
#include <iomanip>
#include <iostream>
#include <unsupported/Eigen/MatrixFunctions>
#include "Alignment.hpp"
#include "EigenSystem.hpp"
#include "Model.hpp"
#include "Msg.hpp"
#include "Node.hpp"
#include "ParameterTree.hpp"
#include "RateMatrix.hpp"
#include "threads.hpp"
#include "TransitionProbabilities.hpp"
#include "Tree.hpp"
#include "UserSettings.hpp"

void padeTransitionProbabilities(Tree* t, const StateMatrix_t& Q, const std::vector<StateMatrix_t*>& probs, const std::vector<StateMatrix_t*>& m);



TransitionProbabilities::TransitionProbabilities(void) {

    isInitialized = false;
    needsUpdate = true;
    activeProbs = 0;
}

TransitionProbabilities::~TransitionProbabilities(void) {

    for (std::map<RbBitSet,TransitionProbabilitiesPair>::iterator it = transProbs.begin(); it != transProbs.end(); it++)
        {
        for (int s=0; s<2; s++)
            {
            for (int n=0; n<it->second.probs[s].size(); n++)
                {
                delete it->second.probs[s][n];
                }
            }
        for (int n=0; n<it->second.aMat.size(); n++)
            {
            delete it->second.aMat[n];
            }
        }
}

void TransitionProbabilities::flipActive(void) {

    if (activeProbs == 0)
        activeProbs = 1;
    else
        activeProbs = 0;
}

std::vector<StateMatrix_t*>& TransitionProbabilities::getTransitionProbabilities(RbBitSet& bs) {

    std::map<RbBitSet,TransitionProbabilitiesPair>::iterator it = transProbs.find(bs);
    if (it == transProbs.end())
        Msg::error("Could not find transition probability vector for mask " + bs.bitString());
    return it->second.probs[activeProbs];
}

StateMatrix_t* TransitionProbabilities::getTransitionProbabilities(RbBitSet& bs, int nodeIdx) {

    std::map<RbBitSet,TransitionProbabilitiesPair>::iterator it = transProbs.find(bs);
    if (it == transProbs.end())
        Msg::error("Could not find transition probability vector for mask " + bs.bitString());
    return it->second.probs[activeProbs][nodeIdx];
}

void TransitionProbabilities::initialize(Model* m, thread_pool* p, std::vector<Alignment*>& alns, int nn, int ns, int sm) {

    if (isInitialized == true)
        {
        Msg::warning("Transition probabilities can only be initialized once");
        return;
        }
                
    UserSettings& settings = UserSettings::userSettings();
    modelPtr = m;
    threadPool = p;
    numNodes = nn;
    numStates = ns;
    substitutionModel = sm;
    numRateCategories = settings.getNumRateCategories();
    std::cout << "   * Number of states = " << numStates << std::endl;
    std::cout << "   * Number of gamma rate categories = " << numRateCategories << std::endl;

    for (int i=0; i<alns.size(); i++)
        {
        // get the taxon mask
        std::vector<bool> m = alns[i]->getTaxonMask();
        RbBitSet mask(m);

        // if the mask is not found in the map, insert it
        std::map<RbBitSet,TransitionProbabilitiesPair>::iterator it = transProbs.find(mask);
        if (it == transProbs.end())
            {
            Tree* t = modelPtr->getTree(mask);
            if (t == NULL)
                Msg::error("Could not find tree for mask " + mask.bitString() + " when initializing transition probabilities");
            TransitionProbabilitiesPair pair;
            for (int s=0; s<2; s++)
                {
                pair.probs[s].resize(t->getNumNodes());
                for (int n=0; n<pair.probs[s].size(); n++)
                    {
                    pair.probs[s][n] = new StateMatrix_t;
                    pair.probs[s][n]->resize(numStates,numStates);
                    }
                }
            pair.aMat.resize(t->getNumNodes());
            for (int n=0; n<pair.aMat.size(); n++)
                {
                pair.aMat[n] = new StateMatrix_t;
                pair.aMat[n]->resize(numStates,numStates);
                }
            
            transProbs.insert( std::make_pair(mask,pair) );
            }
        }
        
    stationaryFreqs[0].resize(numStates);
    stationaryFreqs[1].resize(numStates);
        
    isInitialized = true;
    
    //print();
}

void TransitionProbabilities::print(void) {

    std::cout << std::fixed << std::setprecision(5);
    for (std::map<RbBitSet,TransitionProbabilitiesPair>::iterator it = transProbs.begin(); it != transProbs.end(); it++)
        {
        for (int n=0; n<it->second.probs[activeProbs].size(); n++)
            {
            std::cout << "Transition probabilities for node " << n << " (" << it->first.bitString() << ")" << std::endl;
            for (int i=0; i<numStates; i++)
                {
                for (int j=0; j<numStates; j++)
                    std::cout << (*it->second.probs[activeProbs][n])(i,j) << " ";
                std::cout << std::endl;
                }
            }
        }
}

void TransitionProbabilities::setTransitionProbabilities(void) {

    if (needsUpdate == false)
        return;
        
    if (substitutionModel == jc69)
        {
        // the transition probabilities can be calculated analytically (and quickly)
        setTransitionProbabilitiesJc69();
        }
    else
        {
        // the transition probabilities are calculated using either the Eigen
        // system of the rate matrix or using the Pade approximation
        RateMatrix& rmat = RateMatrix::rateMatrix();
        if (rmat.getUseEigenSystem() == true)
            setTransitionProbabilitiesUsingEigenSystem();
        else
            setTransitionProbabilitiesUsingPadeMethod();
        }
        
    needsUpdate = false;
}

void TransitionProbabilities::setTransitionProbabilitiesJc69(void) {

    // calculate transition probabilities under the Jukes-Cantor (1969) model
    for (std::map<RbBitSet,TransitionProbabilitiesPair>::iterator it = transProbs.begin(); it != transProbs.end(); it++)
        {
        Tree* t = modelPtr->getTree(it->first);
        if (t == NULL)
            Msg::error("Could not find tree for mask " + it->first.bitString());
                
        std::vector<StateMatrix_t*>& probs = it->second.probs[activeProbs];
        
        std::vector<Node*>& traversalSeq = t->getDownPassSequence();
        for (int n=0; n<traversalSeq.size(); n++)
            {
            Node* p = traversalSeq[n];
            StateMatrix_t* tp = probs[p->getIndex()];
            double v = p->getBranchLength();
            
            double x = -((double)numStates/(numStates-1));
            double pChange = (1.0/numStates) - (1.0/numStates) * exp(x * v);
            double pNoChange = (1.0/numStates) + ((double)(numStates-1)/numStates) * exp(x * v);
            for (int i=0; i<numStates; i++)
                {
                for (int j=0; j<numStates; j++)
                    {
                    if (i == j)
                        (*tp)(i,j) = pNoChange;
                    else
                        (*tp)(i,j) = pChange;
                    }
                }
            }
        
        }
            
    double sf = 1.0 / numStates;
    for (int i=0; i<numStates; i++)
        stationaryFreqs[activeProbs][i] = sf;
}

void TransitionProbabilities::setTransitionProbabilitiesUsingEigenSystem(void) {

    EigenSystem& eigs = EigenSystem::eigenSystem();
    Eigen::Matrix<std::complex<double>, Eigen::Dynamic, 1>& ceigenvalue = eigs.getEigenValues();
    std::complex<double>* ccIjk = eigs.getCijk();
    std::vector<std::complex<double> > ceigValExp(numStates);

    for (std::map<RbBitSet,TransitionProbabilitiesPair>::iterator it = transProbs.begin(); it != transProbs.end(); it++)
        {
        Tree* t = modelPtr->getTree(it->first);
        if (t == NULL)
            Msg::error("Could not find tree for mask " + it->first.bitString());
                
        std::vector<StateMatrix_t*>& probs = it->second.probs[activeProbs];
        
        std::vector<Node*>& traversalSeq = modelPtr->getTree(it->first)->getDownPassSequence();
        for (int n=0; n<traversalSeq.size(); n++)
            {
            Node* p = traversalSeq[n];
            StateMatrix_t* tp = probs[p->getIndex()];
            
            double v = p->getBranchLength();
            for (int s=0; s<numStates; s++)
                ceigValExp[s] = exp(ceigenvalue[s] * v);

            std::complex<double>* ptr = ccIjk;
            for (int i=0; i<numStates; i++)
                {
                for (int j=0; j<numStates; j++)
                    {
                    std::complex<double> sum = std::complex<double>(0.0, 0.0);
                    for(int s=0; s<numStates; s++)
                        sum += (*ptr++) * ceigValExp[s];
                    (*tp)(i,j) = (sum.real() < 0.0) ? 0.0 : sum.real();
                    }
                }
            }
            
        }
        
    RateMatrix& rmat = RateMatrix::rateMatrix();
    stationaryFreqs[activeProbs] = rmat.getEquilibriumFrequencies();
}

void TransitionProbabilities::setTransitionProbabilitiesUsingPadeMethod(void) {

#   if 1

    // threaded version
    RateMatrix& rmat = RateMatrix::rateMatrix();
    const StateMatrix_t& Q = rmat.getRateMatrix();
    
    // update the main tree
    for (std::map<RbBitSet,TransitionProbabilitiesPair>::iterator it = transProbs.begin(); it != transProbs.end(); it++)
        {
        Tree* t = modelPtr->getTree(it->first);
        if (t == NULL)
            Msg::error("Could not find tree for mask " + it->first.bitString());
        const std::vector<StateMatrix_t*>& probs = it->second.probs[activeProbs];
        const std::vector<StateMatrix_t*>& m     = it->second.aMat;
        
        threadPool->push_task(padeTransitionProbabilities, t, Q, probs, m);
        }
        
    threadPool->wait_for_tasks();
   
#   else

    // serial version
    RateMatrix& rmat = RateMatrix::rateMatrix();
    Eigen::MatrixXd M(numStates,numStates);
    StateMatrix_t& Q = rmat.getRateMatrix();
    
    // update the main tree
    for (std::map<RbBitSet,TransitionProbabilitiesPair>::iterator it = transProbs.begin(); it != transProbs.end(); it++)
        {
        Tree* t = modelPtr->getTree(it->first);
        if (t == NULL)
            Msg::error("Could not find tree for mask " + it->first.bitString());
        std::vector<StateMatrix_t*>& probs = it->second.probs[activeProbs];
        
        std::vector<Node*>& traversalSeq = t->getDownPassSequence();
        std::vector<double> branchLengths(traversalSeq.size());
        for (int n=0; n<traversalSeq.size(); n++)
            {
            Node* p = traversalSeq[n];
            StateMatrix_t* tp = probs[p->getIndex()];
            double v = p->getBranchLength();
            branchLengths[p->getIndex()] = v;
            
            M = Q * v;
            (*tp) = M.exp();
            }
        }

#   endif
        
    std::vector<double>& rmatFreqs = rmat.getEquilibriumFrequencies();
    for (int i=0; i<numStates; i++)
        stationaryFreqs[activeProbs][i] = rmatFreqs[i];
}

void padeTransitionProbabilities(Tree* t, const StateMatrix_t& Q, const std::vector<StateMatrix_t*>& probs, const std::vector<StateMatrix_t*>& m) {

    std::vector<Node*>& traversalSeq = t->getDownPassSequence();
    for (int n=0; n<traversalSeq.size(); n++)
        {
        Node* p = traversalSeq[n];
        int idx = p->getIndex();
        StateMatrix_t* tp = probs[idx];
        StateMatrix_t* M = m[idx];
        double v = p->getBranchLength();
        (*M) = Q * v;
        (*tp) = M->exp();
        }
}
