#include <cmath>
#include <iomanip>
#include <iostream>
#include "EigenSystem.hpp"
#include "Model.hpp"
#include "Msg.hpp"
#include "Node.hpp"
#include "TransitionProbabilities.hpp"
#include "Tree.hpp"



TransitionProbabilities::TransitionProbabilities(void) {

    isInitialized = false;
    needsUpdate = true;
    activeProbs = 0;
}

TransitionProbabilities::~TransitionProbabilities(void) {

    for (int s=0; s<2; s++)
        {
        for (int n=0; n<probs[s].size(); n++)
            {
            delete [] probs[s][n][0];
            delete [] probs[s][n];
            }
        }
}

void TransitionProbabilities::flipActive(void) {

    if (activeProbs == 0)
        activeProbs = 1;
    else
        activeProbs = 0;
}

void TransitionProbabilities::initialize(Model* m, int nn, int ns, std::string sm) {

    if (isInitialized == true)
        {
        Msg::warning("Transition probabilities can only be initialized once");
        return;
        }
        
    modelPtr = m;
    numNodes = nn;
    numStates = ns;
    substitutionModel = sm;
    int numStatesSquared = numStates * numStates;
    
    
    for (int s=0; s<2; s++)
        {
        probs[s].resize(numNodes);
        for (int n=0; n<probs[s].size(); n++)
            {
            probs[s][n] = new double*[numStates];
            probs[s][n][0] = new double[numStatesSquared];
            for (int i=1; i<numStates; i++)
                probs[s][n][i] = probs[s][n][i-1] + numStates;
            for (int i=0; i<numStates; i++)
                for (int j=0; j<numStates; j++)
                    probs[s][n][i][j] = 0.0;
            }
        }
        
    stationaryFreqs[0].resize(numStates);
    stationaryFreqs[1].resize(numStates);
        
    isInitialized = true;
}

void TransitionProbabilities::print(void) {

    std::cout << std::fixed << std::setprecision(5);
    for (int n=0; n<probs[activeProbs].size(); n++)
        {
        std::cout << "Transition probabilities for node " << n << std::endl;
        for (int i=0; i<numStates; i++)
            {
            for (int j=0; j<numStates; j++)
                std::cout << probs[activeProbs][n][i][j] << " ";
            std::cout << std::endl;
            }
        }
}

void TransitionProbabilities::setTransitionProbabilities(void) {

    if (needsUpdate == false)
        return;
 
    if (substitutionModel == "GTR")
        {
        // calculate transition probabilities for GTR model
        EigenSystem& eigs = EigenSystem::eigenSystem();
        Eigen::Matrix<std::complex<double>, Eigen::Dynamic, 1>& ceigenvalue = eigs.getEigenValues();
        std::complex<double>* ccIjk = eigs.getCijk();

        std::vector<Node*>& traversalSeq = modelPtr->getTree()->getDownPassSequence();
        for (int n=0; n<traversalSeq.size(); n++)
            {
            Node* p = traversalSeq[n];
            double** tp = probs[activeProbs][p->getIndex()];
            double v = p->getBranchLength();
            
            std::complex<double> ceigValExp[numStates];
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
                    tp[i][j] = (sum.real() < 0.0) ? 0.0 : sum.real();
                    }
                }
            }
            
        stationaryFreqs[activeProbs] = eigs.getStationaryFrequencies();
        }
    else
        {
        // calculate transition probabilities under the Jukes-Cantor (1969) model
        std::vector<Node*>& traversalSeq = modelPtr->getTree()->getDownPassSequence();
        for (int n=0; n<traversalSeq.size(); n++)
            {
            Node* p = traversalSeq[n];
            double** tp = probs[activeProbs][p->getIndex()];
            double v = p->getBranchLength();
            
            double x = -((double)numStates/(numStates-1));
            double pChange = (1.0/numStates) - (1.0/numStates) * exp(x * v);
            double pNoChange = (1.0/numStates) + ((double)(numStates-1)/numStates) * exp(x * v);
            for (int i=0; i<numStates; i++)
                {
                for (int j=0; j<numStates; j++)
                    {
                    if (i == j)
                        tp[i][j] = pNoChange;
                    else
                        tp[i][j] = pChange;
                    }
                }
            }
            
        double sf = 1.0 / numStates;
        for (int i=0; i<numStates; i++)
            stationaryFreqs[activeProbs][i] = sf;
        }
        
    needsUpdate = false;
}

