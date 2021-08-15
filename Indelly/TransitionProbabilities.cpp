#include <cmath>
#include <iomanip>
#include <iostream>
#include "EigenSystem.hpp"
#include "Model.hpp"
#include "Msg.hpp"
#include "Node.hpp"
#include "RateMatrix.hpp"
#include "TransitionProbabilities.hpp"
#include "Tree.hpp"
#include "UserSettings.hpp"



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

void TransitionProbabilities::initialize(Model* m, int nn, int ns, int sm) {

    if (isInitialized == true)
        {
        Msg::warning("Transition probabilities can only be initialized once");
        return;
        }
        
    UserSettings& settings = UserSettings::userSettings();
    modelPtr = m;
    numNodes = nn;
    numStates = ns;
    substitutionModel = sm;
    numRateCategories = settings.getNumRateCategories();
    int numStatesSquared = numStates * numStates;
    std::cout << "   * Number of states = " << numStates << std::endl;
    std::cout << "   * Number of gamma rate categories = " << numRateCategories << std::endl;
        
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

void TransitionProbabilities::setTransitionProbabilitiesUsingEigenSystem(void) {

    EigenSystem& eigs = EigenSystem::eigenSystem();
    Eigen::Matrix<std::complex<double>, Eigen::Dynamic, 1>& ceigenvalue = eigs.getEigenValues();
    std::complex<double>* ccIjk = eigs.getCijk();
    std::vector<std::complex<double> > ceigValExp(numStates);

    std::vector<Node*>& traversalSeq = modelPtr->getTree()->getDownPassSequence();
    for (int n=0; n<traversalSeq.size(); n++)
        {
        Node* p = traversalSeq[n];
        double** tp = probs[activeProbs][p->getIndex()];
        
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
                tp[i][j] = (sum.real() < 0.0) ? 0.0 : sum.real();
                }
            }
        }
        
    RateMatrix& rmat = RateMatrix::rateMatrix();
    stationaryFreqs[activeProbs] = rmat.getEquilibriumFrequencies();
}

void TransitionProbabilities::setTransitionProbabilitiesUsingPadeMethod(void) {

    Msg::error("Pade method is not yet implemented");

#   if 0
	int dim = A.dim1();
	if (dim != A.dim2())
		return (1);
	
	// create identity matrices
	MbMatrix<double> D(dim,dim,0.0);
	for (int i=0; i<dim; i++)
		D[i][i] = 1.0;
	MbMatrix<double> N(D.copy()), X(D.copy());

	// create uninitialized matrix
	MbMatrix<double> cX(dim, dim);
	
	// We assume that we have a rate matrix where rows sum to zero
	// Then the infinity-norm is twice the maximum absolute value
	// of the diagonal cells.
	double normA = 0.0;
	for (int i=0; i<dim; i++) {
		double x = fabs (A[i][i]);
		if (x > normA)
			normA = x;
	}
	normA *= 2.0;

	// Calculate 1 + floor (log2(normA))
	int y;
	frexp(normA, &y);	// this will give us the floor(log2(normA)) part in y
	y++;

	// Get max(0,y)
	int j = 0;
	if (y > 0)
		j = y;

	// divide A by scalar 2^j
	A /= ldexp (1.0, j);
	
	double c = 1.0;
	for (int k=1; k<=qValue; k++)
        {
		c = c * (qValue - k + 1.0) / ((2.0 * qValue - k + 1.0) * k);

		/* X = AX */
		X = A * X;

		/* N = N + cX */
		cX = c * X;
		N = N + cX;

		/* D = D + (-1)^k*cX */
		if (k % 2 == 0)
			D = D + cX;
		else
			D = D - cX;
		}

	MbMath::gaussianElimination(D, N, F);

	for (int k=0; k<j; k++)
		F = F * F;
	
	for (int i=0; i<dim; i++)
		{
		for (j=0; j<dim; j++)
			{
			if (F[i][j] < 0.0)
				F[i][j] = fabs(F[i][j]);
			}
		}
	return (0);
#   endif
}
