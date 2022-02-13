#include <cmath>
#include <iomanip>
#include <iostream>
#include "Alignment.hpp"
#include "MatrixMath.hpp"
#include "Model.hpp"
#include "Msg.hpp"
#include "Node.hpp"
#include "ParameterTree.hpp"
#include "RateMatrix.hpp"
#include "Threads.hpp"
#include "TransitionProbabilities.hpp"
#include "Tree.hpp"
#include "UserSettings.hpp"

void computeMatrixExponential(DoubleMatrix* Q, int qValue, double v, DoubleMatrix* A, DoubleMatrix* P, DoubleMatrix* D, DoubleMatrix* N, DoubleMatrix* X, DoubleMatrix* cX, int numStates, DoubleMatrix* scratch1, DoubleMatrix* scratch2, double* scratchVec);
double factorial(int x);
int logBase2Plus1(double x);
int setQvalue(double tol);



class TransitionProbabilitiesTask: public ThreadTask {

    public:
        TransitionProbabilitiesTask(void) {
            numStates  = 0;
            Tree       = NULL;
            Q          = NULL;
            A          = NULL;
            P          = NULL;
            D          = NULL;
            N          = NULL;
            X          = NULL;
            cX         = NULL;
            scratch1   = NULL;
            scratch2   = NULL;
            scratchVec = NULL;
        }
        
        void init(Tree* tree, DoubleMatrix* q, TransitionProbabilitiesInfo& info, int activeProbs) {
        
            numStates  = info.numStates;
            Tree       = tree;
            Q          = q;
            A          = info.a_mat;
            P          = info.probs[activeProbs];
            D          = info.d_mat;
            N          = info.n_mat;
            X          = info.x_mat;
            cX         = info.cX_mat;
            scratch1   = info.scratch_mat1;
            scratch2   = info.scratch_mat2;
            scratchVec = info.scratch_vec;
        }

        virtual void Run(void) {
        
            int qValue = setQvalue(10e-7);
            std::vector<Node*>& traversalSeq = Tree->getDownPassSequence();
            for (int n = 0; n < traversalSeq.size(); n++)
                {
                Node* p = traversalSeq[n];
                if (p->getAncestor() != NULL)
                    {
                    int pIdx = p->getIndex();
                    double v = p->getBranchLength();
                    computeMatrixExponential(Q, qValue, v, A, P[pIdx], D, N, X, cX, numStates, scratch1, scratch2, scratchVec);
                    }
                }
            }

    private:
        int             numStates;
        Tree*           Tree;
        DoubleMatrix*   Q;
        DoubleMatrix*   A;
        DoubleMatrix**  P;
        DoubleMatrix*   D;
        DoubleMatrix*   N;
        DoubleMatrix*   X;
        DoubleMatrix*   cX;
        DoubleMatrix*   scratch1;
        DoubleMatrix*   scratch2;
        double*         scratchVec;
};



TransitionProbabilities::TransitionProbabilities(void) {

    modelPtr          = NULL;
    numNodes          = 0;
    numRateCategories = 0;
    numStates         = 0;
    substitutionModel = NULL;
    threadPool        = NULL;
    isInitialized     = false;
    needsUpdate       = true;
    activeProbs       = 0;
}

TransitionProbabilities::~TransitionProbabilities(void) {

    for (std::map<RbBitSet,TransitionProbabilitiesInfo>::iterator it = transProbs.begin(); it != transProbs.end(); it++)
        {
        for (int s=0; s<2; s++)
            delete [] it->second.probs[s];
        delete it->second.a_mat;
        delete it->second.d_mat;
        delete it->second.n_mat;
        delete it->second.x_mat;
        delete it->second.cX_mat;
        delete it->second.scratch_mat1;
        delete it->second.scratch_mat2;
        delete it->second.scratch_vec;
        }
}

bool TransitionProbabilities::areTransitionProbabilitiesValid(double tolerance) {

    bool allGood = true;
    for (std::map<RbBitSet,TransitionProbabilitiesInfo>::iterator it = transProbs.begin(); it != transProbs.end(); it++)
        {
        for (int n=0; n<it->second.numMatrices; n++)
            {
            for (int i=0; i<numStates; i++)
                {
                double sum = 0.0;
                for (int j=0; j<numStates; j++)
                    sum += (*it->second.probs[activeProbs][n])(i,j);
                if (fabs(1.0 - sum) > tolerance)
                    allGood = false;
                }
            }
        }
    return allGood;
}

void TransitionProbabilities::flipActive(void) {

    activeProbs ^= 1;
}

DoubleMatrix** TransitionProbabilities::getTransitionProbabilities(RbBitSet& bs) {

    std::map<RbBitSet,TransitionProbabilitiesInfo>::iterator it = transProbs.find(bs);
    if (it == transProbs.end())
        Msg::error("Could not find transition probability vector for mask " + bs.bitString());
    return it->second.probs[activeProbs];
}

DoubleMatrix& TransitionProbabilities::getTransitionProbabilities(RbBitSet& bs, int nodeIdx) {

    std::map<RbBitSet,TransitionProbabilitiesInfo>::iterator it = transProbs.find(bs);
    if (it == transProbs.end())
        Msg::error("Could not find transition probability vector for mask " + bs.bitString());
    return *it->second.probs[activeProbs][nodeIdx];
}

void TransitionProbabilities::initialize(Model* m, ThreadPool* p, std::vector<Alignment*>& alns, int nn, int ns, int sm) {

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
        std::vector<bool> bm = alns[i]->getTaxonMask();
        RbBitSet mask(bm);

        // if the mask is not found in the map, insert it
        std::map<RbBitSet,TransitionProbabilitiesInfo>::iterator it = transProbs.find(mask);
        if (it == transProbs.end())
            {
            Tree* t = modelPtr->getTree(mask);
            if (t == NULL)
                Msg::error("Could not find tree for mask " + mask.bitString() + " when initializing transition probabilities");
            int nNodes = t->getNumNodes();
            TransitionProbabilitiesInfo info;
            info.numMatrices = nNodes;
            info.numStates = numStates;

            for (int s=0; s<2; s++)
                {
                info.probs[s] = new DoubleMatrix*[nNodes];
                for (int mi=0; mi<nNodes; mi++)
                    {
                    info.probs[s][mi] = new DoubleMatrix(numStates, numStates);
                    info.probs[s][mi]->setIdentity();
                    }
                }
            info.a_mat = new DoubleMatrix(numStates, numStates);
            info.d_mat = new DoubleMatrix(numStates, numStates);
            info.n_mat = new DoubleMatrix(numStates, numStates);
            info.x_mat = new DoubleMatrix(numStates, numStates);
            info.cX_mat = new DoubleMatrix(numStates, numStates);
            info.scratch_mat1 = new DoubleMatrix(numStates, numStates);
            info.scratch_mat2 = new DoubleMatrix(numStates, numStates);
            info.scratch_vec = new double[numStates];
            
            transProbs.insert( std::make_pair(mask,info) );
            }
        }
        
    stationaryFreqs[0].resize(numStates);
    stationaryFreqs[1].resize(numStates);
        
    isInitialized = true;
    
    //print();
}

void TransitionProbabilities::print(void) {

    std::cout << std::fixed << std::setprecision(5);
    for (std::map<RbBitSet,TransitionProbabilitiesInfo>::iterator it = transProbs.begin(); it != transProbs.end(); it++)
        {
        for (int n=0; n<it->second.numMatrices; n++)
            {
            std::cout << "Transition probabilities for node " << n << " (" << it->first.bitString() << ")" << std::endl;
            it->second.probs[activeProbs][n]->print();
            }
        }
}

void TransitionProbabilities::getStationaryFrequencies(std::vector<double>& f) {

    for (int i=0; i<numStates; i++)
        f[i] = stationaryFreqs[activeProbs][i];
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
        setTransitionProbabilitiesUsingPadeMethod();
        }
        
    needsUpdate = false;
}

void TransitionProbabilities::setTransitionProbabilitiesJc69(void) {

    // calculate transition probabilities under the Jukes-Cantor (1969) model
    for (std::map<RbBitSet,TransitionProbabilitiesInfo>::iterator it = transProbs.begin(); it != transProbs.end(); it++)
        {
        Tree* t = modelPtr->getTree(it->first);
        if (t == NULL)
            Msg::error("Could not find tree for mask " + it->first.bitString());
                
        DoubleMatrix** probs = it->second.probs[activeProbs];
        
        std::vector<Node*>& traversalSeq = t->getDownPassSequence();
        for (int n=0; n<traversalSeq.size(); n++)
            {
            Node* p = traversalSeq[n];
            DoubleMatrix& tp = *probs[p->getIndex()];
            double v = p->getBranchLength();
            
            double x = -((double)numStates/(numStates-1));
            double pChange = (1.0/numStates) - (1.0/numStates) * exp(x * v);
            double pNoChange = (1.0/numStates) + ((double)(numStates-1)/numStates) * exp(x * v);
            for (int i=0; i<numStates; i++)
                {
                for (int j=0; j<numStates; j++)
                    {
                    if (i == j)
                        tp(i,j) = pNoChange;
                    else
                        tp(i,j) = pChange;
                    }
                }
            }
        
        }
            
    double sf = 1.0 / numStates;
    for (int i=0; i<numStates; i++)
        stationaryFreqs[activeProbs][i] = sf;
}

void TransitionProbabilities::setTransitionProbabilitiesUsingPadeMethod(void) {

    RateMatrix& rmat = RateMatrix::rateMatrix();
    DoubleMatrix& Q = rmat.getRateMatrix();

    auto tasks = new TransitionProbabilitiesTask[transProbs.size()];
    auto task = tasks;
    
    // update the main tree
    for (auto it = transProbs.begin(); it != transProbs.end(); it++)
        {
        Tree* t = modelPtr->getTree(it->first);
        if (t == NULL)
            Msg::error("Could not find tree for mask " + it->first.bitString());
        task->init(t, &Q, it->second, activeProbs);
        threadPool->PushTask(task);
        ++task;
        }
        
    threadPool->Wait();
    delete[] tasks;
        
    auto& rmatFreqs = rmat.getEquilibriumFrequencies();
    for (int i=0; i<numStates; i++)
        stationaryFreqs[activeProbs][i] = rmatFreqs[i];
}

/* The method approximates the matrix exponential, P = e^A, using
   the algorithm 11.3.1, described in:

   Golub, G. H., and C. F. Van Loan. 1996. Matrix Computations, Third Edition.
      The Johns Hopkins University Press, Baltimore, Maryland.

   The method has the advantage of error control. The error is controlled by
   setting qValue appropriately (using the function SetQValue). */
void computeMatrixExponential(DoubleMatrix* Q, int qValue, double v, DoubleMatrix* A, DoubleMatrix* P, DoubleMatrix* D, DoubleMatrix* N, DoubleMatrix* X, DoubleMatrix* cX, int numStates, DoubleMatrix* scratch1, DoubleMatrix* scratch2, double* scratchVec) {
    
    // A is the matrix Q * v and p = exp(a)
    MatrixMath::multiplicationByScalar(Q, v, A);

	// set up identity matrices
    D->setIdentity();
    N->setIdentity();
    X->setIdentity();

	double maxAValue = 0.0;
	for (int i=0; i<numStates; i++)
		maxAValue = ((maxAValue > (*A)(i,i) ) ? maxAValue : (*A)(i,i) );

	int y = logBase2Plus1(maxAValue);
	int j = (( 0 > y ) ? 0 : y);
	
	MatrixMath::divideMatrixByPowerOfTwo(A, j);

	double c = 1.0;
	for (int k=1; k<=qValue; k++)
		{
		c = c * (qValue - k + 1.0) / ((2.0 * qValue - k + 1.0) * k);

		/* X = AX */
		MatrixMath::multiplyTwoMatrices(A, X, X, scratch1);

		/* N = N + cX */
		MatrixMath::multiplicationByScalar(X, c, cX);
        N->add(*cX);

		/* D = D + (-1)^k*cX */
		int negativeFactor = (k % 2 == 0 ? 1 : -1);
		if ( negativeFactor == -1 )
			MatrixMath::multiplicationByScalar(cX, negativeFactor, cX);
		D->add(*cX);
		}

	MatrixMath::gaussianElimination(D, N, P, scratch1, scratch2, scratchVec);

	for (int k=0; k<j; k++)
		MatrixMath::multiplyTwoMatrices(P, P, P, scratch1);
	
	for (int i=0; i<numStates; i++)
		{
		for (j=0; j<numStates; j++)
			{
			if ((*P)(i,j) < 0.0)
				(*P)(i,j) = fabs( (*P)(i,j) );
			}
		}
}

int logBase2Plus1(double x) {

	int j = 0;
	while(x > 1.0 - 1.0e-07)
		{
		x /= 2.0;
		j++;
		}
	return (j);
}

int setQvalue(double tol) {
	
	double x = pow(2.0, 3.0 - (0 + 0)) * factorial(0) * factorial(0) / (factorial(0+0) * factorial(0+0+1));
	int qV = 0;
	while (x > tol)
		{
		qV++;
		x = pow(2.0, 3.0 - (qV + qV)) * factorial(qV) * factorial(qV) / (factorial(qV+qV) * factorial(qV+qV+1));
		}
	return (qV);
}

double factorial(int x) {
	
	double fac = 1.0;
	for (int i=0; i<x; i++)
		fac *= (i+1);
	return (fac);
}
