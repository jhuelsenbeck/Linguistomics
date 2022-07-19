#ifndef Parameter_hpp
#define Parameter_hpp

#include <string>
class Model;
class RandomVariable;



class Parameter {

    public:
                            Parameter(void) = delete;
                            Parameter(RandomVariable* r, Model* m, std::string n);
        virtual            ~Parameter(void) { }
        virtual void        accept(void) = 0;
        virtual void        fillParameterValues(double* x, int& start, int maxNumValues) = 0;
        virtual std::string getJsonString(void) = 0;
        virtual std::string getHeader(void) = 0;
        std::string         getName(void) { return parmName; }
        std::string         getLastUpdate(void) { return lastUpdateType; }
        virtual int         getNumValues(void) = 0;
        double              getProposalProbability(void) { return proposalProbability; }
        virtual std::string getString(void) = 0;
        bool                getUpdateChangesRateMatrix(void) { return updateChangesRateMatrix; }
        bool                getUpdateChangesTransitionProbabilities(void) { return updateChangesTransitionProbabilities; }
        virtual double      lnPriorProbability(void) = 0;
        virtual void        print(void) = 0;
        virtual void        reject(void) = 0;
        void                setLastUpdate(std::string s) { lastUpdateType = s; }
        void                setProposalProbability(double x) { proposalProbability = x; }
        void                setUpdateChangesRateMatrix(bool tf) { updateChangesRateMatrix = tf; }
        void                setUpdateChangesTransitionProbabilities(bool tf) { updateChangesTransitionProbabilities = tf; }
        virtual double      update(int iter) = 0;
        
    protected:
        std::string         parmName;
        Model*              modelPtr;
        RandomVariable*     rv;
        std::string         lastUpdateType;
        double              proposalProbability;
        bool                updateChangesRateMatrix;
        bool                updateChangesTransitionProbabilities;
};

#endif
