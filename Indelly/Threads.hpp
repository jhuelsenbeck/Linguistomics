#ifndef threads_hpp
#define threads_hpp

#include <thread>
#include <mutex> 
#include <queue>
#include "Container.hpp"


class ThreadCache {
    public:
        ThreadCache(int numstates);
        ~ThreadCache();

        DoubleMatrix* scratch1;
        DoubleMatrix* scratch2;
        double* scratchVec;
};


class ThreadTask {
    public:
        ThreadTask();
        virtual ~ThreadTask() {};
        virtual void Run(ThreadCache& cache);
};


class ThreadPool {
    public:
        int  ThreadCount;

             explicit ThreadPool(int numstates);
             ~ThreadPool();

        void PushTask(ThreadTask* task);
        void Wait();

    private:
        std::atomic<size_t>     TaskCount;
        std::atomic<bool>       Running;
        std::mutex              TaskMutex,
                                WaitMutex,
                                CheckMutex;
        std::condition_variable WaitCondition,
                                CheckCondition;
        std::queue<ThreadTask*> Tasks;
        std::thread*            Threads;
        int                     numStates;

        void                    Worker();
        ThreadTask*             PopTask();
};


#endif
