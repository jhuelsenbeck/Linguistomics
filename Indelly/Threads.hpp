#ifndef threads_hpp
#define threads_hpp

#include <thread>
#include <mutex> 
#include <queue>


class ThreadTask {
    public:
        ThreadTask();
        virtual ~ThreadTask();
        virtual void Run();
};

class ThreadPool {
    public:
        int  ThreadCount;

             explicit ThreadPool();
             ~ThreadPool();
        void PushTask(ThreadTask* task);
        void Wait();

    private:
        std::atomic<bool>       Running;
        std::mutex              TaskMutex,
                                WaitMutex,
                                CheckMutex;
        std::condition_variable WaitCondition,
                                CheckCondition;
        std::queue<ThreadTask*> Tasks;
        std::thread*            Threads;
        int                     TaskCount;

        void              Worker();
        ThreadTask*       PopTask();
};


#endif
