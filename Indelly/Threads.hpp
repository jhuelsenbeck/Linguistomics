#ifndef threads_hpp
#define threads_hpp

#include <thread>
#include <mutex> 
#include <queue> 


class ThreadTask {
    public:
        ThreadTask() {
        }

        virtual void run() = 0;
};

class ThreadPool {
    public:
        int  thread_count;

            ThreadPool();
            ~ThreadPool();
        void push_task(ThreadTask* task);
        void wait_for_tasks();

    private:
        int  get_tasks_running();
        void sleep_or_yield();
        void worker();
        ThreadTask* pop_task();
    
        std::chrono::microseconds  sleepinterval;
        std::thread*               threads;
        std::atomic<bool>          running;
        std::mutex                 queue_mutex = {};
        std::queue<ThreadTask*>    tasks = {};
        int                        tasks_total;
};

#endif
