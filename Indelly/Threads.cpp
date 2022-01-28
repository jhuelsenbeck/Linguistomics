#include <atomic> 
#include "Threads.hpp"

ThreadTask::ThreadTask() {
}

ThreadTask::~ThreadTask() {
}

void ThreadTask::Run() {
}


#if 1
// Threaded version

ThreadPool::ThreadPool():
    ThreadCount(std::thread::hardware_concurrency()),
    Running(true),
    Active(false),
    TaskCount(0),
    SleepInterval(std::chrono::microseconds(1000)),
    Threads(new std::thread[ThreadCount])
{
    for (int i = 0; i < ThreadCount; i++)
        Threads[i] = std::thread(&ThreadPool::Worker, this);
}

void ThreadPool::PushTask(ThreadTask* task) {
    const std::scoped_lock lock(TaskMutex);
    ++TaskCount;
    Tasks.push(task);
}

ThreadTask* ThreadPool::PopTask() {
    std::scoped_lock lock(TaskMutex);
    if (Tasks.empty()) 
        return NULL;
    else {
        auto task = Tasks.front();
        Tasks.pop();
        return task;
    }
}

void ThreadPool::Wait() {
    for (;;) {
        //      std::unique_lock lock(WaitMutex);
        //      WaitCondition.wait(lock);

        int count;
        {
          std::scoped_lock lock2(TaskMutex);
          count = TaskCount;
        }
        if (count == 0)
            break;
        else
            std::this_thread::sleep_for(SleepInterval);
    }
}

void ThreadPool::Worker() {
    while (Running)
    {
        ThreadTask* task = PopTask();
        if (task)
        {
            task->Run();
            delete task;
            std::scoped_lock lock(TaskMutex);
            --TaskCount;
        }
        else
            std::this_thread::sleep_for(SleepInterval);
        //        std::this_thread::yield();
    }
}

#else
// Serial version

ThreadPool::ThreadPool():
    ThreadCount(1),
    Running(false),
    Active(false),
    TaskCount(0),
    Threads(NULL)
{
}

void ThreadPool::PushTask(ThreadTask* task) {
    task->Run();
    delete task;
}

ThreadTask* ThreadPool::PopTask() {
    return NULL;
}

void ThreadPool::Wait() {
}

void ThreadPool::Worker() {
}
#endif

ThreadPool::~ThreadPool() {
    Wait();
    Running = false;
    if (Threads) 
        {
        for (auto* t = Threads; t < Threads + ThreadCount; ++t)
            t->join();
        delete[] Threads;
        }
}


