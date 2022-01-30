#include "Threads.hpp"

ThreadTask::ThreadTask() {
}

void ThreadTask::Run() {
}


#if 1
// Threaded version

ThreadPool::ThreadPool():
    ThreadCount(std::thread::hardware_concurrency()),
    TaskCount(0),
    Running(true),
    Threads(new std::thread[ThreadCount])
{
    for (int i = 0; i < ThreadCount; i++)
        Threads[i] = std::thread(&ThreadPool::Worker, this);
}

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

void ThreadPool::PushTask(ThreadTask* task) {
    const std::scoped_lock lock(TaskMutex);
    ++TaskCount;
    Tasks.push(task);
    CheckCondition.notify_one();
}

ThreadTask* ThreadPool::PopTask() {
    {
        std::unique_lock mlock(CheckMutex);
        CheckCondition.wait(mlock, [this]{return TaskCount > 0;});
    }

    std::scoped_lock lock(TaskMutex);
    if (Tasks.empty())
        return NULL;
    else 
        {
        auto task = Tasks.front();
        Tasks.pop();
        return task;
        }
}

void ThreadPool::Wait() {
    for (;;) 
        {
        // This WaitCondition is signaled when all tasks are completed
        std::unique_lock lock(WaitMutex);
        WaitCondition.wait(lock);

        // According to the documentation, this needs to be double-checked because of "spurious" wakeups

        int count;
        {
            std::scoped_lock lock2(TaskMutex);
            count = TaskCount;
        }
        if (count == 0)
            break;
        else
            std::this_thread::yield();
        }
}

void ThreadPool::Worker() {
    while (Running)
        {
        ThreadTask* task = PopTask();
        if (task)
            {
            task->Run();

            int count;
            {
              std::scoped_lock lock(TaskMutex);
              count = --TaskCount;
            }
            if (count == 0)
                WaitCondition.notify_one();
            }
        else
            std::this_thread::yield();
    }
}

#else
// Serial version

ThreadPool::ThreadPool():
    ThreadCount(1),
    Running(false),
    TaskCount(0),
    Threads(NULL)
{
}

ThreadPool::~ThreadPool() {
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

