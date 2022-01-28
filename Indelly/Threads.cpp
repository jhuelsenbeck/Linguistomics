#include <atomic> 
#include "threads.hpp"

ThreadTask::ThreadTask() {
}

ThreadTask::~ThreadTask() {
}

void ThreadTask::Run() {
}



ThreadPool::ThreadPool():
    ThreadCount(std::thread::hardware_concurrency()),
    Running(true),
    TaskCount(0)
{
    Threads = new std::thread[ThreadCount];
    for (int i = 0; i < ThreadCount; i++)
    {
        Threads[i] = std::thread(&ThreadPool::Worker, this);
    }
}

ThreadPool::~ThreadPool() {
    Wait();
    Running = false;
    for (auto* t = Threads; t < Threads + ThreadCount; ++t)
        t->join();
    delete[] Threads;
}

void ThreadPool::Wait() {
    while (TaskCount > 0)
      std::this_thread::yield();
}

void ThreadPool::PushTask(ThreadTask* task) {
    const std::scoped_lock lock(TaskMutex);
    TaskCount++;
    Tasks.push(task);
}

ThreadTask* ThreadPool::PopTask() {
    const std::scoped_lock lock(TaskMutex);
    if (Tasks.empty())
        return NULL;
    else
    {
        auto task = Tasks.front();
        Tasks.pop();
        return task;
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
            TaskCount--;
        }
        else
          std::this_thread::yield();
    }
}
