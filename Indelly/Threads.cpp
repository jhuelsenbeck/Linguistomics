#include <atomic>      // std::atomic
#include <future>      // std::future, std::promise
#include <utility>     // std::move
#include "threads.hpp" // std::move


ThreadPool::ThreadPool() {
    thread_count  = std::thread::hardware_concurrency();
    sleepinterval = std::chrono::microseconds(1000);
    running       = true;
    tasks_total   = 0;
    threads       = new std::thread[thread_count];
    for (int i = 0; i < thread_count; i++)
        {
        threads[i] = std::thread(&ThreadPool::worker, this);
        }
}

ThreadPool::~ThreadPool() {
    wait_for_tasks();
    running = false;
    for (auto* t = threads; t < threads + thread_count; ++t)
        {
        t->join();
        }
    delete[] threads;
}

int ThreadPool::get_tasks_running() {
    const std::scoped_lock lock(queue_mutex);
    return tasks_total - (int)tasks.size();
}


void ThreadPool::wait_for_tasks() {
    while (true)
        {
        if (tasks_total == 0)
            break;
        sleep_or_yield();
        }
}

void ThreadPool::push_task(ThreadTask* task) {
    const std::scoped_lock lock(queue_mutex);
    tasks_total++;
    tasks.push(task);
}

ThreadTask* ThreadPool::pop_task() {
    const std::scoped_lock lock(queue_mutex);
    if (tasks.empty())
        return NULL;
    else
        {
            auto task = std::move(tasks.front());
            tasks.pop();
            return task;
        }
}

void ThreadPool::sleep_or_yield() {
    std::this_thread::sleep_for(sleepinterval);
}

void ThreadPool::worker() {
    while (running)
        {
        ThreadTask* task = pop_task();
        if (task)
            {
            task->run();
            tasks_total--;
            }
        else
            {
            sleep_or_yield();
            }
        }
}
