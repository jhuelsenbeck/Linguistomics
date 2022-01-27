#include <atomic>      // std::atomic
#include <functional>  // std::function
#include <future>      // std::future, std::promise
#include <utility>     // std::move
#include "threads.hpp" // std::move


thread_pool::thread_pool()
    {
    thread_count  = std::thread::hardware_concurrency();
    sleepinterval = std::chrono::microseconds(1000);
    running       = true;
    tasks_total   = 0;
    threads       = new std::thread[thread_count];
    for (int i = 0; i < thread_count; i++)
        threads[i] = std::thread(&thread_pool::worker, this);
}

thread_pool::~thread_pool()
{
    wait_for_tasks();
    running = false;
    for (auto *t = threads; t < threads + thread_count; ++t)
        t->join();
    delete[] threads;
}

int thread_pool::get_tasks_running()
{
    const std::scoped_lock lock(queue_mutex);
    return tasks_total - (int)tasks.size();
}


void thread_pool::wait_for_tasks()
{
    while (true)
    {
        if (tasks_total == 0)
            break;
        sleep_or_yield();
    }
}

bool thread_pool::pop_task(std::function<void()> &task)
{
    const std::scoped_lock lock(queue_mutex);
    if (tasks.empty())
        return false;
    else
    {
        task = std::move(tasks.front());
        tasks.pop();
        return true;
    }
}

void thread_pool::sleep_or_yield()
{
    std::this_thread::sleep_for(sleepinterval);
}

void thread_pool::worker()
{
    while (running)
    {
        std::function<void()> task;
        if (pop_task(task))
        {
            task();
            tasks_total--;
        }
        else
        {
            sleep_or_yield();
        }
    }
}
