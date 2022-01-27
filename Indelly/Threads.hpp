#pragma once

#include <thread>
#include <mutex> 
#include <queue> 

class thread_pool {
  public:
      int         thread_count;

                  thread_pool();
                  ~thread_pool();
                  int  get_tasks_running();
                  void sleep_or_yield();
                  void worker();
                  void wait_for_tasks();
                  bool pop_task(std::function<void()>& task);

                  template <typename F>
                  void push_task(const F& task)
                  {
                    const std::scoped_lock lock(queue_mutex);
                    tasks_total++;
                    tasks.push(std::function<void()>(task));
                  }

                  template <typename F, typename... A>
                  void push_task(const F& task, const A &...args)
                  {
                      push_task([task, args...]
                          { task(args...); });
                  }


private:
                  std::chrono::microseconds         sleepinterval;
                  std::thread*                      threads;
                  std::atomic<bool>                 running;
                  std::mutex                        queue_mutex = {};
                  std::queue<std::function<void()>> tasks = {};
                  int                               tasks_total;
};

