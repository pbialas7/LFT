#pragma once

#include <condition_variable>
#include <functional>
#include <future>
#include <mutex>
#include <queue>
#include <stdexcept>
#include <thread>
#include <vector>

class ThreadPool {
public:
    // Construct pool with `num_threads` worker threads (defaults to hardware concurrency)
    explicit ThreadPool(size_t num_threads = std::thread::hardware_concurrency())
        : stop_(false)
    {
        if (num_threads == 0)
            num_threads = 1;

        workers_.reserve(num_threads);
        for (size_t i = 0; i < num_threads; ++i) {
            workers_.emplace_back([this] { worker_loop(); });
        }
    }

    // Non-copyable, non-movable
    ThreadPool(const ThreadPool&)            = delete;
    ThreadPool& operator=(const ThreadPool&) = delete;

    // Destructor: finish all queued tasks, then join all threads
    ~ThreadPool() {
        {
            std::unique_lock<std::mutex> lock(queue_mutex_);
            stop_ = true;
        }
        cv_.notify_all();
        for (std::thread& t : workers_)
            t.join();
    }

    // Submit a callable with any arguments.
    // Returns a std::future holding the return value (or exception).
    //
    // Usage:
    //   auto fut = pool.submit(my_func, arg1, arg2);
    //   auto result = fut.get();   // blocks until done
    template<typename F, typename... Args>
    auto submit(F&& func, Args&&... args)
        -> std::future<std::invoke_result_t<F, Args...>>
    {
        using return_type = std::invoke_result_t<F, Args...>;

        auto task = std::make_shared<std::packaged_task<return_type()>>(
            std::bind(std::forward<F>(func), std::forward<Args>(args)...)
        );

        std::future<return_type> fut = task->get_future();

        {
            std::unique_lock<std::mutex> lock(queue_mutex_);
            if (stop_)
                throw std::runtime_error("submit() called on a stopped ThreadPool");

            task_queue_.emplace([task]() { (*task)(); });
        }
        cv_.notify_one();
        return fut;
    }

    // Number of worker threads
    size_t size() const { return workers_.size(); }

    // Number of tasks currently waiting in the queue
    size_t queue_size() const {
        std::unique_lock<std::mutex> lock(queue_mutex_);
        return task_queue_.size();
    }

private:
    // Each worker runs this loop until stop_ is set and the queue is empty
    void worker_loop() {
        while (true) {
            std::function<void()> task;
            {
                std::unique_lock<std::mutex> lock(queue_mutex_);
                cv_.wait(lock, [this] {
                    return stop_ || !task_queue_.empty();
                });

                if (stop_ && task_queue_.empty())
                    return;

                task = std::move(task_queue_.front());
                task_queue_.pop();
            }
            task();   // execute outside the lock
        }
    }

    std::vector<std::thread>          workers_;
    std::queue<std::function<void()>> task_queue_;
    mutable std::mutex                queue_mutex_;
    std::condition_variable           cv_;
    bool                              stop_;
};
