#ifndef CORRELATION_HPP
#define CORRELATION_HPP

#include <iostream>
#include <cmath>
#include <queue>
#include <cstring>
#include <vector>
#include <chrono>
#include <future>
#include "graph_fast.hpp"

#define JOB_STEALING false
typedef float c_real;

void cmc_part_rec(std::vector<Graph>& x, c_real* y, bool* b,
    c_real& best_corr_global, c_real* best_corr_local, uint64_t* local_b,
    double delta, int64_t* max_queue, int64_t* total_processed);

void cmc_part_rec_subcall(std::vector<Graph>& x, c_real* y, bool* b,
    DFSCode& code, std::vector<PDFS>& candidates,
    std::vector<uint32_t>& rmpath, std::vector<History>& histories,
    EdgeList& edges, c_real& best_corr_global, c_real* best_corr_local,
    uint64_t* local_b, double delta, c_real curr_mu,
    int64_t& max_queue, int64_t sub_max_queue, int64_t& total_processed);

// ThreadPool will be initialized with parameterized number of threads,
// but currently not all might be active during the execution.
// To get max_corr, min(n_threads, features) are used
class ThreadPool {

private:
  
  struct FindMaxJobRec {
    std::vector<Graph>& x;
    c_real* y;
    bool* b;
    c_real& best_corr_global;
    c_real* best_corr_local;
    uint64_t* local_b;
    double delta;
    int64_t* max_queue;
    int64_t* total_processed;
    void call() {
      cmc_part_rec(x, y, b, best_corr_global, best_corr_local,
                   local_b, delta, max_queue, total_processed);
    }
    //~ #if JOB_STEALING
    //~ void launder(uint64_t* local_b, c_real* best_corr_thread,
                //~ int64_t* max_queue, int64_t* total_processed) {
      //~ this->local_b = local_b;
      //~ this->best_corr_thread = best_corr_thread;
      //~ this->max_queue = max_queue;
      //~ this->total_processed = total_processed;
    //~ }
    //~ #endif
    FindMaxJobRec(std::vector<Graph>& x, c_real* y, bool* b,
                  c_real& best_corr_global, c_real* best_corr_local,
                  uint64_t* local_b, double delta, int64_t* max_queue,
                  int64_t* total_processed):
                      x(x), y(y), b(b),
                      best_corr_global(best_corr_global),
                      best_corr_local(best_corr_local),
                      local_b(local_b), delta(delta),
                      max_queue(max_queue),
                      total_processed(total_processed) {}
        
  };

  int32_t n_threads;
  int32_t jobs_done;
  bool* job_done;
  std::vector<std::thread> threads;
  uint64_t* local_bs;
  bool running;
  c_real* best_corr_threads;

  std::queue<FindMaxJobRec*>* jobs_rec;
  
  c_real** mus;
  bool stopped;
  
  int64_t* max_queues;
  int64_t* total_processeds;
  
  int32_t nt_max_corr;
  
  #if JOB_STEALING
  struct SimpleLock {
    std::atomic<bool> locked;
    SimpleLock() {
      locked = false;
    }
    void lock() {
      bool expect;
      do {
        expect = false;
      } while(!locked.compare_exchange_strong(expect, true, std::memory_order_acq_rel));
    }
    void unlock() {
      bool expect;
      do {
        expect = true;
      } while(!locked.compare_exchange_strong(expect, false, std::memory_order_acq_rel));
    }
  };
  
  //~ std::vector<std::mutex> qmutex;
  std::vector<SimpleLock> qmutex;
  #endif
  
  
  
public:

  bool is_running() {
    return running;
  }

  c_real cmc_full(std::vector<Graph>& x, c_real* y, bool* b, bool* fb,
      double delta, int64_t& max_queue, int64_t& total_processed) {

    size_t n = x.size();
    c_real best_corr = 0;
    uint32_t featbits = 1; //trash
    
    // playing it safe for now
    nt_max_corr = std::min(n_threads, 1);

    jobs_done = 0;
    
    if (!local_bs) {
      local_bs = new uint64_t[nt_max_corr * featbits];
      std::memset(local_bs, 0, sizeof(uint64_t) * nt_max_corr * featbits);
    }
    
    for (int32_t i = 0; i < nt_max_corr; ++i) {
      job_done[i] = false;
      best_corr_threads[i] = 0.0;
      max_queues[i] = 0;
      total_processeds[i] = 0;
    }
    std::vector<std::future<void>> results;
    std::vector<FindMaxJobRec*> local_jobs_rec;
    for (uint32_t i = 0; i < n; ++i) {
      int32_t thread_no = i % nt_max_corr;
      FindMaxJobRec* j = new FindMaxJobRec(x, y, b, best_corr,
          best_corr_threads + thread_no, local_bs + thread_no * featbits,
          delta, max_queues + thread_no, total_processeds + thread_no);
      local_jobs_rec.push_back(j);
      jobs_rec[thread_no].push(j);
    }
    for (int32_t i = 0; i < nt_max_corr; ++i) {
      results.push_back(std::async([=]{ return run_find_max_rec(i); }));
    }
    
    do {
      c_real new_best_corr = 0.0;
      for (int32_t i = 0; i < nt_max_corr; ++i) {
        c_real bct = best_corr_threads[i];
        if (std::fabs(bct) > std::fabs(new_best_corr)) {
          new_best_corr = bct;
        }
      }
      if (std::fabs(best_corr) < std::fabs(new_best_corr)) {
        best_corr = new_best_corr;
      }
      jobs_done = 0;
      for (int32_t i = 0; i < nt_max_corr; ++i) {
        jobs_done += job_done[i];
      }
      std::this_thread::sleep_for(std::chrono::nanoseconds(100));
    } while (jobs_done != nt_max_corr);
    
    for (int32_t i = 0; i < nt_max_corr; ++i) {
      results[i].get();
    }
    
    c_real new_best_corr = 0.0;
    for (int32_t i = 0; i < nt_max_corr; ++i) {
      c_real bct = best_corr_threads[i];
      if (std::fabs(bct) > std::fabs(new_best_corr)) {
        new_best_corr = bct;
      }
    }
    if (std::fabs(best_corr) < std::fabs(new_best_corr)) {
      best_corr = new_best_corr;
    }

    for (int32_t i = 0; i < nt_max_corr; ++i) {
      if (best_corr_threads[i] == best_corr) {
        std::memcpy(b, local_bs + i * featbits, featbits * sizeof(uint64_t));
      }
    }
    
        //~ std::cout << "BEST: " << best_corr << "\n";
        
    //todo: subgraph check in this one... somehow
    
    //~ for (int64_t i = 0; i < samples; i++) {
      //~ int32_t j = 0;
      //~ for (; j < featbits; ++j) {
        //~ if ((x[i * rowbits + j] & b[j]) != b[j]) {
          //~ break;
        //~ }
      //~ }
        //~ fb[i] = (j == featbits);
              //~ std::cout << fb[i];
    //~ }
    //~ std::cout << std::endl;
    
    //~ for (int64_t i = 0; i < samples; ++i) {
      //~ std::cout << y[i] << " ";
    //~ }
    //~ std::cout << "\n";
    
    //~ #if GET_STATS
    //~ int64_t new_max_queue = 0;
    //~ int64_t new_total_processed = 0;
    //~ for (int32_t i = 0; i < nt_max_corr; ++i) {
      //~ new_max_queue = std::max(new_max_queue, max_queues[i]);
      //~ new_total_processed += total_processeds[i];
    //~ }
    //~ max_queue = new_max_queue;
    //~ total_processed = new_total_processed;
    //~ #endif
    
    // will this suffer much from false sharing? ---> could be, this is actually pretty slow
    //~ int32_t spt = (samples + nt_update - 1) / nt_update;
    //~ std::vector<std::future<void>> results_update;
    //~ for (int32_t i = 0; i < nt_update; ++i) {
      //~ results_update.push_back(std::async([=]{
          //~ return update_freq_range(x, b, fb, i * spt, std::min(i * spt + spt, samples), featbits);
      //~ }));
    //~ }
    //~ for (int32_t i = 0; i < nt_update; ++i) {
      //~ results_update[i].get();
    //~ }
    
    for (FindMaxJobRec* job : local_jobs_rec) {
      delete job;
    }
    
    return best_corr;
  }
  
  c_real cmc_full(std::vector<Graph>& x, c_real* y, bool* b, bool* fb) {
    int64_t dummy_max;
    int64_t dummy_total;
    return cmc_full(x, y, b, fb, 0.0, dummy_max, dummy_total);
  }

  void run_find_max_rec(int32_t id) {
    FindMaxJobRec* function = nullptr;
    while (!jobs_rec[id].empty()) {
      #if JOB_STEALING
      qmutex[id].lock();
      if (!jobs_rec[id].empty()) {
      #endif
        function = jobs_rec[id].front();
        jobs_rec[id].pop();
      #if JOB_STEALING
      }
      qmutex[id].unlock();
      #endif
      if (function) {
        function->call();
      }
    }
    job_done[id] = true;
    #if JOB_STEALING
    bool ndone = true;
    while (ndone) {
      ndone = false;
      for (int32_t i = 0; i < nt_max_corr; ++i) {
        if (i != id) {
          qmutex[i].lock();
          if (!jobs_rec[i].empty()) {
            function = jobs_rec[i].front();
            jobs_rec[i].pop();
          }
          qmutex[i].unlock();
          if (function) {
            ndone = true;
            function->launder(local_bs + id * featbits, best_corr_threads + id,
                 max_queues + id, total_processeds + id);
            function->call();
            function = nullptr;
          }
        }
      }
    }
    #endif
  }
  
  void stop() {
    if (!stopped) {
      delete[] best_corr_threads;
      delete[] jobs_rec;
      delete[] job_done;
      delete[] max_queues;
      delete[] total_processeds;
      if (local_bs) {
        delete[] local_bs;
        local_bs = nullptr;
      }
    }
    stopped = true;
  }
  
  ~ThreadPool() {
    stop();
  }

  ThreadPool(int32_t n) : n_threads(n) {
    best_corr_threads = new c_real[n_threads];
    jobs_rec = new std::queue<FindMaxJobRec*>[n_threads];
    job_done = new bool[n_threads];
    local_bs = nullptr;
    running = true;
    stopped = false;
    max_queues = new int64_t[n_threads];
    total_processeds = new int64_t[n_threads];
  };
  

};

#endif
