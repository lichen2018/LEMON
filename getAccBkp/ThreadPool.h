#ifndef SEQLIB_THREAD_POOL_H
#define SEQLIB_THREAD_POOL_H

#include <stdexcept>
#include "htslib/thread_pool.h"

typedef struct {

	struct hts_tpool *pool; // The shared thread pool itself

	int qsize;    // Size of I/O queue to use for this fp

}htsThreadPool1;

class ThreadPool {

 public: 
  
 ThreadPool() : nthreads(1) { p.pool = NULL; }

 ThreadPool(int n) : nthreads(1) {
    p.pool = NULL;
    if (n < 1)
      throw std::invalid_argument( "n threads must be > 0");
    if (!(p.pool = hts_tpool_init(n))) 
      throw std::runtime_error( "Error creating thread pool");
  }

  bool IsOpen() { return p.pool != NULL; }

  htsThreadPool1 p;
  size_t nthreads;

};

#endif
