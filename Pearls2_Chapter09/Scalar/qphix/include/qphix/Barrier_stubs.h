#ifndef QPHIX_BARRIER_MIC_STUBS_H
#define QPHIX_BARRIER_MIC_STUBS_H

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iostream>
using namespace std;

namespace QPhiX {

class Barrier {
public:

  Barrier(int num_threads, int threads_per_core) {}
  void init(int tid) { }
  ~Barrier() {}
  void wait(int tid) {}
 private:
  int dummy;
};

}; // end namespace

#endif
