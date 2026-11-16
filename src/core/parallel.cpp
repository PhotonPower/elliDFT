#include "core/parallel.hpp"

#ifdef _OPENMP
#include <omp.h>
#endif

namespace ellidft {

#ifdef _OPENMP
bool parallel_enabled() { return true; }
int max_threads() { return omp_get_max_threads(); }
void set_num_threads(int n) { omp_set_num_threads(n < 1 ? 1 : n); }
int thread_index() { return omp_get_thread_num(); }
#else
bool parallel_enabled() { return false; }
int max_threads() { return 1; }
void set_num_threads(int) {}
int thread_index() { return 0; }
#endif

}  // namespace ellidft
