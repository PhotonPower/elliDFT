#pragma once

// Thin wrapper around the OpenMP runtime, so that the rest of the code and
// the tests do not depend on omp.h. Without OpenMP everything runs in one
// thread.

namespace ellidft {

// True if the library was built with OpenMP.
bool parallel_enabled();

// Number of threads that parallel regions will use.
int max_threads();

// Sets the number of threads for subsequent parallel regions (at least 1).
// No effect without OpenMP.
void set_num_threads(int n);

// Index of the calling thread inside a parallel region, 0 outside.
int thread_index();

}  // namespace ellidft
