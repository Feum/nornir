#if !defined(UTILS_HPP_)
#define UTILS_HPP_

#include <sys/time.h>

static size_t elapsed_time(size_t us)
{
  struct timeval t;
  gettimeofday(&t, NULL);
  return t.tv_sec * 1000000 + t.tv_usec - us;
}

static double compute_mean(const double * times)
{
  int i;
  double sum = 0;

  for (i = 0; i < NUM_ITERS; ++i) {
    sum += times[i];
  }

  return sum / NUM_ITERS;
}

static double compute_cov(const double * times, double mean)
{
  double variance = 0;

  int i;
  for (i = 0; i < NUM_ITERS; ++i) {
    variance += (times[i] - mean) * (times[i] - mean);
  }

  variance /= NUM_ITERS;

  double cov = sqrt(variance);;
  cov /= mean;
  return cov;
}

static size_t reduce_min(long val, int id, int nprocs)
{
  static long buffer[MAX_PROCS];

  buffer[id] = val;
  //pthread_barrier_wait(&barrier);

  long min = LONG_MAX;
  int i;
  for (i = 0; i < nprocs; ++i) {
    if (buffer[i] < min) min = buffer[i];
  }

  return min;
}

#endif /* UTILS_HPP_ */
