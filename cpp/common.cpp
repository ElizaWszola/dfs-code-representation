#include "common.hpp"

int64_t get_time_difference(clocktime *start, clocktime *end) {
  int64_t seconds_diff = (int64_t)(end->tv_sec)
      - (int64_t)(start->tv_sec);
  int64_t end_nano = end->tv_nsec;
  int64_t start_nano = start->tv_nsec;
  if (start_nano <= end_nano)
    return seconds_diff * GIGA + end_nano - start_nano;
  return (seconds_diff - 1) * GIGA + GIGA - (start_nano - end_nano);
}

int get_time(clocktime *time) {
    return clock_gettime(CLOCK_REALTIME, time);
}