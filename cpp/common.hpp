#ifndef COMMON_DFSR_HPP
#define COMMON_DFSR_HPP

#include <ctime>
#include <cstdlib>
#include <limits>


#define GIGA 1000000000
typedef struct timespec clocktime;

int64_t get_time_difference(clocktime *start, clocktime *end);
int get_time(clocktime *time);

#endif
