#include "time.h"

#ifndef BENCH_LOG_SIZE
#define BENCH_LOG_SIZE 128
#endif

__uint64_t bench_itt = 0;
double bench_times[BENCH_LOG_SIZE];

#define START_BENCH struct timespec start, end;\
                    clock_gettime(CLOCK_MONOTONIC_RAW, &start);\

#define LOG_BENCH   clock_gettime(CLOCK_MONOTONIC_RAW, &end);\
                    double sec = (double)(end.tv_sec - start.tv_sec) + 1e-9 * (double)(end.tv_nsec - start.tv_nsec);\
                    bench_times[bench_itt % BENCH_LOG_SIZE] = sec * 1e3;\
                    bench_itt++;

#define PRINT_BENCH {  \
                            double mean = 0.;\
                            double std = 0.;\
                            for(int i = 0; i < BENCH_LOG_SIZE; i++){\
                                mean += bench_times[i]; \
                                std += bench_times[i]*bench_times[i];\
                            }\
                            mean /= BENCH_LOG_SIZE;\
                            std /= BENCH_LOG_SIZE;\
                            std = sqrt(std - mean*mean);\
                            printf("bench time : %9.4lf ms (± %9.4lf) - iteration : %lu\n", mean, std, bench_itt);\
                        }

#define PRINT_BENCH_PERIOD(t) if(bench_itt % t == 0) PRINT_BENCH
