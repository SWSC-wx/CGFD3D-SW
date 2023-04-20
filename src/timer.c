#include <stdio.h>
#include <gptl.h>
#include <mpi.h>
#include "timer.h"

struct timeval t1[MAXNUM], t2[MAXNUM];
long timer_begin[MAXNUM], timer_end[MAXNUM];
long timer_total[MAXNUM];
long timer_count[MAXNUM];

int Gstart(char *s)
{
    int ret;
    // MPI_Barrier(MPI_COMM_WORLD);
    ret = GPTLstart(s);
    return ret;
}

int Gstop(char *s)
{
    int ret;
    // MPI_Barrier(MPI_COMM_WORLD);
    ret = GPTLstop(s);
    return ret;
}

int __attribute__((noinline)) __wrap_GPTLstart(char *s) 
{
    int ret;
    // MPI_Barrier(MPI_COMM_WORLD);
    ret = __real_GPTLstart(s);
    return ret;
}

int  __attribute__((noinline)) __wrap_GPTLstop(char *s) 
{
    int ret;
    ret = __real_GPTLstop(s);
    return ret;
}

inline void time_init() {
  int i;
  for(i = 0 ; i < MAXNUM ; i ++) {
    timer_total[i] = 0;
    timer_count[i] = 0;
  }
}

inline void time_begin(int i) {
  timer_count[i] ++;
  gettimeofday(&t1[i], NULL);
}

inline void time_end1(int i, char str[]) {
  gettimeofday(&t2[i], NULL);
  timer_total[i] += TIME(t1[i], t2[i]);
}

inline void time_end2(int i, char str[]) {
  gettimeofday(&t2[i], NULL);
  timer_total[i] += TIME(t1[i], t2[i]);
  printf("Timer timer_index: %d, time: %lf, timer_total: %ld, %s\n", i, TIME(t1[i], t2[i]), timer_total[i], str);
}

inline void time_end3(int i, char str[]) {
  long t = timer_total[i] - COST_TIME_TEST * timer_count[i];
  printf("Timer timer_index: %d, total_time: %fms, total_time_real: %ldms, count: %ld, average_time: %fms %s", i, timer_total[i] * 1.0 / 1e3, t, timer_count[i], t * 1.0 / timer_count[i] / 1e3, str);
}

inline void time_end4(int i, int times, char str[]) {
  long t = timer_total[i] - COST_TIME_TEST * timer_count[i];
  printf("Timer timer_index: %d, total_time: %fms, total_time_real: %fms, count: %ld, average_time: %fms, time_per_step: %fms, %s\n", i, timer_total[i] * 1.0 / 1e3, t * 1.0 / 1e3, timer_count[i], t * 1.0 / timer_count[i] / 1e3, t * times * 1.0 / timer_count[i] / 1e3, str);
}

inline void rpcc_init() {
  int i;
  for(i = 0 ; i < MAXNUM ; i ++) {
    timer_total[i] = 0;
    timer_count[i] = 0;
  }
}

inline void rpcc_begin(int i) {
  long rpcc;
#ifdef ARCH_SW
  asm volatile("rtc %0": "=r" (rpcc) : );
#elif ARCH_MPE
  asm volatile("rtc %0": "=r" (rpcc) : );
#else
  long low32, high32;
  __asm__ ("RDTSC" : "=a"(low32),"=d"(high32));
  rpcc = (high32 << 32) + low32;
#endif
  timer_begin[i] = rpcc;
  timer_count[i] ++;
}

inline void rpcc_end1(int i, char str[]) {
  long rpcc;
#ifdef ARCH_SW
  asm volatile("rtc %0": "=r" (rpcc) : );
#elif ARCH_MPE
  asm volatile("rtc %0": "=r" (rpcc) : );
#else
  long low32, high32;
  __asm__ volatile("RDTSC" : "=a"(low32),"=d"(high32));
  rpcc = (high32 << 32) + low32;
#endif
  timer_total[i] += rpcc - timer_begin[i];
}

inline void rpcc_end2(int i, char str[]) {
  long rpcc;
#ifdef ARCH_SW
  asm volatile("rtc %0": "=r" (rpcc) : );
#elif ARCH_MPE
  asm volatile("rtc %0": "=r" (rpcc) : );
#else
  long low32, high32;
  __asm__ volatile("RDTSC" : "=a"(low32),"=d"(high32));
  rpcc = (high32 << 32) + low32;
#endif
  timer_end[i] = rpcc;
  timer_total[i] += timer_end[i] - timer_begin[i];
  printf("Timer timer_index: %d, rpcc: %ld, total_rpcc: %ld, count: %ld, average_rpcc: %f,  average_time: %f, %s\n", i, timer_end[i] - timer_begin[i], timer_total[i], timer_count[i], timer_total[i] * 1.0 / timer_count[i], timer_total[i] * 1.0 / timer_count[i] / RPCC_TIME, str);
}

inline void rpcc_end3(int i, char str[]) {
  printf("Timer timer_index: %d, total_rpcc: %ld, count: %ld, average_rpcc: %f,  average_time: %f, %s\n", i, timer_total[i], timer_count[i], timer_total[i] * 1.0 / timer_count[i], timer_total[i] * 1.0 / RPCC_TIME, str);
}

inline void rpcc_end4(int i, int times, char str[]) {
  long t = timer_total[i] - COST_RPCC_TEST * timer_count[i];
  printf("Timer timer_index: %d, total_rpcc: %ld, total_rpcc_real: %ld, count: %ld, average_rpcc: %f,  total_time: %f, time_per_step: %fms, %s\n", i, timer_total[i], t, timer_count[i], t * 1.0 / timer_count[i], t * 1.0 / RPCC_TIME, t * times * 1000.0 / timer_count[i] / RPCC_TIME, str);
}

inline void rpcc_cost() {
  int count = 0;
  int i;
  for(i = 0 ; i < MAXNUM ; i ++) {
    count += timer_count[i];
  }
  //printf("Timer cost rpcc: %d, count: %d, average_time: %f, %s\n", count * (15 + 15), count, count * (15 + 15) / RPCC_TIME);
  printf("Timer cost rpcc: %d, count: %d, average_time: %f, %s\n", count * (106 + 167), count, count * (106 + 167) / RPCC_TIME, "timer cost");
}
//master: asm("rtc %0": "=r" (rpcc) : );
//slave: asm volatile("rcsr  %0, 4":"=r"(rpcc));
//sw: 1.5G cycles / 1 second
//sw small machine : 1.3G cycles / 1 second
