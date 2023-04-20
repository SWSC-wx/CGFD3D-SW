#ifndef __TIMER_H__
#define __TIMER_H__


#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>


#define TIME(a,b) (((b).tv_sec-(a).tv_sec)*1e6+((b).tv_usec-(a).tv_usec))
#define MAXNUM 128
#define RPCC_TIME 1.5e9
#define COST_TIME_BEGIN 17000 / 8888.0
#define COST_TIME_END1 37000 / 8888.0
#define COST_TIME_TEST (COST_TIME_BEGIN + COST_TIME_END1)
#define COST_RPCC_BEGIN 106
#define COST_RPCC_END1 167
#define COST_RPCC_TEST (COST_RPCC_BEGIN + COST_RPCC_END1)
extern struct timeval t1[MAXNUM], t2[MAXNUM];
extern long timer_begin[MAXNUM], timer_end[MAXNUM];
extern long timer_total[MAXNUM];
extern long timer_count[MAXNUM];
void time_init();
void time_begin(int i);
void time_end1(int i, char str[]);
void time_end2(int i, char str[]);
void time_end3(int i, char str[]);
void time_end4(int i, int times, char str[]);
void rpcc_init();
void rpcc_begin(int i);
void rpcc_end1(int i, char str[]);
void rpcc_end2(int i, char str[]);
void rpcc_end3(int i, char str[]);
void rpcc_end4(int i, int times, char str[]);
void rpcc_cost();
int Gstart(char *s);
int Gstop(char *s);
int __real_GPTLstart(char *s);
int __real_GPTLstop(char *s);


#endif