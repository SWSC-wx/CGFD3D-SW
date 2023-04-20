#include"sort_int.h"
#include<stdio.h>
#include<stdlib.h>
#include<time.h>
int comp(const Data *a, const Data *b) {
  return (*a).x > (*b).x;
}

void sort_int(int a[], int index[], int n){
  Data *dat = (Data *)malloc(sizeof(Data) * n);
  int i;
  for(i = 0 ; i < n ; i ++) {
    dat[i].x = a[i];
    dat[i].y = index[i];
  }
  qsort(&dat[0], n, sizeof(dat[0]), comp);
  for(i = 0 ; i < n ; i ++) {
    a[i] = dat[i].x;
    index[i] = dat[i].y;
  }
  free(dat);
}

/*
int main() {
  data dat[N];
  srand(unsigned(time(0)));
  int a[N], index[N];
  for(int i = 0 ; i < N ; i ++) {
    a[i] = rand();
    index[i] = i;
  }
  bubble_sort_int(a, index, N);
  for(int i = 0 ; i < N ; i ++) {
    printf("a[%d] = %d, index[%d] = %d\n", i, a[i], i, index[i]);
  }
}
*/
