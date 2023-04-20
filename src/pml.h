#ifndef H_PML
#define H_PML

inline float cal_pml_R(int N);
inline float cal_pml_dmax(float L, float Vp, float Rpp);
inline float cal_pml_amax(float fc);
inline float cal_pml_d(float x, float L, float dmax);
inline float cal_pml_a(float x, float L, float amax);
inline float cal_pml_b(float x, float L, float bmax);
void abs_init(float *pml);
void Alloc_aux(struct aux *Aux, int isx1, int isx2, int isy1, int isy2, int isz1, int isz2, int this_rank);
float *Alloc_pml();
int coef_surface(float *M, float *matVx2Vz, float *matVy2Vz, float *pml);
#endif
