#pragma once
#include <stdio.h>
#include "abc_001_config.h"
#include "abc_blas_lapack_lib.h"
#include "abc_mem.h"
#include "abc_datatype.h"
#if M_INTERFACE==1
	#define  r_printf(...) mexPrintf(__VA_ARGS__)
	#define  r_error(x)    mexErrMsgTxt(x)
	#define  r_malloc(x)   mxMalloc(x) 
	#define  r_free(x)     mxFree(x)
#elif R_INTERFACE==1
	#define  r_printf(...)  Rprintf(__VA_ARGS__)
	#define  r_error(x)     error(x)
	#define  r_malloc(x)    Calloc(x,char)  
	#define  r_free(x)      Free(x) 
#endif
extern void    fill_1_to_N(rFLTPTR p,int N);
extern void    quickSortD(float * _restrict arr,int32_t * _restrict index,int32_t low,int32_t high);
extern void    quickSortA(float * _restrict arr,int32_t * _restrict index,int32_t low,int32_t high);
extern void    transpose_inplace(float * _restrict m,int w,int h);
extern int32_t strcicmp(char const * _restrict a,char const * _restrict b);
extern float   determine_period(FLTPTR data,int32_t N,float  omissionValue);
extern int32_t find_changepoint(FLTPTR prob,FLTPTR mem,float threshold,INTPTR cpt,FLTPTR cptCI,int32_t N,int32_t minSepDist,int32_t maxCptNumber);
extern float fastlog(float x);
extern float sum_log_diag(rFLTPTR p,rINT K);
extern float fastexp(float x);
extern float fast_sqrt(float x);
static INLINE void normalize(rFLTPTR ptr,int N)
{
	float mean,std;
	r_ippsMeanStdDev_32f(ptr,N,&mean,&std,ippAlgHintAccurate);
	r_ippsSubC_32f_I(mean,ptr,N);
	r_cblas_sscal(N,1.f/std,ptr,1L);
}
static INLINE void normalize_x_factor(float * _restrict ptr,int N,float factor)
{
	float mean,std;
	r_ippsMeanStdDev_32f(ptr,N,&mean,&std,ippAlgHintAccurate);
	r_ippsSubC_32f_I(mean,ptr,N);
	r_cblas_sscal(N,1/std*factor,ptr,1);
}
