#pragma once
#include "abc_001_config.h"
#include <stdio.h>
#include "abc_mem.h"
#include "math.h"
#include "abc_datatype.h"
typedef void(*GFUNC)(FLOATPTR theta,FLOATPTR G,int32_t N,FLOATPTR gParam,FLOATPTR mem);
typedef struct
{
	FLOATPTR boundL;
	FLOATPTR boundG;
	GFUNC    gfun;
	FLOATPTR THETA;
	FLOATPTR GAP;
	FLOATPTR THETA01;
	FLOATPTR GAP1;
	FLOATPTR GAP0;
	FLOATPTR THETAfrac;
	FLOATPTR FRAC;
	int32_t  Nfrac;
	FLOATPTR WEIGHT;
	int32_t LEN;
	FLOATPTR mem;
	int32_t NBIN;
	FLOAT   dx;
	int32_t Npar;
	FLOAT   Gparam_init[2];
} PAR_STRUCT;
typedef struct OPT
{
	FLOATPTR THETA;
	FLOATPTR GAP;
	int32_t  NPIXEL;
	int32_t  NBIN;
	int32_t  ite;
	int32_t  Nfrac;
	uint8_t  algorithm;
	int32_t  GQ_knotNumber;
} OPT;
typedef struct RES
{
	void * LAD_TYPE;
	FLOATPTR LAI_HA57;
	FLOATPTR LAI_MILLER;
	FLOATPTR LAI_BNR;
	FLOATPTR AIC_BNR;
	FLOATPTR LAD1_BNR;
	FLOATPTR LAD2_BNR;
	FLOATPTR LAI_L2F;
	FLOATPTR AIC_L2F;
	FLOATPTR LAD1_L2F;
	FLOATPTR LAD2_L2F;
	FLOATPTR LAI_L1F;
	FLOATPTR AIC_L1F;
	FLOATPTR LAD1_L1F;
	FLOATPTR LAD2_L1F;
	FLOATPTR LAI_L2LF;
	FLOATPTR AIC_L2LF;
	FLOATPTR LAD1_L2LF;
	FLOATPTR LAD2_L2LF;
	FLOATPTR LAI_L1LF;
	FLOATPTR AIC_L1LF;
	FLOATPTR LAD1_L1LF;
	FLOATPTR LAD2_L1LF;
	FLOATPTR LAI_L2WL;
	FLOATPTR AIC_L2WL;
	FLOATPTR LAD1_L2WL;
	FLOATPTR LAD2_L2WL;
	FLOATPTR LAI_L1WL;
	FLOATPTR AIC_L1WL;
	FLOATPTR LAD1_L1WL;
	FLOATPTR LAD2_L1WL; 
} RES;
extern void logitinv(FLOATPTR x,FLOATPTR y,FLOATPTR bnd,int32_t N);
extern void logitfwd(FLOATPTR x,FLOATPTR y,FLOATPTR bnd,int32_t N);
extern FLOAT log_prob_beerlaw(FLOATPTR X,FLOATPTR df,int32_t N,void *par);
extern FLOAT alg_BNR(FLOATPTR X,FLOATPTR df,int32_t N,void *par);
extern FLOAT alg_L2F(FLOATPTR X,FLOATPTR df,int32_t N,void *par);
extern FLOAT alg_L1F(FLOATPTR X,FLOATPTR df,int32_t N,void *par);
extern FLOAT alg_L1F_weight(FLOATPTR X,int32_t N,void *par);
extern FLOAT alg_L2LF(FLOATPTR X,FLOATPTR df,int32_t N,void *par);
extern FLOAT alg_L1LF(FLOATPTR X,FLOATPTR df,int32_t N,void *par);
extern FLOAT alg_L1LF_weight(FLOATPTR X,int32_t N,void *par);
extern FLOAT alg_L2WL(FLOATPTR X,FLOATPTR df,int32_t N,void *par);
extern FLOAT alg_L1WL(FLOATPTR X,FLOATPTR df,int32_t N,void *par);
extern FLOAT alg_L1WL_weight(FLOATPTR X,int32_t N,void *par);
