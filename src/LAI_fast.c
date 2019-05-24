#include "abc_000_macro.h"
#if defined(__GNUC__)||defined(__CLANG__) 
DISABLE_WARNING(unused-variable,unused-variable,NOT_USED)
DISABLE_WARNING(unused-function,unused-function,NOT_USED)
#endif
#if MSVC==1
#include <windows.h>               
#include "intrin.h"                
#endif
#include <inttypes.h>
#include <stdio.h>	               
#include <string.h>	               
#include <time.h>
#include <math.h>
#include <float.h>
#include "abc_001_config.h"
#include "abc_mem.h"              
#include "abc_common.h"           
#if M_INTERFACE==1
#include "hemiphoto_lib_matlab.h"
#elif R_INTERFACE==1
#include "hemiphoto_lib_R.h"
#endif
#if MYRAND_LIBRARY==1
#include "abc_rand_pcg.h"
#elif MKLRAND_LIBRARY==1
#include "abc_rand_mkl.h"
VSLStreamStatePtr stream;  
#endif
static uint64_t elapsedTime;
static int32_t arrayMax(int32_t * arr,int N)
{
	rINT max=-9999;
	for (rINT i=0; i < N; i++)
	{
		if (arr[i] > max) max=arr[i];
	}	
	return max;
}
#if MSVC==1
static LARGE_INTEGER t1,t2;
static LARGE_INTEGER Frequency;
#else
static struct timespec t1,t2;
#endif
#include "abc_optimize.h"
#include "gfun.h"
#include "abc_gaussquad.h"
FLOATPTR COS;
FLOAT    DX;
#if M_INTERFACE==1
#define pY    prhs[0]
#define pOpt  prhs[2]
#define NULL_RET NULL
void DllExport LAI_estimate(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[])
#else 
extern int nlhs;
#define NULL_RET R_NilValue
SEXP DllExport LAI_estimate(SEXP pY1,SEXP pY2,SEXP pOpt)
#endif
{
#if MKL_LIBRARY==1
#endif
#if MSVC==1
#else	
#endif
	FLOAT nan;
	nan=1.0e300;
	nan=nan*nan*0.f;
	MemPointers MEM;
	MEM=(MemPointers) {.init=mem_init};
	MEM.init(&MEM);
	OPT	opt;
#if  R_INTERFACE==1	
	get_options(&opt,pY1,pY2,pOpt,&MEM);
#elif M_INTERFACE==1
	get_options(&opt,nrhs,prhs,&MEM);
#endif
	GFUNC gfunList[]={
		gfunc_uni,
		gfunc_sph,
		gfunc_hor,
		gfunc_vtc,
		gfunc_erc,
		gfunc_pln,
		gfunc_plg,
		gfunc_ext,
		gfunc_bet,
		gfunc_elt,
		gfunc_rgm,
		gfunc_dic,
		gfunc_lin,
		gfunc_es1,
		gfunc_es2,
		gfunc_res,
		gfunc_sts,
		gfunc_lan,
		gfunc_vhf
	};
	char *gfunName[]={
		"uni",
		"sph",
		"hor",
		"vtc",
		"erc",
		"pln",
		"plg",
		"ext",
		"bet",
		"elt",
		"rgm",
		"dic",
		"lin",
		"es1",
		"es2",
		"res",
		"sts",
		"lan",
		"vhf",
	};
	int32_t NumOfGfunc=sizeof(gfunList)/sizeof(gfunList[0]);
#if M_INTERFACE==1
#define ANS plhs[0]
#else
	SEXP    ANS=NULL;
#endif
	RES res={ .LAI_HA57=NULL,};
	int32_t  nprt=0;
	if (nlhs==1){
#if R_INTERFACE==1
		ANS=allocate_output(&res,NumOfGfunc,&nprt,gfunName);
#elif M_INTERFACE==1
		ANS=allocate_output(&res,NumOfGfunc,&nprt,gfunName);
#endif
	}
	int32_t Npar_max;
	Npar_max=1L+2L ; 
	PAR_STRUCT par;
	par.THETA=opt.THETA;
	par.GAP=opt.GAP;
	par.LEN=opt.NPIXEL;
	par.dx=0.00001;
	par.boundL=MyALLOC(MEM,2 * Npar_max,FLOAT,0);
	par.boundL[0]=0;
	par.boundL[1]=15;
	par.boundG=par.boundL+2 * 1L;	
	par.Nfrac=opt.Nfrac;
	{
		rFLOAT fraction;
		rFLOAT LAI;
		rFLOAT hAngle=57.5/180 * PI;
		rINT N=par.LEN;
		rINT n1=0;
		rINT n01=0;		
		for (rINT i=0; i < N; i++)
		{ 
			rFLOAT t=par.THETA[i];
			if (t !=t) continue;
			if (fabs(t - hAngle) < 10./180 * PI)
			{
				n01++;
				if (par.GAP[i]==1) n1++;				
			}
		}
		fraction=(FLOAT)n1/n01;
		LAI=-cos(hAngle)/0.5*log(fraction);
		for (rINT i=0; i < NumOfGfunc; i++)
		{
			res.LAI_HA57[i]=LAI;
		}
	}
	FLOATPTR memChunk;
	{
		int memSizeArr[5];
		memSizeArr[0]=(opt.GQ_knotNumber+1)*(opt.GQ_knotNumber+1); 
		memSizeArr[1]=opt.NBIN *5;   
		memSizeArr[2]=opt.Nfrac *3;  
		int N=5;
		memSizeArr[3]=2*N+(N+1)*(N+1)+2*N;  
		memChunk=MyALLOC(MEM,arrayMax(memSizeArr,4),FLOAT,0);
	}
	{
		rINT  N=5;
		rFLOAT LAI=nan;
		while (N >=3 && (LAI !=LAI||LAI >=FLT_MAX||LAI < -FLT_MAX))
		{
			FLOATPTR lroots=memChunk;
			FLOATPTR weight=lroots+N;
			FLOATPTR coef_MEM=weight+N;
			lege_roots(lroots,weight,coef_MEM,N);
			rFLOAT c1=(PI/2 - 0)/2;
			rFLOAT c2=(PI/2+0)/2;
			for (rINT i=0; i < N; i++)
			{
				lroots[i]=c1*lroots[i]+c2;
			}
			FLOATPTR bndValues=weight+N;
			for (rINT i=0; i < (N - 1); i++)
			{
				bndValues[i]=(lroots[i]+lroots[i+1])/2;
			}
			bndValues[N - 1]=PI/2;
			FLOATPTR Ngap=bndValues+N;
			FLOATPTR Ntotal=Ngap+N;
			memset(Ngap,0,sizeof(FLOAT)* 2 * N);
			for (rINT i=0; i <par.LEN; i++)
			{
				rFLOAT t=par.THETA[i];				
				if (t !=t) continue;
				rINT idx=N;
				for (rINT j=1; j <=N; j++)
				{
					if (t < bndValues[j - 1])
					{
						idx=j;
						break;
					}
				}
				idx--;
				++Ntotal[idx];
				if (par.GAP[i]==1)++Ngap[idx];
			}
			LAI=0;
			for (rINT i=0; i <N; i++)
			{
				rFLOAT x=lroots[i];
				LAI+=log(Ngap[i]/Ntotal[i])*cos(x)*sin(x)*weight[i];
			}
			LAI=c1*LAI*-2;
			N--;
		}
		for (rINT i=0; i < NumOfGfunc; i++)
		{
			res.LAI_MILLER[i]=LAI;
		}
	}
	{
		NQ=opt.GQ_knotNumber;
		lroots=MyALLOC(MEM,2 * NQ,FLOAT,0);
		weight=lroots+NQ; 
		lege_roots(lroots,weight,memChunk,NQ);
	}
	int32_t  ITE=opt.ite;
	int32_t  NBIN=opt.NBIN;
	par.NBIN=NBIN;
	par.THETA01=memChunk;
	par.GAP1=memChunk+NBIN;
	par.GAP0=memChunk+2*NBIN;
	memset(par.GAP1,0,sizeof(FLOAT)* 2 * NBIN);
	rFLOAT bin=(PI/2)/NBIN;
	for (rINT i=0; i < NBIN; i++)
	{ 
		par.THETA01[i]=i*bin+bin*0.5;	 
	}
	for (rINT i=0; i < par.LEN; i++)
	{
		rFLOAT t=par.THETA[i];
		if (t !=t) continue;
		rINT idx=(t/bin);
		idx=idx >=NBIN ? NBIN - 1 : idx;
		if (par.GAP[i]==1)
			par.GAP1[idx]++;
		else
			par.GAP0[idx]++;	
	}
	par.mem=MyALLOC(MEM,(2 * Npar_max+2 * par.LEN),FLOAT,0);
	FLOATPTR X=MyALLOC(MEM,2 * Npar_max,FLOAT,0);
	FLOATPTR df=X+Npar_max;	
	FLOATPTR mem4Mininize=MyALLOC(MEM,5 * Npar_max,FLOAT,0);
	int32_t  Npar;
	for (rINT i=0; i <NumOfGfunc; i++)
	{
		{
			selectGfunc(gfunName[i],&par);
			Npar=1L+par.Npar;
			X[0]=2L;
			memcpy(X+1L,par.Gparam_init,sizeof(FLOAT)*par.Npar);
			logitfwd(X,X,par.boundL,Npar);
		}
		miminize(alg_BNR,X,Npar,&par,ITE,mem4Mininize);
		FLOAT f,AIC;
		f=alg_BNR(X,df,Npar,&par);
		AIC=2 * f+2 * (1+par.Npar);
		logitinv(X,X,par.boundL,Npar);		
		res.LAI_BNR[i]=X[0];
		res.AIC_BNR[i]=AIC;
		if (par.Npar==0)
			res.LAD1_BNR[i]=nan,res.LAD2_BNR[i]=nan;
		else if (par.Npar==1)
			res.LAD1_BNR[i]=X[1],res.LAD2_BNR[i]=nan;
		else
			res.LAD1_BNR[i]=X[1],res.LAD2_BNR[i]=X[2];
		r_printf("%s LAI:%.4f AIC:%.4f\n",gfunName[i],X[0],AIC) ;
	}
	NBIN=opt.Nfrac;
	par.THETAfrac=memChunk;
	par.GAP1=memChunk+NBIN;
	par.GAP0=memChunk+2*NBIN;
	par.FRAC=memChunk+3*NBIN;	
	memset(par.GAP1,0,sizeof(FLOAT)*3*NBIN);
	bin=(PI/2)/NBIN;
	for (rINT i=0; i < NBIN; i++)
		par.THETAfrac[i]=i*bin+bin*0.5;
	for (rINT i=0; i < par.LEN; i++)
	{
		rFLOAT t=par.THETA[i];
		rINT idx=(t/bin);
		idx=idx >=NBIN ? NBIN - 1 : idx;
		if (par.GAP[i]==1)
			par.GAP1[idx]++;
		else if (par.GAP[i]==0)
			par.GAP0[idx]++; 
	}
	int32_t newNbin=0;
	for (rINT i=0; i < NBIN; i++)
	{
		rFLOAT sum=(par.GAP1[i]+par.GAP0[i]);
		if (sum>0.1)
		{
			par.FRAC[newNbin]=par.GAP1[i]/sum;
			par.THETAfrac[newNbin]=par.THETAfrac[i];
			newNbin++;
		}		
	}
	par.Nfrac=newNbin;
	for (rINT i=0; i <NumOfGfunc; i++)
	{
		{
			selectGfunc(gfunName[i],&par);
			Npar=1L+par.Npar;
			X[0]=2L;
			memcpy(X+1L,par.Gparam_init,sizeof(FLOAT)*par.Npar);
			logitfwd(X,X,par.boundL,Npar);
		}
		miminize(alg_L2F,X,Npar,&par,ITE,mem4Mininize);
		FLOAT f,AIC;
		f=alg_L2F(X,df,Npar,&par);
		AIC=par.Nfrac*log( 2 * f)+2 * (1+par.Npar);
		logitinv(X,X,par.boundL,Npar);
		res.LAI_L2F[i]=X[0];
		res.AIC_L2F[i]=AIC;
		if (par.Npar==0)
			res.LAD1_L2F[i]=nan,res.LAD2_L2F[i]=nan;
		else if (par.Npar==1)
			res.LAD1_L2F[i]=X[1],res.LAD2_L2F[i]=nan;
		else
			res.LAD1_L2F[i]=X[1],res.LAD2_L2F[i]=X[2];
		r_printf("%s LAI:%.4f AIC:%.4f\n",gfunName[i],X[0],AIC);
	}
	NBIN=opt.Nfrac;
	par.THETAfrac=memChunk;
	par.GAP1=memChunk+NBIN;
	par.GAP0=memChunk+2 * NBIN;
	par.FRAC=memChunk+3 * NBIN;
	par.WEIGHT=memChunk+4 * NBIN;
	par.Nfrac=newNbin;
	for (rINT i=0; i < par.Nfrac; i++)
		par.WEIGHT[i]=1.;
	for (rINT i=0; i <NumOfGfunc; i++)
	{
		{
			selectGfunc(gfunName[i],&par);
			Npar=1L+par.Npar;
			X[0]=2L;
			memcpy(X+1L,par.Gparam_init,sizeof(FLOAT)*par.Npar);
			logitfwd(X,X,par.boundL,Npar);
		}
		for (rINT I=0; I < 10; I++)
		{
			miminize(alg_L1F,X,Npar,&par,ITE,mem4Mininize);
			alg_L1F_weight(X,Npar,&par);
		}		
		FLOAT f,AIC;
		f=alg_L1F(X,df,Npar,&par);
		AIC=par.Nfrac*log(2 * f)+2 * (1+par.Npar);
		logitinv(X,X,par.boundL,Npar);
		res.LAI_L1F[i]=X[0];
		res.AIC_L1F[i]=AIC;
		if (par.Npar==0)
			res.LAD1_L1F[i]=nan,res.LAD2_L1F[i]=nan;
		else if (par.Npar==1)
			res.LAD1_L1F[i]=X[1],res.LAD2_L1F[i]=nan;
		else
			res.LAD1_L1F[i]=X[1],res.LAD2_L1F[i]=X[2];
		r_printf("%s LAI:%.4f AIC:%.4f\n",gfunName[i],X[0],AIC);
	}
	NBIN=opt.Nfrac;
	par.Nfrac=NBIN;
	bin=(PI/2)/NBIN;
	for (rINT i=0; i < NBIN; i++)
		par.THETAfrac[i]=i*bin+bin*0.5;
	 newNbin=0;
	for (rINT i=0; i < NBIN; i++)
	{
		rFLOAT sum=(par.GAP1[i]+par.GAP0[i]);
		if (sum > 1 && par.GAP1[i] > 0)
		{
			par.FRAC[newNbin]=log(par.GAP1[i]/sum);
			par.THETAfrac[newNbin]=par.THETAfrac[i];
			newNbin++;
		}
	}
	par.Nfrac=newNbin;
	for (rINT i=0; i <NumOfGfunc; i++)
	{
		{
			selectGfunc(gfunName[i],&par);
			Npar=1L+par.Npar;
			X[0]=2L;
			memcpy(X+1L,par.Gparam_init,sizeof(FLOAT)*par.Npar);
			logitfwd(X,X,par.boundL,Npar);
		}
		miminize(alg_L2LF,X,Npar,&par,ITE,mem4Mininize);
		FLOAT f,AIC;
		f=alg_L2LF(X,df,Npar,&par);
		AIC=par.Nfrac*log(2 * f)+2 * (1+par.Npar);
		logitinv(X,X,par.boundL,Npar);
		res.LAI_L2LF[i]=X[0];
		res.AIC_L2LF[i]=AIC;
		if (par.Npar==0)
			res.LAD1_L2LF[i]=nan,res.LAD2_L2LF[i]=nan;
		else if (par.Npar==1)
			res.LAD1_L2LF[i]=X[1],res.LAD2_L2LF[i]=nan;
		else
			res.LAD1_L2LF[i]=X[1],res.LAD2_L2LF[i]=X[2];
		r_printf("%s LAI:%.4f AIC:%.4f\n",gfunName[i],X[0],AIC);
	}
	NBIN=opt.Nfrac;
	par.THETAfrac=memChunk;
	par.GAP1=memChunk+NBIN;
	par.GAP0=memChunk+2 * NBIN;
	par.FRAC=memChunk+3 * NBIN;
	par.WEIGHT=memChunk+4 * NBIN;
	par.Nfrac=newNbin;
	for (rINT i=0; i < par.Nfrac; i++)
		par.WEIGHT[i]=1.;
	for (rINT i=0; i <NumOfGfunc; i++)
	{
		{
			selectGfunc(gfunName[i],&par);
			Npar=1L+par.Npar;
			X[0]=2L;
			memcpy(X+1L,par.Gparam_init,sizeof(FLOAT)*par.Npar);
			logitfwd(X,X,par.boundL,Npar);
		}
		for (rINT I=0; I < 10; I++)
		{
			miminize(alg_L1LF,X,Npar,&par,ITE,mem4Mininize);
			alg_L1LF_weight(X,Npar,&par);
		}
		FLOAT f,AIC;
		f=alg_L1LF(X,df,Npar,&par);
		AIC=par.Nfrac*log(2 * f)+2 * (1+par.Npar);
		logitinv(X,X,par.boundL,Npar);
		res.LAI_L1LF[i]=X[0];
		res.AIC_L1LF[i]=AIC;
		if (par.Npar==0)
			res.LAD1_L1LF[i]=nan,res.LAD2_L1LF[i]=nan;
		else if (par.Npar==1)
			res.LAD1_L1LF[i]=X[1],res.LAD2_L1LF[i]=nan;
		else
			res.LAD1_L1LF[i]=X[1],res.LAD2_L1LF[i]=X[2];
		r_printf("%s LAI:%.4f AIC:%.4f\n",gfunName[i],X[0],AIC);
	}
	NBIN=opt.Nfrac;
	par.Nfrac=NBIN;
	memset(par.GAP1,0,sizeof(FLOAT)* 3 * NBIN);
	bin=(PI/2)/NBIN;
	for (rINT i=0; i < NBIN; i++)
		par.THETAfrac[i]=i*bin+bin*0.5;
	for (rINT i=0; i < par.LEN; i++)
	{
		rFLOAT t=par.THETA[i];
		rINT idx=(t/bin);
		idx=idx >=NBIN ? NBIN - 1 : idx;
		if (par.GAP[i]==1)
			par.GAP1[idx]++;
		else if (par.GAP[i]==0)
			par.GAP0[idx]++;
	}
	newNbin=0;
	for (rINT i=0; i < NBIN; i++)
	{
		rFLOAT sum=(par.GAP1[i]+par.GAP0[i]);
		if (sum > 1 && par.GAP1[i] > 0)
		{
			par.THETAfrac[newNbin]=par.THETAfrac[i];
			par.FRAC[newNbin]=cos(par.THETAfrac[newNbin])*log(par.GAP1[i]/sum);
			newNbin++;
		}
	}
	par.Nfrac=newNbin;
	for (rINT i=0; i <NumOfGfunc; i++)
	{
		{
			selectGfunc(gfunName[i],&par);
			Npar=1L+par.Npar;
			X[0]=2L;
			memcpy(X+1L,par.Gparam_init,sizeof(FLOAT)*par.Npar);
			logitfwd(X,X,par.boundL,Npar);
		}
		miminize(alg_L2WL,X,Npar,&par,ITE,mem4Mininize);
		FLOAT f,AIC;
		f=alg_L2WL(X,df,Npar,&par);
		AIC=par.Nfrac*log(2 * f)+2 * (1+par.Npar);
		logitinv(X,X,par.boundL,Npar);
		res.LAI_L2WL[i]=X[0];
		res.AIC_L2WL[i]=AIC;
		if (par.Npar==0)
			res.LAD1_L2WL[i]=nan,res.LAD2_L2WL[i]=nan;
		else if (par.Npar==1)
			res.LAD1_L2WL[i]=X[1],res.LAD2_L2WL[i]=nan;
		else
			res.LAD1_L2WL[i]=X[1],res.LAD2_L2WL[i]=X[2];
		r_printf("%s LAI:%.4f AIC:%.4f\n",gfunName[i],X[0],AIC);
	}
	NBIN=opt.Nfrac;
	par.THETAfrac=memChunk;
	par.GAP1=memChunk+NBIN;
	par.GAP0=memChunk+2 * NBIN;
	par.FRAC=memChunk+3 * NBIN;
	par.WEIGHT=memChunk+4 * NBIN;
	par.Nfrac=newNbin;
	for (rINT i=0; i < par.Nfrac; i++)
		par.WEIGHT[i]=1.;
	for (rINT i=0; i <NumOfGfunc; i++)
	{
		{
			selectGfunc(gfunName[i],&par);
			Npar=1L+par.Npar;
			X[0]=2L;
			memcpy(X+1L,par.Gparam_init,sizeof(FLOAT)*par.Npar);
			logitfwd(X,X,par.boundL,Npar);
		}
		for (rINT I=0; I < 10; I++)
		{
			miminize(alg_L1WL,X,Npar,&par,ITE,mem4Mininize);
			alg_L1WL_weight(X,Npar,&par);
		}
		FLOAT f,AIC;
		f=alg_L1WL(X,df,Npar,&par);
		AIC=par.Nfrac*log(2 * f)+2 * (1+par.Npar);
		logitinv(X,X,par.boundL,Npar);
		res.LAI_L1WL[i]=X[0];
		res.AIC_L1WL[i]=AIC;
		if (par.Npar==0)
			res.LAD1_L1WL[i]=nan,res.LAD2_L1WL[i]=nan;
		else if (par.Npar==1)
			res.LAD1_L1WL[i]=X[1],res.LAD2_L1WL[i]=nan;
		else
			res.LAD1_L1WL[i]=X[1],res.LAD2_L1WL[i]=X[2];
		r_printf("%s LAI:%.4f AIC:%.4f\n",gfunName[i],X[0],AIC);
	}
	rFLOAT minAIC;
	FLOATPTR aicArr;
	aicArr=res.AIC_BNR;
	minAIC=1e200;
	for (rINT i=0; i < NumOfGfunc; i++)
	{
		if (aicArr[i] < minAIC)	minAIC=aicArr[i];
	}
	for (rINT i=0; i < NumOfGfunc; i++) aicArr[i] -=minAIC;
	aicArr=res.AIC_L2F;
	minAIC=1e200;
	for (rINT i=0; i < NumOfGfunc; i++)
	{
		if (aicArr[i] < minAIC)	minAIC=aicArr[i];
	}
	for (rINT i=0; i < NumOfGfunc; i++) aicArr[i] -=minAIC;
	aicArr=res.AIC_L1F;
	minAIC=1e200;
	for (rINT i=0; i < NumOfGfunc; i++)
	{
		if (aicArr[i] < minAIC)	minAIC=aicArr[i];
	}
	for (rINT i=0; i < NumOfGfunc; i++) aicArr[i] -=minAIC;
	aicArr=res.AIC_L2LF;
	minAIC=1e200;
	for (rINT i=0; i < NumOfGfunc; i++)
	{
		if (aicArr[i] < minAIC)	minAIC=aicArr[i];
	}
	for (rINT i=0; i < NumOfGfunc; i++) aicArr[i] -=minAIC;
	aicArr=res.AIC_L1LF;
	minAIC=1e200;
	for (rINT i=0; i < NumOfGfunc; i++)
	{
		if (aicArr[i] < minAIC)	minAIC=aicArr[i];
	}
	for (rINT i=0; i < NumOfGfunc; i++) aicArr[i] -=minAIC;
	aicArr=res.AIC_L2WL;
	minAIC=1e200;
	for (rINT i=0; i < NumOfGfunc; i++)
	{
		if (aicArr[i] < minAIC)	minAIC=aicArr[i];
	}
	for (rINT i=0; i < NumOfGfunc; i++) aicArr[i] -=minAIC;
	aicArr=res.AIC_L1WL;
	minAIC=1e200;
	for (rINT i=0; i < NumOfGfunc; i++)
	{
		if (aicArr[i] < minAIC)	minAIC=aicArr[i];
	}
	for (rINT i=0; i < NumOfGfunc; i++) aicArr[i] -=minAIC;
#if R_INTERFACE==1	
	UNPROTECT(nprt);	
	MEM.free_all(&MEM);
	return ANS;
#endif
	MEM.free_all(&MEM);
	return NULL_RET;
} 
#if defined(__GNUC__)||defined(__CLANG__) 
ENABLE_WARNING(unused-variable,unused-variable,NOT_USED)
ENABLE_WARNING(unused-function,unused-function,NOT_USED)
#endif
