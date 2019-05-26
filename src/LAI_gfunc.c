#include "abc_000_macro.h"
#if defined(__GNUC__)||defined(__CLANG__) 
DISABLE_WARNING(unused - variable,unused - variable,NOT_USED)
DISABLE_WARNING(unused - function,unused - function,NOT_USED)
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
#include "abc_gaussquad.h"
#include "gfun.h"
#include "abc_rand_pcg.h"
FLOATPTR COS;
FLOAT    DX;
#if M_INTERFACE==1
#define pY     prhs[0]
#define pOpt   prhs[2]
#define NULL_RET NULL
void DllExport LAI_gfunc(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[])
#else 
extern int nlhs;
#define NULL_RET R_NilValue
SEXP DllExport LAI_gfunc(SEXP pY1,SEXP pY2,SEXP pOpt)
#endif
{
#if MSVC==1
#else	
#endif
	FLOAT nan;
	nan=1.0e300;
	nan=nan*nan*0.f;
	MemPointers MEM;
	MEM=(MemPointers) { .init=mem_init };
	MEM.init(&MEM);
	OPT	opt;
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
	SEXP    ansListNames;
#endif
	PAR_STRUCT par;
	int Npar_max=3;
	par.boundL=MyALLOC(MEM,2 * Npar_max,FLOAT,0);
	par.boundL[0]=0;
	par.boundL[1]=15;
	par.boundG=par.boundL+2 * 1L;	
	FLOATPTR THETA;
	FLOATPTR G;
	FLOATPTR PAR;
	int32_t  N=200;
	int32_t  idxGfunc;	
	FLOATPTR inPar;	
	int32_t  isInputRandomNumber;
#if M_INTERFACE==1
	mxArray *tmp;
	if (mxIsChar(prhs[0]))
	{
		char *tmpStr=mxArrayToString(prhs[0]);
		idxGfunc=gfuncSelect(tmpStr,gfunName,NumOfGfunc)+1;
		mxFree(tmpStr);
	}
	else
	{
		idxGfunc=(int32_t)mxGetScalar(prhs[0]);
	}
	inPar=(FLOATPTR) mxGetData(  prhs[1] ) ;
	tmp=mxGetField(prhs[2],0,"scaledParameter");
	if (tmp==NULL)
		isInputRandomNumber=1;
	else
		isInputRandomNumber=(int32_t)mxGetScalar(tmp);
	size_t dim[2]={1,1};
	char *fldNames[]={"THETA","G","PAR"};
	plhs[0]=mxCreateStructArray(2L,dim,3,fldNames);
	tmp=mxCreateNumericMatrix(N,1,mxDOUBLE_CLASS,mxREAL);
	mxSetField(plhs[0],0,"THETA",tmp);
	THETA=(FLOATPTR)mxGetData(tmp);
	tmp=mxCreateNumericMatrix(N,1,mxDOUBLE_CLASS,mxREAL);
	mxSetField(plhs[0],0,"G",tmp);
	G=(FLOATPTR)mxGetData(tmp);
	tmp=mxCreateNumericMatrix(2,1,mxDOUBLE_CLASS,mxREAL);
	mxSetField(plhs[0],0,"PAR",tmp);
	PAR=(FLOATPTR)mxGetData(tmp);
#elif  R_INTERFACE==1
	SEXP tmpSEXP;
	int32_t nprt=0;
	if (isString(pY1))
	{
		char *tmpStr=CHAR(STRING_ELT(pY1,0));		
		idxGfunc=gfuncSelect(tmpStr,gfunName,NumOfGfunc)+1;
	}
	else if (isReal(pY1))
	{
		idxGfunc=(int32_t) asReal(pY1);
	}
	else if (isInteger(pY1))
	{
		idxGfunc=asInteger(pY1);
	}
	inPar=(FLOATPTR)REAL(pY2);
	if (!isNewList(pOpt))
	{
		MEM.free_all(&MEM);
		r_error("Error: The input parameter OPTION should be a List variable!\n");
	}
	tmpSEXP=getListElement(pOpt,"scaledParameter");
	isInputRandomNumber=(tmpSEXP !=NULL) ? (PROTECT(tmpSEXP=coerceVector(tmpSEXP,REALSXP)),nprt++,(int)asReal(tmpSEXP)) : 1;
	size_t dim[2]={ 1,1 };
	char *fldNames[]={ "THETA","G","PAR" };
	int32_t elementNum=sizeof(fldNames)/sizeof(fldNames[0]);
	PROTECT(ANS=allocVector(VECSXP,elementNum));++nprt;
	PROTECT(ansListNames=allocVector(STRSXP,elementNum));++nprt;	  
	int32_t idx;
	idx=0;
	PROTECT(tmpSEXP=allocVector(REALSXP,N));++nprt;
	THETA=(FLOATPTR)REAL(tmpSEXP); SET_VECTOR_ELT(ANS,idx,tmpSEXP); SET_STRING_ELT(ansListNames,idx,mkChar(fldNames[idx]));
	idx++;
	PROTECT(tmpSEXP=allocVector(REALSXP,N));++nprt;
	G=(FLOATPTR)REAL(tmpSEXP); SET_VECTOR_ELT(ANS,idx,tmpSEXP); SET_STRING_ELT(ansListNames,idx,mkChar(fldNames[idx]));
	idx++;
	PROTECT(tmpSEXP=allocVector(REALSXP,2));++nprt;
	PAR=(FLOATPTR)REAL(tmpSEXP); SET_VECTOR_ELT(ANS,idx,tmpSEXP); SET_STRING_ELT(ansListNames,idx,mkChar(fldNames[idx]));
	setAttrib(ANS,R_NamesSymbol,ansListNames);
#endif
	selectGfunc(gfunName[idxGfunc - 1],&par);
	r_printf("%s\n",gfunName[idxGfunc - 1]);
	if (isInputRandomNumber==0)
	{
		if (par.Npar==0)
			PAR[0]=PAR[1]=nan;
		else if (par.Npar==1)
			PAR[0]=inPar[0],
			PAR[1]=nan;
		else if (par.Npar==2)
			PAR[0]=inPar[0],
			PAR[1]=inPar[1];
	}
	else
	{
		if (par.Npar==0)
			PAR[0]=PAR[1]=nan;
		else if (par.Npar==1)
			PAR[0]=(par.boundG[1] - par.boundG[0])*inPar[0]+par.boundG[0],
			PAR[1]=nan;
		else if (par.Npar==2)
			PAR[0]=(par.boundG[1] - par.boundG[0])*inPar[0]+par.boundG[0],
			PAR[1]=(par.boundG[3] - par.boundG[2])*inPar[1]+par.boundG[2];
	}
	FLOAT   bin=PI/2/N;
	for (rINT i=0; i < N; i++)
	{
		THETA[i]=(i+0.5)*bin;
	}
	{
		NQ=21;
		FLOATPTR  memChunk=MyALLOC(MEM,(NQ+1)*(NQ+1),FLOAT,0);
		lroots=MyALLOC(MEM,2 * NQ,FLOAT,0);
		weight=lroots+NQ;
		lege_roots(lroots,weight,memChunk,NQ);
	}
	par.gfun(THETA,G,N,PAR,NULL);
#if R_INTERFACE==1	
	UNPROTECT(nprt);
	MEM.free_all(&MEM);
	return ANS;
#endif
	MEM.free_all(&MEM);
	return NULL_RET;
} 
#if defined(__GNUC__)||defined(__CLANG__) 
ENABLE_WARNING(unused - variable,unused - variable,NOT_USED)
ENABLE_WARNING(unused - function,unused - function,NOT_USED)
#endif
