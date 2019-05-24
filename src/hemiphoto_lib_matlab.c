#include "abc_000_macro.h"
static char fileID='c';
#if M_INTERFACE==1
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <float.h>
#include "abc_001_config.h"
#include "abc_blas_lapack_lib.h"
#include "abc_common.h"
#include "hemiphoto_lib_matlab.h"
#include "hemiphoto.h"
bool  get_options(OPT * _restrict opt,int32_t nrhs,const mxArray * _restrict prhs[],MemPointers * _restrict MEM)
{
	if (nrhs < 2)
	return false;
	opt->THETA=mxGetData(prhs[0]);
	opt->GAP=mxGetData(prhs[1]);
	opt->NPIXEL=mxGetM(prhs[0])* mxGetN(prhs[0]);
	void *ptr=NULL;
	if (nrhs < 3)
	{
		opt->ite=-100;
		opt->NBIN=200;
		opt->Nfrac=8;
		opt->GQ_knotNumber=21;
		opt->algorithm=1;
	}
	else if (mxIsStruct(prhs[2]))
	{
		ptr=mxGetField(prhs[2],0,"ite");
		opt->ite=(ptr !=NULL) ? mxGetScalar((mxArray*)ptr) : 50;
		opt->ite=abs(opt->ite) < 30 ? 30 : opt->ite;
		ptr=mxGetField(prhs[2],0,"nbin");
		opt->NBIN=(ptr !=NULL) ? mxGetScalar((mxArray*)ptr) : 200;
		opt->NBIN=opt->NBIN  < 100 ? 100 : opt->NBIN;
		ptr=mxGetField(prhs[2],0,"nfrac");
		opt->Nfrac=(ptr !=NULL) ? mxGetScalar((mxArray*)ptr) : 8;
		opt->Nfrac=opt->Nfrac  > 10 ? 10 : opt->Nfrac;
		ptr=mxGetField(prhs[2],0,"gq_knotNum");
		opt->GQ_knotNumber=(ptr !=NULL) ? mxGetScalar((mxArray*)ptr) : 21;
		opt->GQ_knotNumber=opt->GQ_knotNumber  < 10 ? 21 : opt->Nfrac;
	}
	else if (mxIsChar(prhs[2])) 
	{
		opt->ite=-100;
		opt->NBIN=200;
		opt->Nfrac=8;
		opt->GQ_knotNumber=21;
		ptr=(void *)mxArrayToString(prhs[2]);
		if (strcicmp(ptr,"mle")==0)
			opt->algorithm=1;
		else
			opt->algorithm=2;
		r_free(ptr);
	}
	else
	{
		opt->ite=-100;
		opt->NBIN=200;
		opt->Nfrac=8;
		opt->GQ_knotNumber=21;
		opt->algorithm=1;
	}
	return true;
}
mxArray *allocate_output(RES * _restrict result,int32_t N,int * nptr,char **gfuncName)
{
	mxArray * _restrict out;
	mxArray * _restrict mxPointer;
	mwSize	 dim2[2]={ 1,1 };	
	mwSize 	 dim3[3]={ N,2,11 };
	char *names[]={
		"LAD_TYPE","LAI_HA57","LAI_MILLER",
		"LAI_BNR","AIC_BNR","LAD1_BNR","LAD2_BNR",
		"LAI_L2F","AIC_L2F","LAD1_L2F","LAD2_L2F",
		"LAI_L1F","AIC_L1F","LAD1_L1F","LAD2_L1F",
		"LAI_L2LF","AIC_L2LF","LAD1_L2LF","LAD2_L2LF",
		"LAI_L1LF","AIC_L1LF","LAD1_L1LF","LAD2_L1LF",
		"LAI_L2WL","AIC_L2WL","LAD1_L2WL","LAD2_L2WL",
		"LAI_L1WL","AIC_L1WL","LAD1_L1WL","LAD2_L1WL"
	};
	int  fieldNumber=sizeof(names)/sizeof(names[0]);
	out=mxCreateStructArray(2,dim2,fieldNumber,names);
	int idx=0;
	mxPointer=mxCreateCharMatrixFromStrings(N,gfuncName); mxSetField(out,0,names[idx],mxPointer);
	idx++;
	mxPointer=mxCreateNumericMatrix(N,1,mxDOUBLE_CLASS,mxREAL); 	mxSetField(out,0,names[idx],mxPointer);
	result->LAI_HA57=mxGetData(mxPointer); idx++;
	mxPointer=mxCreateNumericMatrix(N,1,mxDOUBLE_CLASS,mxREAL); 	mxSetField(out,0,names[idx],mxPointer);
	result->LAI_MILLER=mxGetData(mxPointer); idx++;
	mxPointer=mxCreateNumericMatrix(N,1,mxDOUBLE_CLASS,mxREAL); 	mxSetField(out,0,names[idx],mxPointer);
	result->LAI_BNR=mxGetData(mxPointer); idx++;
	mxPointer=mxCreateNumericMatrix(N,1,mxDOUBLE_CLASS,mxREAL); 	mxSetField(out,0,names[idx],mxPointer);
	result->AIC_BNR=mxGetData(mxPointer); idx++;
	mxPointer=mxCreateNumericMatrix(N,1,mxDOUBLE_CLASS,mxREAL); 	mxSetField(out,0,names[idx],mxPointer);
	result->LAD1_BNR=mxGetData(mxPointer); idx++;
	mxPointer=mxCreateNumericMatrix(N,1,mxDOUBLE_CLASS,mxREAL); 	mxSetField(out,0,names[idx],mxPointer);
	result->LAD2_BNR=mxGetData(mxPointer); idx++;
	mxPointer=mxCreateNumericMatrix(N,1,mxDOUBLE_CLASS,mxREAL); 	mxSetField(out,0,names[idx],mxPointer);
	result->LAI_L2F=mxGetData(mxPointer); idx++;
	mxPointer=mxCreateNumericMatrix(N,1,mxDOUBLE_CLASS,mxREAL); 	mxSetField(out,0,names[idx],mxPointer);
	result->AIC_L2F=mxGetData(mxPointer); idx++;
	mxPointer=mxCreateNumericMatrix(N,1,mxDOUBLE_CLASS,mxREAL); 	mxSetField(out,0,names[idx],mxPointer);
	result->LAD1_L2F=mxGetData(mxPointer); idx++;
	mxPointer=mxCreateNumericMatrix(N,1,mxDOUBLE_CLASS,mxREAL); 	mxSetField(out,0,names[idx],mxPointer);
	result->LAD2_L2F=mxGetData(mxPointer); idx++;
	mxPointer=mxCreateNumericMatrix(N,1,mxDOUBLE_CLASS,mxREAL); 	mxSetField(out,0,names[idx],mxPointer);
	result->LAI_L1F=mxGetData(mxPointer); idx++;
	mxPointer=mxCreateNumericMatrix(N,1,mxDOUBLE_CLASS,mxREAL); 	mxSetField(out,0,names[idx],mxPointer);
	result->AIC_L1F=mxGetData(mxPointer); idx++;
	mxPointer=mxCreateNumericMatrix(N,1,mxDOUBLE_CLASS,mxREAL); 	mxSetField(out,0,names[idx],mxPointer);
	result->LAD1_L1F=mxGetData(mxPointer); idx++;
	mxPointer=mxCreateNumericMatrix(N,1,mxDOUBLE_CLASS,mxREAL); 	mxSetField(out,0,names[idx],mxPointer);
	result->LAD2_L1F=mxGetData(mxPointer); idx++;
	mxPointer=mxCreateNumericMatrix(N,1,mxDOUBLE_CLASS,mxREAL); 	mxSetField(out,0,names[idx],mxPointer);
	result->LAI_L2LF=mxGetData(mxPointer); idx++;
	mxPointer=mxCreateNumericMatrix(N,1,mxDOUBLE_CLASS,mxREAL); 	mxSetField(out,0,names[idx],mxPointer);
	result->AIC_L2LF=mxGetData(mxPointer); idx++;
	mxPointer=mxCreateNumericMatrix(N,1,mxDOUBLE_CLASS,mxREAL); 	mxSetField(out,0,names[idx],mxPointer);
	result->LAD1_L2LF=mxGetData(mxPointer); idx++;
	mxPointer=mxCreateNumericMatrix(N,1,mxDOUBLE_CLASS,mxREAL); 	mxSetField(out,0,names[idx],mxPointer);
	result->LAD2_L2LF=mxGetData(mxPointer); idx++;
	mxPointer=mxCreateNumericMatrix(N,1,mxDOUBLE_CLASS,mxREAL); 	mxSetField(out,0,names[idx],mxPointer);
	result->LAI_L1LF=mxGetData(mxPointer); idx++;
	mxPointer=mxCreateNumericMatrix(N,1,mxDOUBLE_CLASS,mxREAL); 	mxSetField(out,0,names[idx],mxPointer);
	result->AIC_L1LF=mxGetData(mxPointer); idx++;
	mxPointer=mxCreateNumericMatrix(N,1,mxDOUBLE_CLASS,mxREAL); 	mxSetField(out,0,names[idx],mxPointer);
	result->LAD1_L1LF=mxGetData(mxPointer); idx++;
	mxPointer=mxCreateNumericMatrix(N,1,mxDOUBLE_CLASS,mxREAL); 	mxSetField(out,0,names[idx],mxPointer);
	result->LAD2_L1LF=mxGetData(mxPointer); idx++;
	mxPointer=mxCreateNumericMatrix(N,1,mxDOUBLE_CLASS,mxREAL); 	mxSetField(out,0,names[idx],mxPointer);
	result->LAI_L2WL=mxGetData(mxPointer); idx++;
	mxPointer=mxCreateNumericMatrix(N,1,mxDOUBLE_CLASS,mxREAL); 	mxSetField(out,0,names[idx],mxPointer);
	result->AIC_L2WL=mxGetData(mxPointer); idx++;
	mxPointer=mxCreateNumericMatrix(N,1,mxDOUBLE_CLASS,mxREAL); 	mxSetField(out,0,names[idx],mxPointer);
	result->LAD1_L2WL=mxGetData(mxPointer); idx++;
	mxPointer=mxCreateNumericMatrix(N,1,mxDOUBLE_CLASS,mxREAL); 	mxSetField(out,0,names[idx],mxPointer);
	result->LAD2_L2WL=mxGetData(mxPointer); idx++;
	mxPointer=mxCreateNumericMatrix(N,1,mxDOUBLE_CLASS,mxREAL); 	mxSetField(out,0,names[idx],mxPointer);
	result->LAI_L1WL=mxGetData(mxPointer); idx++;
	mxPointer=mxCreateNumericMatrix(N,1,mxDOUBLE_CLASS,mxREAL); 	mxSetField(out,0,names[idx],mxPointer);
	result->AIC_L1WL=mxGetData(mxPointer); idx++;
	mxPointer=mxCreateNumericMatrix(N,1,mxDOUBLE_CLASS,mxREAL); 	mxSetField(out,0,names[idx],mxPointer);
	result->LAD1_L1WL=mxGetData(mxPointer); idx++;
	mxPointer=mxCreateNumericMatrix(N,1,mxDOUBLE_CLASS,mxREAL); 	mxSetField(out,0,names[idx],mxPointer);
	result->LAD2_L1WL=mxGetData(mxPointer); idx++;
	return out;
}
DllExport
#define pY   prhs[0]
#define pOpt prhs[2] 
void DllExport LAI_estimate(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[]);
void DllExport LAI_gfunc(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[]);
void DllExport LAI_simulate(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[]);
void DllExport mexFunction(int nlhs,mxArray * _restrict plhs[],int nrhs,const mxArray * _restrict prhs[])
{
	if (!mxIsStruct(pOpt))
	{
		r_error("Error (mexFunction):The input parameter OPTION should be a struct variable!\n");
		return ;
	}
	mxArray *tmpPointer=mxGetField(pOpt,0,"algorithm");
	if (tmpPointer==NULL)
	{
		r_error("Error (mexFunction): In the option,the algorithm parameter must be specified!\n");
		return;
	}
	char *algName=mxIsChar(tmpPointer) ? mxArrayToString(tmpPointer) : NULL;
	if (algName !=NULL)
	{
		if (strcicmp(algName,"LAI_estimate")==0)
		{
		LAI_estimate(nlhs,plhs,nrhs,prhs);
		}	
		else if (strcicmp(algName,"LAI_gfunc")==0)
		{
			LAI_gfunc(nlhs,plhs,nrhs,prhs);
		}
		else if (strcicmp(algName,"LAI_simulate")==0)
		{
			LAI_simulate(nlhs,plhs,nrhs,prhs);
		}		
		r_free(algName);
	}
	return;
}
#endif
