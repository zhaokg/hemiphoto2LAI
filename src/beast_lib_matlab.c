#include "abc_000_macro.h"
DISABLE_MANY_WARNINGS
static char fileID UNUSED_DECORATOR='M';
#if M_INTERFACE==1
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <float.h>
#include "abc_001_config.h"
#include "abc_blas_lapack_lib.h"
#include "abc_common.h"
#include "beast_common.h"
#include "beast_lib_matlab.h"
extern float determine_period(F32PTR data,int32_t N,float  omissionValue);
bool  get_options(Options * _restrict opt,const mxArray * _restrict pY,const mxArray * _restrict pOpt,MemPointers * _restrict MEM)
{
	void	*tmpPointer=NULL;
	int32_t  N=0;
	if (!mxIsStruct(pOpt))
	{
		MEM->free_all(MEM);
		r_error("Error:The input parameter OPTION should be a struct variable!\n");
		return false;
	}
	opt->inputFromDisk=false;
	opt->period=-999; 
	tmpPointer=mxGetField(pOpt,0,"period");
	if (tmpPointer !=NULL)
	{
		opt->period=(float)mxGetScalar((mxArray*)tmpPointer);
	}
	tmpPointer=mxGetField(pOpt,0,"omissionValue");
	opt->omissionValue=(tmpPointer !=NULL) ? (float)mxGetScalar((mxArray*)tmpPointer) : -999;
	if (mxIsChar(pY)) 
	{
		opt->isSingleyInput=true;
		opt->inputFromDisk=true;
		tmpPointer=(void *)mxArrayToString(pY);
		opt->inputFile=MEM->alloc(MEM,sizeof(char)*(strlen(tmpPointer)+1),0);
		strcpy(opt->inputFile,(char *)tmpPointer);
		r_free(tmpPointer);
		N=opt->N=(int)mxGetScalar(mxGetField(pOpt,0,"lengthPerTimeSeries_infile"));
		tmpPointer=fopen(opt->inputFile,"rb");
		if (opt->period < 0)
		{
			float *data;
			data=r_malloc(sizeof(float)*N);
			fread(data,sizeof(float),N,(FILE *)tmpPointer);
			opt->period=determine_period(data,N,opt->omissionValue);
			r_free(data);
			r_printf("THE PERIOD PARAMETER IS NOT SET. A BEST GUESS of IT IS%d AND USED IN THE DECOMPOSITION.\n",(int)opt->period);
		}
		fseek((FILE *)tmpPointer,0,SEEK_END);
		int64_t fileSizeinByte=ftell((FILE *)tmpPointer);
		fclose((FILE *)tmpPointer);
		opt->totalPixelNum=(uint32_t)(fileSizeinByte/(4 * N));
		r_printf("File Size%d ; totalPixelNum%d \n",fileSizeinByte,opt->totalPixelNum);
	}
	else
	{
		opt->isSingleyInput=(mxGetClassID(pY)==mxSINGLE_CLASS);
		opt->inputFromDisk=false;
		if ((int)min(mxGetM(pY),mxGetN(pY))==1)
			N=opt->N=(int)max(mxGetM(pY),mxGetN(pY)),opt->totalPixelNum=1;
		else
			N=opt->N=(int)mxGetM(pY),opt->totalPixelNum=(int)mxGetN(pY);
		opt->yInputData=(float *)mxGetData(pY);
		if (opt->period < 0) 
		{
			float *data=r_malloc(sizeof(float)*N);
			if (opt->isSingleyInput)
			{
				memcpy(data,opt->yInputData,sizeof(float)*N);
			}
			else
			{
				rF32PTR * yInput=opt->yInputData ;
				for (int i=1; i <=N; i++)
					*data++=(float)(*((double *)yInput)++);
				data=data - N;
			}
			opt->period=determine_period(data,N,opt->omissionValue);
			r_printf("THE PERIOD PARAMETER IS NOT SET. A BEST GUESS of IT IS%d AND USED IN THE DECOMPOSITION.\n",(int)opt->period);
			r_free(data);
		}
	}
	tmpPointer=mxGetField(pOpt,0,"startTime");
	opt->startTime=(tmpPointer !=NULL) ? (float)mxGetScalar((mxArray*)tmpPointer) : 1;
	tmpPointer=mxGetField(pOpt,0,"timeInterval");
	opt->timeInterval=(tmpPointer !=NULL) ? (float)mxGetScalar((mxArray*)tmpPointer) : 1;
	tmpPointer=mxGetField(pOpt,0,"minSeasonOrder");
	opt->minSeasonOrder=(tmpPointer !=NULL) ? (unsigned char)mxGetScalar((mxArray*)tmpPointer) : 1;
	tmpPointer=mxGetField(pOpt,0,"maxSeasonOrder");
	opt->maxSeasonOrder=(tmpPointer !=NULL) ? (unsigned char)mxGetScalar((mxArray*)tmpPointer) : (unsigned char)(opt->period/2) - 1;
	opt->maxSeasonOrder=(unsigned char)min(opt->period/2,opt->maxSeasonOrder);
	opt->maxSeasonOrder=max(opt->maxSeasonOrder,opt->minSeasonOrder);
	tmpPointer=mxGetField(pOpt,0,"minTrendOrder");
	opt->minTrendOrder=(tmpPointer !=NULL) ? (unsigned char)mxGetScalar((mxArray*)tmpPointer) : 0;
	tmpPointer=mxGetField(pOpt,0,"maxTrendOrder");
	opt->maxTrendOrder=(tmpPointer !=NULL) ? (unsigned char)mxGetScalar((mxArray*)tmpPointer) : 1;
	opt->maxTrendOrder=max(opt->maxTrendOrder,opt->minTrendOrder);
	tmpPointer=mxGetField(pOpt,0,"minSepDist_Trend");
	opt->minSepDist_Trend=(tmpPointer !=NULL) ? (uint16_t)mxGetScalar((mxArray*)tmpPointer) : (uint16_t)(opt->period/2);
	opt->minSepDist_Trend=max(opt->minSepDist_Trend,1);
	tmpPointer=mxGetField(pOpt,0,"minSepDist_Season");
	opt->minSepDist_Season=(tmpPointer !=NULL) ? (uint16_t)mxGetScalar((mxArray*)tmpPointer) : (uint16_t)(opt->period/2);
	opt->minSepDist_Season=max(opt->minSepDist_Season,1);
	tmpPointer=mxGetField(pOpt,0,"maxKnotNum_Trend");
	opt->maxKnotNum_Trend=(tmpPointer !=NULL) ? (uint16_t)mxGetScalar((mxArray*)tmpPointer) : (uint16_t)floor(N/(opt->minSepDist_Trend+1) - 1);
	opt->maxKnotNum_Trend=(uint16_t)min(floor(opt->N/(opt->minSepDist_Trend+1) - 1),opt->maxKnotNum_Trend);
	tmpPointer=mxGetField(pOpt,0,"maxKnotNum_Season");
	opt->maxKnotNum_Season=(tmpPointer !=NULL) ? (uint16_t)mxGetScalar((mxArray*)tmpPointer) : (uint16_t)floor(N/(opt->minSepDist_Season+1) - 1);
	opt->maxKnotNum_Season=(uint16_t)min(floor(opt->N/(opt->minSepDist_Season+1) - 1),opt->maxKnotNum_Season);
	tmpPointer=mxGetField(pOpt,0,"maxMoveStepSize");
	opt->maxMoveStepSize=(tmpPointer !=NULL) ? (uint16_t)mxGetScalar((mxArray*)tmpPointer) : (uint16_t)(opt->minSepDist_Trend/3+1);
	tmpPointer=mxGetField(pOpt,0,"samples");
	opt->samples=(tmpPointer !=NULL) ? (uint32_t)mxGetScalar((mxArray*)tmpPointer) : 10000;
	tmpPointer=mxGetField(pOpt,0,"thinningFactor");
	opt->thinningFactor=(tmpPointer !=NULL) ? (uint16_t)mxGetScalar((mxArray*)tmpPointer) : 1;
	tmpPointer=mxGetField(pOpt,0,"burnin");
	opt->burnin=(tmpPointer !=NULL) ? (uint32_t)mxGetScalar((mxArray*)tmpPointer) : 500;
	opt->burnin=max(opt->burnin,500);
	tmpPointer=mxGetField(pOpt,0,"chainNumber");
	opt->chainNumber=(tmpPointer !=NULL) ? (uint32_t)mxGetScalar((mxArray*)tmpPointer) : 5;
	opt->chainNumber=max(opt->chainNumber,1);
	tmpPointer=mxGetField(pOpt,0,"resamplingTrendOrderProb");
	opt->resamplingTrendOrderProb=(tmpPointer !=NULL) ? (float)mxGetScalar((mxArray*)tmpPointer) : 0.1f;
	tmpPointer=mxGetField(pOpt,0,"resamplingSeasonOrderProb");
	opt->resamplingSeasonOrderProb=(tmpPointer !=NULL) ? (float)mxGetScalar((mxArray*)tmpPointer) : 0.17f;
	tmpPointer=mxGetField(pOpt,0,"seed");
	opt->seed=(tmpPointer !=NULL) ? (uint64_t)mxGetScalar((mxArray*)tmpPointer) : 0;
	tmpPointer=mxGetField(pOpt,0,"outputToDisk");
	opt->outputToDisk=(tmpPointer !=NULL) ? (bool)mxGetScalar((mxArray*)tmpPointer) : false;
	if (opt->outputToDisk)
	{
		tmpPointer=mxGetField(pOpt,0,"outputFolder");
		tmpPointer=mxIsChar((mxArray *)tmpPointer) ? mxArrayToString((mxArray *)tmpPointer) : NULL;
		if (tmpPointer !=NULL)
			opt->outputFolder=MEM->alloc(MEM,sizeof(char)*(strlen(tmpPointer)+20),0),
			strcpy(opt->outputFolder,(char *)tmpPointer),
			r_free((char *)tmpPointer);
		else
			opt->outputToDisk=false;
	}
	tmpPointer=mxGetField(pOpt,0,"printToScreen");
	opt->printToScreen=(tmpPointer !=NULL) ? (bool)mxGetScalar((mxArray*)tmpPointer) : true;
	if (opt->printToScreen)
	{
		tmpPointer=mxGetField(pOpt,0,"printCharLen");
		opt->printCharLen=(tmpPointer !=NULL) ? (int16_t)mxGetScalar((mxArray*)tmpPointer) : 150;
	}
	tmpPointer=mxGetField(pOpt,0,"computeCredible");
	opt->computeCredible=(tmpPointer !=NULL) ? (bool)mxGetScalar((mxArray*)tmpPointer) : false;
	if (opt->computeCredible)
	{
		tmpPointer=mxGetField(pOpt,0,"fastCIComputation");
		opt->fastCIComputation=(tmpPointer !=NULL) ? (bool)mxGetScalar((mxArray*)tmpPointer) : true;
		tmpPointer=mxGetField(pOpt,0,"alphaLevel");
		opt->alphaLevel=(tmpPointer !=NULL) ? (float)mxGetScalar((mxArray*)tmpPointer) : 0.95;
	}
	tmpPointer=mxGetField(pOpt,0,"computeSlopeSign");
	opt->computeSlopeSign=(tmpPointer !=NULL) ? (bool)mxGetScalar((mxArray*)tmpPointer) : false;
	tmpPointer=mxGetField(pOpt,0,"ridgeFactor");
	opt->ridgeFactor=(tmpPointer !=NULL) ? (float)mxGetScalar((mxArray*)tmpPointer) : 0.00001f;
	tmpPointer=mxGetField(pOpt,0,"computeHarmonicOrder");
	opt->computeHarmonicOrder=(tmpPointer !=NULL) ? (bool)mxGetScalar((mxArray*)tmpPointer) : false;
	tmpPointer=mxGetField(pOpt,0,"computeTrendOrder");
	opt->computeTrendOrder=(tmpPointer !=NULL) ? (bool)mxGetScalar((mxArray*)tmpPointer) : false;
	tmpPointer=mxGetField(pOpt,0,"computeChangepoints");
	opt->computeChangepoints=(tmpPointer !=NULL) ? (bool)mxGetScalar((mxArray*)tmpPointer) : false;
	opt->Npad=(int32_t)ceil((float)opt->N/8.0f) * 8;
	opt->Npad16=(int32_t)ceil((float)opt->N/16.0f) * 16;
	opt->separator='.';
	return true;
}
mxArray *allocate_output(RESULT * _restrict matOutput,Options * _restrict opt,int * nptr)
{
	mxArray * _restrict out;
	mxArray * _restrict mxPointer;
	mwSize	 dim2[2]={ 1,1 };
	mwSize   N=opt->N;
	mwSize   K=opt->totalPixelNum;
	mwSize 	 dim3[3]={ N,2,K };
	char *fieldNames[]={ "time","sN","tN","sNProb","tNProb","sProb","tProb","s","sCI","sSD",\
		"t","tCI","tSD","b","bCI","bSD","marg_lik","sig2","bsign","horder","torder","tcp","scp","tcpCI","scpCI" };
	int  fieldNumber=25;
	out=mxCreateStructArray(2,dim2,fieldNumber,fieldNames);
	mxPointer=mxCreateNumericMatrix(N,1,mxSINGLE_CLASS,mxREAL); 		mxSetField(out,0,"time",mxPointer); 	matOutput->time=mxGetData(mxPointer);
	mxPointer=mxCreateNumericMatrix(1,K,mxSINGLE_CLASS,mxREAL); 		mxSetField(out,0,"sN",mxPointer); 	matOutput->sN=mxGetData(mxPointer);
	mxPointer=mxCreateNumericMatrix(1,K,mxSINGLE_CLASS,mxREAL); 		mxSetField(out,0,"tN",mxPointer);	matOutput->tN=mxGetData(mxPointer);
	mxPointer=mxCreateNumericMatrix(opt->maxKnotNum_Season+1,K,mxSINGLE_CLASS,mxREAL); 	mxSetField(out,0,"sNProb",mxPointer); 	matOutput->sNProb=mxGetData(mxPointer);
	mxPointer=mxCreateNumericMatrix(opt->maxKnotNum_Trend+1,K,mxSINGLE_CLASS,mxREAL); 	mxSetField(out,0,"tNProb",mxPointer);	matOutput->tNProb=mxGetData(mxPointer);
	mxPointer=mxCreateNumericMatrix(N,K,mxSINGLE_CLASS,mxREAL); 		mxSetField(out,0,"sProb",mxPointer);	matOutput->sProb=mxGetData(mxPointer);
	mxPointer=mxCreateNumericMatrix(N,K,mxSINGLE_CLASS,mxREAL); 		mxSetField(out,0,"tProb",mxPointer); matOutput->tProb=mxGetData(mxPointer);
	mxPointer=mxCreateNumericMatrix(N,K,mxSINGLE_CLASS,mxREAL); 		mxSetField(out,0,"s",mxPointer); 	matOutput->s=mxGetData(mxPointer);
	if (opt->computeCredible) mxPointer=mxCreateNumericArray(3,dim3,mxSINGLE_CLASS,mxREAL),mxSetField(out,0,"sCI",mxPointer),matOutput->sCI=mxGetData(mxPointer);
	mxPointer=mxCreateNumericMatrix(N,K,mxSINGLE_CLASS,mxREAL); 		mxSetField(out,0,"sSD",mxPointer); 	matOutput->sSD=mxGetData(mxPointer);
	mxPointer=mxCreateNumericMatrix(N,K,mxSINGLE_CLASS,mxREAL); 		mxSetField(out,0,"t",mxPointer); 	matOutput->t=mxGetData(mxPointer);
	if (opt->computeCredible) 	mxPointer=mxCreateNumericArray(3,dim3,mxSINGLE_CLASS,mxREAL),mxSetField(out,0,"tCI",mxPointer),matOutput->tCI=mxGetData(mxPointer);
	mxPointer=mxCreateNumericMatrix(N,K,mxSINGLE_CLASS,mxREAL); 		mxSetField(out,0,"tSD",mxPointer); 	matOutput->tSD=mxGetData(mxPointer);
	mxPointer=mxCreateNumericMatrix(N,K,mxSINGLE_CLASS,mxREAL); 		mxSetField(out,0,"b",mxPointer); 	matOutput->b=mxGetData(mxPointer);
	if (opt->computeCredible)  mxPointer=mxCreateNumericArray(3,dim3,mxSINGLE_CLASS,mxREAL),mxSetField(out,0,"bCI",mxPointer),matOutput->bCI=mxGetData(mxPointer);
	mxPointer=mxCreateNumericMatrix(N,K,mxSINGLE_CLASS,mxREAL); 		mxSetField(out,0,"bSD",mxPointer); 	matOutput->bSD=mxGetData(mxPointer);
	mxPointer=mxCreateNumericMatrix(1,K,mxSINGLE_CLASS,mxREAL); 		mxSetField(out,0,"marg_lik",mxPointer); matOutput->marg_lik=mxGetData(mxPointer);
	mxPointer=mxCreateNumericMatrix(1,K,mxSINGLE_CLASS,mxREAL); 		mxSetField(out,0,"sig2",mxPointer);    matOutput->sig2=mxGetData(mxPointer);
	if (opt->computeSlopeSign) 	   mxPointer=mxCreateNumericMatrix(N,K,mxSINGLE_CLASS,mxREAL),mxSetField(out,0,"bsign",mxPointer),matOutput->bsign=mxGetData(mxPointer);
	if (opt->computeHarmonicOrder) mxPointer=mxCreateNumericMatrix(N,K,mxSINGLE_CLASS,mxREAL),mxSetField(out,0,"horder",mxPointer),matOutput->horder=mxGetData(mxPointer);
	if (opt->computeTrendOrder) mxPointer=mxCreateNumericMatrix(N,K,mxSINGLE_CLASS,mxREAL),mxSetField(out,0,"torder",mxPointer),matOutput->torder=mxGetData(mxPointer);
	if (opt->computeChangepoints)
	{
		mwSize 	 dim1_3[3]={ opt->maxKnotNum_Season,2,K };
		mxPointer=mxCreateNumericMatrix(opt->maxKnotNum_Season,K,mxSINGLE_CLASS,mxREAL); 		mxSetField(out,0,"scp",mxPointer);    matOutput->scp=mxGetData(mxPointer);
		mxPointer=mxCreateNumericArray(3,dim1_3,mxSINGLE_CLASS,mxREAL); 						mxSetField(out,0,"scpCI",mxPointer);  matOutput->scpCI=mxGetData(mxPointer);
		mwSize 	 dim2_3[3]={ opt->maxKnotNum_Trend,2,K };
		mxPointer=mxCreateNumericMatrix(opt->maxKnotNum_Trend,K,mxSINGLE_CLASS,mxREAL); 		mxSetField(out,0,"tcp",mxPointer);    matOutput->tcp=mxGetData(mxPointer);
		mxPointer=mxCreateNumericArray(3,dim2_3,mxSINGLE_CLASS,mxREAL); 	 					mxSetField(out,0,"tcpCI",mxPointer);  matOutput->tcpCI=mxGetData(mxPointer);
	}
	fill_1_to_N(matOutput->time,N);
	r_ippsMulC_32f_I(opt->timeInterval,matOutput->time,(const int)N);
	r_ippsSubC_32f_I(-(opt->startTime - opt->timeInterval),matOutput->time,(const int)N);
	return out;
}
bool   get_options_trend(Options * _restrict opt,const mxArray * _restrict pY,const mxArray * _restrict pOpt,MemPointers * _restrict MEM)
{
	void	*tmpPointer;
	int32_t  N=0;
	if (!mxIsStruct(pOpt))
	{
		MEM->free_all(MEM);
		r_error("Error:The input parameter OPTION should be a struct variable!\n");
		return false;
	}
	opt->inputFromDisk=false;
	opt->period=-999; 
	tmpPointer=mxGetField(pOpt,0,"period");
	opt->period=(tmpPointer !=NULL) ? (float)mxGetScalar((mxArray*)tmpPointer) : -999;
	tmpPointer=mxGetField(pOpt,0,"omissionValue");
	opt->omissionValue=(tmpPointer !=NULL) ? (float)mxGetScalar((mxArray*)tmpPointer) : -999;
	if (mxIsChar(pY)) 
	{
		opt->isSingleyInput=true;
		opt->inputFromDisk=true;
		tmpPointer=(void *)mxArrayToString(pY);
		opt->inputFile=MEM->alloc(MEM,sizeof(char)*(strlen(tmpPointer)+1),0);
		strcpy(opt->inputFile,(char *)tmpPointer);
		r_free(tmpPointer);
		N=opt->N=(int)mxGetScalar(mxGetField(pOpt,0,"lengthPerTimeSeries_infile"));
		tmpPointer=fopen(opt->inputFile,"rb");
		fseek((FILE *)tmpPointer,0,SEEK_END);
		int64_t fileSizeinByte=ftell((FILE *)tmpPointer);
		fclose((FILE *)tmpPointer);
		opt->totalPixelNum=(uint32_t)(fileSizeinByte/(4 * N));
		r_printf("File Size%d ; totalPixelNum%d \n",fileSizeinByte,opt->totalPixelNum);
	}
	else
	{
		opt->isSingleyInput=(mxGetClassID(pY)==mxSINGLE_CLASS);
		opt->inputFromDisk=false;
		if ((int)min(mxGetM(pY),mxGetN(pY))==1)
			N=opt->N=(int)max(mxGetM(pY),mxGetN(pY)),opt->totalPixelNum=1;
		else
			N=opt->N=(int)mxGetM(pY),opt->totalPixelNum=(int)mxGetN(pY);
		opt->yInputData=(float *)mxGetData(pY);
	}
	tmpPointer=mxGetField(pOpt,0,"startTime");
	opt->startTime=(tmpPointer !=NULL) ? (float)mxGetScalar((mxArray*)tmpPointer) : 1;
	tmpPointer=mxGetField(pOpt,0,"timeInterval");
	opt->timeInterval=(tmpPointer !=NULL) ? (float)mxGetScalar((mxArray*)tmpPointer) : 1;
	tmpPointer=mxGetField(pOpt,0,"minTrendOrder");
	opt->minTrendOrder=(tmpPointer !=NULL) ? (unsigned char)mxGetScalar((mxArray*)tmpPointer) : 0;
	tmpPointer=mxGetField(pOpt,0,"maxTrendOrder");
	opt->maxTrendOrder=(tmpPointer !=NULL) ? (unsigned char)mxGetScalar((mxArray*)tmpPointer) : 1;
	opt->maxTrendOrder=max(opt->maxTrendOrder,opt->minTrendOrder);
	tmpPointer=mxGetField(pOpt,0,"minSepDist_Trend");
	opt->minSepDist_Trend=(tmpPointer !=NULL) ? (uint16_t)mxGetScalar((mxArray*)tmpPointer) : (uint16_t)(opt->period/2);
	tmpPointer=mxGetField(pOpt,0,"maxKnotNum_Trend");
	opt->maxKnotNum_Trend=(tmpPointer !=NULL) ? (uint16_t)mxGetScalar((mxArray*)tmpPointer) : (uint16_t)floor(N/(opt->minSepDist_Trend+1) - 1);
	opt->maxKnotNum_Trend=(uint16_t)min(floor(opt->N/(opt->minSepDist_Trend+1) - 1),opt->maxKnotNum_Trend);
	tmpPointer=mxGetField(pOpt,0,"maxMoveStepSize");
	opt->maxMoveStepSize=(tmpPointer !=NULL) ? (uint16_t)mxGetScalar((mxArray*)tmpPointer) : (uint16_t)(opt->minSepDist_Trend/3+1);
	tmpPointer=mxGetField(pOpt,0,"samples");
	opt->samples=(tmpPointer !=NULL) ? (uint32_t)mxGetScalar((mxArray*)tmpPointer) : 10000;
	tmpPointer=mxGetField(pOpt,0,"thinningFactor");
	opt->thinningFactor=(tmpPointer !=NULL) ? (uint16_t)mxGetScalar((mxArray*)tmpPointer) : 1;
	tmpPointer=mxGetField(pOpt,0,"burnin");
	opt->burnin=(tmpPointer !=NULL) ? (uint32_t)mxGetScalar((mxArray*)tmpPointer) : 500;
	opt->burnin=max(opt->burnin,500);
	tmpPointer=mxGetField(pOpt,0,"chainNumber");
	opt->chainNumber=(tmpPointer !=NULL) ? (uint32_t)mxGetScalar((mxArray*)tmpPointer) : 5;
	opt->chainNumber=max(opt->chainNumber,1);
	tmpPointer=mxGetField(pOpt,0,"resamplingTrendOrderProb");
	opt->resamplingTrendOrderProb=(tmpPointer !=NULL) ? (float)mxGetScalar((mxArray*)tmpPointer) : 0.1f;
	tmpPointer=mxGetField(pOpt,0,"resamplingSeasonOrderProb");
	opt->resamplingSeasonOrderProb=(tmpPointer !=NULL) ? (float)mxGetScalar((mxArray*)tmpPointer) : 0.17f;
	tmpPointer=mxGetField(pOpt,0,"seed");
	opt->seed=(tmpPointer !=NULL) ? (uint64_t)mxGetScalar((mxArray*)tmpPointer) : 0;
	tmpPointer=mxGetField(pOpt,0,"outputToDisk");
	opt->outputToDisk=(tmpPointer !=NULL) ? (bool)mxGetScalar((mxArray*)tmpPointer) : false;
	if (opt->outputToDisk)
	{
		tmpPointer=mxGetField(pOpt,0,"outputFolder");
		tmpPointer=mxIsChar((mxArray *)tmpPointer) ? mxArrayToString((mxArray *)tmpPointer) : NULL;
		if (tmpPointer !=NULL)
			opt->outputFolder=MEM->alloc(MEM,sizeof(char)*(strlen(tmpPointer)+20),0),
			strcpy(opt->outputFolder,(char *)tmpPointer),
			r_free((char *)tmpPointer);
		else
			opt->outputToDisk=false;
	}
	tmpPointer=mxGetField(pOpt,0,"printToScreen");
	opt->printToScreen=(tmpPointer !=NULL) ? (bool)mxGetScalar((mxArray*)tmpPointer) : true;
	if (opt->printToScreen)
	{
		tmpPointer=mxGetField(pOpt,0,"printCharLen");
		opt->printCharLen=(tmpPointer !=NULL) ? (int16_t)mxGetScalar((mxArray*)tmpPointer) : 150;
	}
	tmpPointer=mxGetField(pOpt,0,"computeCredible");
	opt->computeCredible=(tmpPointer !=NULL) ? (bool)mxGetScalar((mxArray*)tmpPointer) : false;
	if (opt->computeCredible)
	{
		tmpPointer=mxGetField(pOpt,0,"fastCIComputation");
		opt->fastCIComputation=(tmpPointer !=NULL) ? (bool)mxGetScalar((mxArray*)tmpPointer) : true;
		tmpPointer=mxGetField(pOpt,0,"alphaLevel");
		opt->alphaLevel=(tmpPointer !=NULL) ? (float)mxGetScalar((mxArray*)tmpPointer) : 0.95;
	}
	tmpPointer=mxGetField(pOpt,0,"computeSlopeSign");
	opt->computeSlopeSign=(tmpPointer !=NULL) ? (bool)mxGetScalar((mxArray*)tmpPointer) : false;
	tmpPointer=mxGetField(pOpt,0,"computeTrendOrder");
	opt->computeTrendOrder=(tmpPointer !=NULL) ? (bool)mxGetScalar((mxArray*)tmpPointer) : false;
	tmpPointer=mxGetField(pOpt,0,"computeChangepoints");
	opt->computeChangepoints=(tmpPointer !=NULL) ? (bool)mxGetScalar((mxArray*)tmpPointer) : false;
	tmpPointer=mxGetField(pOpt,0,"ridgeFactor");
	opt->ridgeFactor=(tmpPointer !=NULL) ? (float)mxGetScalar((mxArray*)tmpPointer) : 0.00001f;
	opt->Npad=(int32_t)ceil((float)opt->N/8.0f) * 8;
	opt->Npad16=(int32_t)ceil((float)opt->N/16.0f) * 16;
	opt->separator='.';
	return true;
}
mxArray *allocate_output_trend(RESULT * _restrict matOutput,Options * _restrict opt,int * nptr)
{
	mxArray * _restrict out;
	mxArray * _restrict mxPointer;
	mwSize	 dim2[2]={ 1,1 };
	mwSize   N=opt->N;
	mwSize   K=opt->totalPixelNum;
	mwSize 	 dim3[3]={ N,2,K };
	char *fieldNames[]={ "time","tN","tNProb","tProb","t","tCI","tSD","b","bCI","bSD","marg_lik","sig2","bsign","torder",
		"tcp","tcpCI" };
	int  fieldNumber=13+1+2;
	out=mxCreateStructArray(2,dim2,fieldNumber,fieldNames);
	mxPointer=mxCreateNumericMatrix(N,1,mxSINGLE_CLASS,mxREAL); 		mxSetField(out,0,"time",mxPointer); 	matOutput->time=mxGetData(mxPointer);
	mxPointer=mxCreateNumericMatrix(1,K,mxSINGLE_CLASS,mxREAL); 		mxSetField(out,0,"tN",mxPointer);	matOutput->tN=mxGetData(mxPointer);
	mxPointer=mxCreateNumericMatrix(opt->maxKnotNum_Trend+1,K,mxSINGLE_CLASS,mxREAL); 	mxSetField(out,0,"tNProb",mxPointer);	matOutput->tNProb=mxGetData(mxPointer);
	mxPointer=mxCreateNumericMatrix(N,K,mxSINGLE_CLASS,mxREAL); 		mxSetField(out,0,"tProb",mxPointer); matOutput->tProb=mxGetData(mxPointer);
	mxPointer=mxCreateNumericMatrix(N,K,mxSINGLE_CLASS,mxREAL); 		mxSetField(out,0,"t",mxPointer); 	matOutput->t=mxGetData(mxPointer);
	if (opt->computeSlopeSign) mxPointer=mxCreateNumericArray(3,dim3,mxSINGLE_CLASS,mxREAL),mxSetField(out,0,"tCI",mxPointer),matOutput->tCI=mxGetData(mxPointer);
	mxPointer=mxCreateNumericMatrix(N,K,mxSINGLE_CLASS,mxREAL); 		mxSetField(out,0,"tSD",mxPointer); 	matOutput->tSD=mxGetData(mxPointer);
	mxPointer=mxCreateNumericMatrix(N,K,mxSINGLE_CLASS,mxREAL); 		mxSetField(out,0,"b",mxPointer); 	matOutput->b=mxGetData(mxPointer);
	if (opt->computeSlopeSign) mxPointer=mxCreateNumericArray(3,dim3,mxSINGLE_CLASS,mxREAL),mxSetField(out,0,"bCI",mxPointer),matOutput->bCI=mxGetData(mxPointer);
	mxPointer=mxCreateNumericMatrix(N,K,mxSINGLE_CLASS,mxREAL); 		mxSetField(out,0,"bSD",mxPointer); 	matOutput->bSD=mxGetData(mxPointer);
	mxPointer=mxCreateNumericMatrix(1,K,mxSINGLE_CLASS,mxREAL); 		mxSetField(out,0,"marg_lik",mxPointer); matOutput->marg_lik=mxGetData(mxPointer);
	mxPointer=mxCreateNumericMatrix(1,K,mxSINGLE_CLASS,mxREAL); 		mxSetField(out,0,"sig2",mxPointer);    matOutput->sig2=mxGetData(mxPointer);
	if (opt->computeSlopeSign) 	mxPointer=mxCreateNumericMatrix(N,K,mxSINGLE_CLASS,mxREAL),mxSetField(out,0,"bsign",mxPointer),matOutput->bsign=mxGetData(mxPointer);
	if (opt->computeTrendOrder) {
		mxPointer=mxCreateNumericMatrix(N,K,mxSINGLE_CLASS,mxREAL),mxSetField(out,0,"torder",mxPointer),matOutput->torder=mxGetData(mxPointer);
	}
	if (opt->computeChangepoints)
	{
		mwSize 	 dim2_3[3]={ opt->maxKnotNum_Trend,2,K };
		mxPointer=mxCreateNumericMatrix(opt->maxKnotNum_Trend,K,mxSINGLE_CLASS,mxREAL); 		mxSetField(out,0,"tcp",mxPointer);    matOutput->tcp=mxGetData(mxPointer);
		mxPointer=mxCreateNumericArray(3,dim2_3,mxSINGLE_CLASS,mxREAL); 	 					mxSetField(out,0,"tcpCI",mxPointer);  matOutput->tcpCI=mxGetData(mxPointer);
	}
	fill_1_to_N(matOutput->time,N);
	r_ippsMulC_32f_I(opt->timeInterval,matOutput->time,(const int)N);
	r_ippsSubC_32f_I(-(opt->startTime - opt->timeInterval),matOutput->time,(const int)N);
	return out;
}
DllExport
#define pY   prhs[0]
#define pOpt prhs[1] 
void DllExport beastST_multipleChain_fast(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[]);
void DllExport beastTrend_multipleChain_fast(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[]);
void DllExport mvST_multipleChain(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[]);
void DllExport SBM_ST(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[]);
void DllExport SBM_ST_BIC(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[]);
void DllExport WinMainDemoTrend(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[]);
void DllExport WinMainDemoST(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[]);
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
		if (strcicmp(algName,"beastST_old")==0)
		{			
		}
		else if (strcicmp(algName,"beastST_multipleChain_fast")==0)
		{
			beastST_multipleChain_fast(nlhs,plhs,nrhs,prhs);
		}
		else if (strcicmp(algName,"beastTrend_multipleChain_fast")==0)
		{
			beastTrend_multipleChain_fast(nlhs,plhs,nrhs,prhs);
		}
		else if (strcicmp(algName,"mvST_multipleChain")==0)
		{
			mvST_multipleChain(nlhs,plhs,nrhs,prhs);
		}
		else if (strcicmp(algName,"SBM_ST")==0)
		{
			SBM_ST(nlhs,plhs,nrhs,prhs);
		}
		else if (strcicmp(algName,"SBM_ST_BIC")==0)
		{
			SBM_ST_BIC(nlhs,plhs,nrhs,prhs);
		}
		else if (strcicmp(algName,"WinMainDemoTrend")==0)
		{
			WinMainDemoTrend(nlhs,plhs,nrhs,prhs);
		}
		else if (strcicmp(algName,"WinMainDemoST")==0)
		{
			WinMainDemoST(nlhs,plhs,nrhs,prhs);
		}
		r_free(algName);
	}
	if (algName==NULL)
	{
	}
	return;
}
#endif
ENABLE_MANY_WARNINGS
