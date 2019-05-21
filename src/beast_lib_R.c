#include "abc_000_macro.h"
DISABLE_MANY_WARNINGS
char fileID UNUSED_DECORATOR='a';
#if R_INTERFACE==1
#include "beast_lib_R.h"
extern float determine_period(F32PTR data,int32_t N,float  omissionValue);
int nlhs=1;
SEXP getListElement(SEXP list,const char *str)
{
	SEXP elmt=NULL; 
	SEXP names=getAttrib(list,R_NamesSymbol);
	for (int i=0; i < length(list); i++)
	if (strcmp(CHAR(STRING_ELT(names,i)),str)==0) {
		elmt=VECTOR_ELT(list,i);
		break;
	}
	return elmt;
}
bool get_options(Options * _restrict opt,SEXP pY,SEXP pOpt,MemPointers * _restrict MEM)
{
	SEXP	tmpSEXP=NULL;
	void	*tmpPointer;
	int nprt=0;
	int N=-1;
	if (!isNewList(pOpt))
	{
		MEM->free_all(MEM);
		r_error("Error: The input parameter OPTION should be a List variable!\n");
	}
	bool isPeriodParGiven=false;
	opt->period=-9999; 
	tmpSEXP=getListElement(pOpt,"period");
	if (tmpSEXP !=NULL)
	{
		PROTECT(tmpSEXP=coerceVector(tmpSEXP,REALSXP)); nprt++;
		opt->period=(float)asReal(tmpSEXP);
	}
	tmpSEXP=getListElement(pOpt,"omissionValue");
	opt->omissionValue=(tmpSEXP !=NULL) ? (PROTECT(tmpSEXP=coerceVector(tmpSEXP,REALSXP)),nprt++,(float)asReal(tmpSEXP)) : -9999;
	if (TYPEOF(pY)==STRSXP)
	{	
		opt->inputFromDisk=true;
		char const *tmp=CHAR(STRING_ELT(pY,0));
		opt->inputFile=MEM->alloc(MEM,sizeof(char)*(strlen(tmp)+1),0);
		strcpy(opt->inputFile,tmp);
		tmpSEXP=getListElement(pOpt,"lengthPerTimeSeries_infile");
		if (tmpSEXP==NULL)
		{
			MEM->free_all(MEM);
			r_error("Error: If data ara read from a binary file,the 'opt$lengthPerTimeSeries_infile' must be set"
				" to tell how many sample points a time series has (i.e.,the length of time series.");
		}
		PROTECT(tmpSEXP=coerceVector(tmpSEXP,INTSXP));  nprt++;
		N=INTEGER_VALUE(tmpSEXP);
		tmpPointer=fopen(opt->inputFile,"rb");
		if (isPeriodParGiven==false)
		{
			float *data=(float *)r_malloc(sizeof(float)*N);
			rI32 readBytes=fread(data,sizeof(float),N,(FILE *)tmpPointer);
			opt->period=determine_period(data,N,opt->omissionValue);
			r_free(data);
			r_printf("THE PERIOD PARAMETER IS NOT SET. A BEST GUESS of IT IS%d AND USED IN THE DECOMPOSITION.\n",(int)opt->period);
		}
		fseek((FILE *)tmpPointer,0,SEEK_END);
		int64_t fileSizeinByte=ftell((FILE *)tmpPointer);
		fclose((FILE *)tmpPointer);
		opt->totalPixelNum=(uint32_t)(fileSizeinByte/(4 * N));
		r_printf("File Size:%d ; Number of time series/pixels:%d \n",fileSizeinByte,opt->totalPixelNum);
	}
	else if (TYPEOF(pY)==REALSXP)
	{ 
		opt->inputFromDisk=false;
		opt->yInputData=(float *) REAL(pY);
		if (isMatrix(pY))
		{
			SEXP dims=getAttrib(pY,R_DimSymbol);
			N=(int)INTEGER(dims)[0];
			opt->totalPixelNum=(int)INTEGER(dims)[1];
			int totalElement=N*opt->totalPixelNum;
		}
		else if (isReal(pY))
		{
			N=Rf_length(pY);
			opt->totalPixelNum=1;
		}
		if (opt->period <0) 
		{
			float *data=(float *)r_malloc(sizeof(float)*N);
			{
				double *tmpDouble=(double *) opt->yInputData;
				for (int i=1; i <=N; i++)
					*data++=(float)(*tmpDouble++);
				data=data - N;
			}			
			opt->period=determine_period(data,N,opt->omissionValue);			
			r_printf("THE PERIOD PARAMETER IS NOT SET. A BEST GUESS of IT IS%d AND USED IN THE DECOMPOSITION.\n",(int)opt->period);
			r_free(data);
		}
	}
	opt->isSingleyInput=0;
	tmpSEXP=getListElement(pOpt,"startTime");
	opt->startTime=(tmpSEXP !=NULL) ? (PROTECT(tmpSEXP=coerceVector(tmpSEXP,REALSXP)),nprt++,(float)asReal(tmpSEXP)) : 1.f;
	tmpSEXP=getListElement(pOpt,"timeInterval");
	opt->timeInterval=(tmpSEXP !=NULL) ? (PROTECT(tmpSEXP=coerceVector(tmpSEXP,REALSXP)),nprt++,(float)asReal(tmpSEXP)) : 1.f;
	tmpSEXP=getListElement(pOpt,"minSeasonOrder");
	opt->minSeasonOrder=(tmpSEXP !=NULL) ? (PROTECT(tmpSEXP=coerceVector(tmpSEXP,INTSXP)),nprt++,(unsigned char)asInteger(tmpSEXP)) : 1;
	tmpSEXP=getListElement(pOpt,"maxSeasonOrder");
	opt->maxSeasonOrder=(tmpSEXP !=NULL) ? (PROTECT(tmpSEXP=coerceVector(tmpSEXP,INTSXP)),nprt++,(unsigned char)asInteger(tmpSEXP)) : (unsigned char)(opt->period/2) - 1;
	opt->maxSeasonOrder=(unsigned char)min(opt->period/2,opt->maxSeasonOrder);
	opt->maxSeasonOrder=max(opt->maxSeasonOrder,opt->minSeasonOrder);
	tmpSEXP=getListElement(pOpt,"minTrendOrder");
	opt->minTrendOrder=(tmpSEXP !=NULL) ? (PROTECT(tmpSEXP=coerceVector(tmpSEXP,INTSXP)),nprt++,(unsigned char)asInteger(tmpSEXP)) : 0;
	tmpSEXP=getListElement(pOpt,"maxTrendOrder");
	opt->maxTrendOrder=(tmpSEXP !=NULL) ? (PROTECT(tmpSEXP=coerceVector(tmpSEXP,INTSXP)),nprt++,(unsigned char)asInteger(tmpSEXP)) : 1;
	opt->maxTrendOrder=max(opt->maxTrendOrder,opt->minTrendOrder);
	tmpSEXP=getListElement(pOpt,"minSepDist_Trend");
	opt->minSepDist_Trend=(tmpSEXP !=NULL) ? (PROTECT(tmpSEXP=coerceVector(tmpSEXP,INTSXP)),nprt++,(uint16_t)asInteger(tmpSEXP)) : (uint16_t)(opt->period/2);
	tmpSEXP=getListElement(pOpt,"minSepDist_Season");
	opt->minSepDist_Season=(tmpSEXP !=NULL) ? (PROTECT(tmpSEXP=coerceVector(tmpSEXP,INTSXP)),nprt++,(uint16_t)asInteger(tmpSEXP)) : (uint16_t)(opt->period/2);
	tmpSEXP=getListElement(pOpt,"maxKnotNum_Trend");
	opt->maxKnotNum_Trend=(tmpSEXP !=NULL) ? (PROTECT(tmpSEXP=coerceVector(tmpSEXP,INTSXP)),nprt++,(uint16_t)asInteger(tmpSEXP)) : (uint16_t)floor(N/(opt->minSepDist_Trend+1) - 1);
	opt->maxKnotNum_Trend=(uint16_t)min(floor((float)N/(opt->minSepDist_Trend+1) - 1),opt->maxKnotNum_Trend);
	tmpSEXP=getListElement(pOpt,"maxKnotNum_Season");
	opt->maxKnotNum_Season=(tmpSEXP !=NULL) ? (PROTECT(tmpSEXP=coerceVector(tmpSEXP,INTSXP)),nprt++,(uint16_t)asInteger(tmpSEXP)) : (uint16_t)floor(N/(opt->minSepDist_Season+1) - 1);
	opt->maxKnotNum_Season=(uint16_t)min(floor((float)N/(opt->minSepDist_Season+1) - 1),opt->maxKnotNum_Season);
	tmpSEXP=getListElement(pOpt,"maxMoveStepSize");
	opt->maxMoveStepSize=(tmpSEXP !=NULL) ? (PROTECT(tmpSEXP=coerceVector(tmpSEXP,INTSXP)),nprt++,(uint16_t)asInteger(tmpSEXP)) : (uint16_t)(opt->minSepDist_Trend/3+1);
	tmpSEXP=getListElement(pOpt,"samples");
	opt->samples=(tmpSEXP !=NULL) ? (PROTECT(tmpSEXP=coerceVector(tmpSEXP,INTSXP)),nprt++,(mwSize)asInteger(tmpSEXP)) : 10000;
	tmpSEXP=getListElement(pOpt,"thinningFactor");
	opt->thinningFactor=(tmpSEXP !=NULL) ? (PROTECT(tmpSEXP=coerceVector(tmpSEXP,INTSXP)),nprt++,(uint16_t)asInteger(tmpSEXP)) : 1;
	tmpSEXP=getListElement(pOpt,"burnin");
	opt->burnin=(tmpSEXP !=NULL) ? (PROTECT(tmpSEXP=coerceVector(tmpSEXP,INTSXP)),nprt++,(mwSize)asInteger(tmpSEXP)) : 500;
	opt->burnin=max(opt->burnin,500);
	tmpSEXP=getListElement(pOpt,"chainNumber");
	opt->chainNumber=(tmpSEXP !=NULL) ? (PROTECT(tmpSEXP=coerceVector(tmpSEXP,INTSXP)),nprt++,(mwSize)asInteger(tmpSEXP)) : 5;
	opt->chainNumber=max(opt->chainNumber,1);
	tmpSEXP=getListElement(pOpt,"resamplingTrendOrderProb");
	opt->resamplingTrendOrderProb=(tmpSEXP !=NULL) ? (PROTECT(tmpSEXP=coerceVector(tmpSEXP,REALSXP)),nprt++,(float)asReal(tmpSEXP)) : 0.1f;
	tmpSEXP=getListElement(pOpt,"resamplingSeasonOrderProb");
	opt->resamplingSeasonOrderProb=(tmpSEXP !=NULL) ? (PROTECT(tmpSEXP=coerceVector(tmpSEXP,REALSXP)),nprt++,(float)asReal(tmpSEXP)) : 0.17f;
	tmpSEXP=getListElement(pOpt,"seed");
	opt->seed=(tmpSEXP !=NULL) ? (PROTECT(tmpSEXP=coerceVector(tmpSEXP,INTSXP)),nprt++,(uint64_t)asInteger(tmpSEXP)) : 0;
	tmpSEXP=getListElement(pOpt,"outputToDisk");
	opt->outputToDisk=(tmpSEXP !=NULL) ? (PROTECT(tmpSEXP=coerceVector(tmpSEXP,INTSXP)),nprt++,(bool)asInteger(tmpSEXP)) : false;
	if (opt->outputToDisk)
	{
		tmpSEXP=getListElement(pOpt,"outputFolder");
		if (tmpSEXP !=NULL)
		{
			char const *tmpChar=CHAR(asChar(tmpSEXP));
			opt->outputFolder=MEM->alloc(MEM,sizeof(char)*(strlen(tmpChar)+20),0);
			strcpy(opt->outputFolder,tmpChar);
		}
		else
			opt->outputToDisk=false;
	}
	if (!opt->outputToDisk && nlhs==0)
	{
		MEM->free_all(MEM);
		r_error("To run,you must save the output to either a lefthand-side"
			"variable (e.g.,out=func(...) ) or specify the output path \n"
			"in the Option.outputFolder parameter. \n");
		return 0;
	}
	tmpSEXP=getListElement(pOpt,"printToScreen");
	opt->printToScreen=(tmpSEXP !=NULL) ? (PROTECT(tmpSEXP=coerceVector(tmpSEXP,INTSXP)),nprt++,(bool)asInteger(tmpSEXP)) : false;
	if (opt->printToScreen)
	{
		tmpSEXP=getListElement(pOpt,"printCharLen");
		opt->printCharLen=(tmpSEXP !=NULL) ? (PROTECT(tmpSEXP=coerceVector(tmpSEXP,INTSXP)),nprt++,(int16_t)asInteger(tmpSEXP)) : 80;
	}
	tmpSEXP=getListElement(pOpt,"computeCredible");
	opt->computeCredible=(tmpSEXP !=NULL) ? (PROTECT(tmpSEXP=coerceVector(tmpSEXP,INTSXP)),nprt++,(bool)asInteger(tmpSEXP)) : false;
	if (opt->computeCredible)
	{
		tmpSEXP=getListElement(pOpt,"fastCIComputation");
		opt->fastCIComputation=(tmpSEXP !=NULL) ? (PROTECT(tmpSEXP=coerceVector(tmpSEXP,INTSXP)),nprt++,(bool)asInteger(tmpSEXP)) : true;
		tmpSEXP=getListElement(pOpt,"alphaLevel");
		opt->alphaLevel=(tmpSEXP !=NULL) ? (PROTECT(tmpSEXP=coerceVector(tmpSEXP,INTSXP)),nprt++,(float)asReal(tmpSEXP)) : 0.95;
	}
	tmpSEXP=getListElement(pOpt,"computeSlopeSign");
	opt->computeSlopeSign=(tmpSEXP !=NULL) ? (PROTECT(tmpSEXP=coerceVector(tmpSEXP,INTSXP)),nprt++,(bool)asInteger(tmpSEXP)) : false;
	tmpSEXP=getListElement(pOpt,"computeHarmonicOrder");
	opt->computeHarmonicOrder=(tmpSEXP !=NULL) ? (PROTECT(tmpSEXP=coerceVector(tmpSEXP,INTSXP)),nprt++,(bool)asInteger(tmpSEXP)) : false;
	tmpSEXP=getListElement(pOpt,"computeTrendOrder");
	opt->computeTrendOrder=(tmpSEXP !=NULL) ? (PROTECT(tmpSEXP=coerceVector(tmpSEXP,INTSXP)),nprt++,(bool)asInteger(tmpSEXP)) : false;
	tmpSEXP=getListElement(pOpt,"computeChangepoints");
	opt->computeChangepoints=(tmpSEXP !=NULL) ? (PROTECT(tmpSEXP=coerceVector(tmpSEXP,INTSXP)),nprt++,(bool)asInteger(tmpSEXP)) : false;
	tmpSEXP=getListElement(pOpt,"ridgeFactor");
	opt->ridgeFactor=(float) (tmpSEXP !=NULL) ? (PROTECT(tmpSEXP=coerceVector(tmpSEXP,REALSXP)),nprt++,(float)asReal(tmpSEXP)) : 0.0001f;
	UNPROTECT(nprt);
	nprt=0;
	opt->Npad=(int32_t)ceil((float)N/8.0f) * 8;
	opt->Npad16=(int32_t)ceil((float)N/16.f) * 16;
	opt->N=N;
	opt->separator='$';
	return true;
}
SEXP allocate_output( RESULT * _restrict matOutput,Options * _restrict opt,int * _restrict nprt)
{
	SEXP ANS=NULL;
	SEXP ansListNames=NULL;
	int N=opt->N;
	int totalPixelNum=opt->totalPixelNum;
	SEXP tmpSEXP;
	int idx;
	int32_t elementNum=15+(int32_t)opt->computeCredible * 3+(int32_t)opt->computeSlopeSign  \
		+(int32_t)opt->computeHarmonicOrder+(int32_t)opt->computeTrendOrder * 1+(int32_t)opt->computeChangepoints * 4;
	PROTECT(ANS=allocVector(VECSXP,elementNum));++*nprt;
	PROTECT(ansListNames=allocVector(STRSXP,elementNum));++*nprt;
	idx=0; PROTECT(tmpSEXP=allocVector(REALSXP,N));++*nprt; 
	matOutput->time=(float *)REAL(tmpSEXP); SET_VECTOR_ELT(ANS,idx,tmpSEXP); SET_STRING_ELT(ansListNames,idx,mkChar("time"));
	idx++; PROTECT(tmpSEXP=allocVector(REALSXP,totalPixelNum));++*nprt; 
	matOutput->sN=(float *)REAL(tmpSEXP); SET_VECTOR_ELT(ANS,idx,tmpSEXP); SET_STRING_ELT(ansListNames,idx,mkChar("sN"));
	idx++; PROTECT(tmpSEXP=allocVector(REALSXP,totalPixelNum));++*nprt;
	matOutput->tN=(float *)REAL(tmpSEXP); SET_VECTOR_ELT(ANS,idx,tmpSEXP); SET_STRING_ELT(ansListNames,idx,mkChar("tN"));
	idx++; PROTECT(tmpSEXP=allocMatrix(REALSXP,(opt->maxKnotNum_Season+1),totalPixelNum));++*nprt;
	matOutput->sNProb=(int32_t *)REAL(tmpSEXP); SET_VECTOR_ELT(ANS,idx,tmpSEXP); SET_STRING_ELT(ansListNames,idx,mkChar("sNProb"));
	idx++; PROTECT(tmpSEXP=allocMatrix(REALSXP,(opt->maxKnotNum_Trend+1),totalPixelNum));++*nprt;
	matOutput->tNProb=(int32_t *)REAL(tmpSEXP); SET_VECTOR_ELT(ANS,idx,tmpSEXP); SET_STRING_ELT(ansListNames,idx,mkChar("tNProb"));
	idx++; PROTECT(tmpSEXP=allocMatrix(REALSXP,N,totalPixelNum));++*nprt;
	matOutput->sProb=(int32_t *)REAL(tmpSEXP); SET_VECTOR_ELT(ANS,idx,tmpSEXP); SET_STRING_ELT(ansListNames,idx,mkChar("sProb"));
	idx++; PROTECT(tmpSEXP=allocMatrix(REALSXP,N,totalPixelNum));++*nprt;
	matOutput->tProb=(int32_t *)REAL(tmpSEXP); SET_VECTOR_ELT(ANS,idx,tmpSEXP); SET_STRING_ELT(ansListNames,idx,mkChar("tProb"));
	idx++; PROTECT(tmpSEXP=allocMatrix(REALSXP,N,totalPixelNum));++*nprt; 
	matOutput->s=(float *)REAL(tmpSEXP); SET_VECTOR_ELT(ANS,idx,tmpSEXP); SET_STRING_ELT(ansListNames,idx,mkChar("s"));
	if (opt->computeCredible) {
		idx++; PROTECT(tmpSEXP=alloc3DArray(REALSXP,N,2,totalPixelNum));++*nprt;
		matOutput->sCI=(float *)REAL(tmpSEXP); SET_VECTOR_ELT(ANS,idx,tmpSEXP); SET_STRING_ELT(ansListNames,idx,mkChar("sCI"));
	}
	idx++; PROTECT(tmpSEXP=allocMatrix(REALSXP,N,totalPixelNum));++*nprt; 
	matOutput->sSD=(float *)REAL(tmpSEXP); SET_VECTOR_ELT(ANS,idx,tmpSEXP); SET_STRING_ELT(ansListNames,idx,mkChar("sSD"));
	idx++; PROTECT(tmpSEXP=allocMatrix(REALSXP,N,totalPixelNum));++*nprt;
	matOutput->t=(float *)REAL(tmpSEXP); SET_VECTOR_ELT(ANS,idx,tmpSEXP); SET_STRING_ELT(ansListNames,idx,mkChar("t"));
	if (opt->computeCredible) {
		idx++; PROTECT(tmpSEXP=alloc3DArray(REALSXP,N,2,totalPixelNum));++*nprt;
		matOutput->tCI=(float *)REAL(tmpSEXP); SET_VECTOR_ELT(ANS,idx,tmpSEXP); SET_STRING_ELT(ansListNames,idx,mkChar("tCI"));
	}
	idx++; PROTECT(tmpSEXP=allocMatrix(REALSXP,N,totalPixelNum));++*nprt; 
	matOutput->tSD=(float *)REAL(tmpSEXP); SET_VECTOR_ELT(ANS,idx,tmpSEXP); SET_STRING_ELT(ansListNames,idx,mkChar("tSD"));
	idx++; PROTECT(tmpSEXP=allocMatrix(REALSXP,N,totalPixelNum));++*nprt; 
	matOutput->b=(float *)REAL(tmpSEXP); SET_VECTOR_ELT(ANS,idx,tmpSEXP); SET_STRING_ELT(ansListNames,idx,mkChar("b"));
	if (opt->computeCredible) {
		idx++; PROTECT(tmpSEXP=alloc3DArray(REALSXP,N,2,totalPixelNum));++*nprt;
		matOutput->bCI=(float *)REAL(tmpSEXP); SET_VECTOR_ELT(ANS,idx,tmpSEXP); SET_STRING_ELT(ansListNames,idx,mkChar("bCI"));
	}
	idx++; PROTECT(tmpSEXP=allocMatrix(REALSXP,N,totalPixelNum));++*nprt; 
	matOutput->bSD=(float *)REAL(tmpSEXP); SET_VECTOR_ELT(ANS,idx,tmpSEXP); SET_STRING_ELT(ansListNames,idx,mkChar("bSD"));
	idx++; PROTECT(tmpSEXP=allocVector(REALSXP,totalPixelNum));++*nprt; 
	matOutput->marg_lik=(float *)REAL(tmpSEXP); SET_VECTOR_ELT(ANS,idx,tmpSEXP); SET_STRING_ELT(ansListNames,idx,mkChar("marg_lik"));
	idx++; PROTECT(tmpSEXP=allocVector(REALSXP,totalPixelNum));++*nprt; 
	matOutput->sig2=(float *)REAL(tmpSEXP); SET_VECTOR_ELT(ANS,idx,tmpSEXP); SET_STRING_ELT(ansListNames,idx,mkChar("sig2"));
	if (opt->computeSlopeSign)
	{
		idx++; PROTECT(tmpSEXP=allocMatrix(REALSXP,N,totalPixelNum));++*nprt;
		matOutput->bsign=(I32PTR )REAL(tmpSEXP); SET_VECTOR_ELT(ANS,idx,tmpSEXP); SET_STRING_ELT(ansListNames,idx,mkChar("bsign"));
	}
	if (opt->computeHarmonicOrder)
	{
		idx++; PROTECT(tmpSEXP=allocMatrix(REALSXP,N,totalPixelNum));++*nprt; 
		matOutput->horder=(I32PTR)REAL(tmpSEXP); SET_VECTOR_ELT(ANS,idx,tmpSEXP); SET_STRING_ELT(ansListNames,idx,mkChar("horder"));
	}
	if (opt->computeTrendOrder)
	{
		idx++; PROTECT(tmpSEXP=allocMatrix(REALSXP,N,totalPixelNum));++*nprt;
		matOutput->torder=(I32PTR) REAL(tmpSEXP); SET_VECTOR_ELT(ANS,idx,tmpSEXP); SET_STRING_ELT(ansListNames,idx,mkChar("torder"));
	}
	if (opt->computeChangepoints)
	{
		idx++; PROTECT(tmpSEXP=allocMatrix(REALSXP,opt->maxKnotNum_Season,totalPixelNum));++*nprt; 
		matOutput->scp=(float *)REAL(tmpSEXP); SET_VECTOR_ELT(ANS,idx,tmpSEXP); SET_STRING_ELT(ansListNames,idx,mkChar("scp"));
		idx++; PROTECT(tmpSEXP=allocMatrix(REALSXP,opt->maxKnotNum_Trend,totalPixelNum));
		++*nprt;  matOutput->tcp=(float *)REAL(tmpSEXP); SET_VECTOR_ELT(ANS,idx,tmpSEXP); SET_STRING_ELT(ansListNames,idx,mkChar("tcp"));
		idx++; PROTECT(tmpSEXP=alloc3DArray(REALSXP,opt->maxKnotNum_Season,2,totalPixelNum));++*nprt; 
		matOutput->scpCI=(float *)REAL(tmpSEXP); SET_VECTOR_ELT(ANS,idx,tmpSEXP); SET_STRING_ELT(ansListNames,idx,mkChar("scpCI"));
		idx++; PROTECT(tmpSEXP=alloc3DArray(REALSXP,opt->maxKnotNum_Trend,2,totalPixelNum)); 
		++*nprt;  matOutput->tcpCI=(float *)REAL(tmpSEXP); SET_VECTOR_ELT(ANS,idx,tmpSEXP); SET_STRING_ELT(ansListNames,idx,mkChar("tcpCI"));
	}
	setAttrib(ANS,R_NamesSymbol,ansListNames);
	PROTECT(tmpSEXP=mkString("beast"));
	setAttrib(ANS,R_ClassSymbol,tmpSEXP);
	UNPROTECT(1);
	PROTECT(tmpSEXP=mkString("beastST"));
	setAttrib(ANS,install("algorithm"),tmpSEXP);
	UNPROTECT(1);
	double *doublePointer=(double*)matOutput->time;
	for (ptrdiff_t i=(ptrdiff_t)(N - 1); i >=0; i--)
		*(doublePointer+i)=(double)*(matOutput->time+i);
	return  ANS;
}
bool get_options_trend(Options *  _restrict opt,SEXP pY,SEXP pOpt,MemPointers * _restrict MEM)
{
	SEXP	tmpSEXP=NULL;
	void	*tmpPointer;
	float * R_INPUTDATA=NULL;
	int nprt=0;
	int N=-1;
	if (!isNewList(pOpt))
	{
		MEM->free_all(MEM);
		r_error("Error: The input parameter OPTION should be a List variable!\n");
	}
	bool isPeriodParGiven=false;
	opt->period=-9999; 
	tmpSEXP=getListElement(pOpt,"omissionValue");
	opt->omissionValue=(tmpSEXP !=NULL) ? (PROTECT(tmpSEXP=coerceVector(tmpSEXP,REALSXP)),nprt++,(float)asReal(tmpSEXP)) : -9999;
	if (TYPEOF(pY)==STRSXP)
	{	
		opt->inputFromDisk=true;
		char const *tmp=CHAR(STRING_ELT(pY,0));
		opt->inputFile=MEM->alloc(MEM,sizeof(char)*(strlen(tmp)+1),0);
		strcpy(opt->inputFile,tmp);
		tmpSEXP=getListElement(pOpt,"lengthPerTimeSeries_infile");
		if (tmpSEXP==NULL)
		{
			MEM->free_all(MEM);
			r_error("Error: If data ara read from a binary file,the 'opt$lengthPerTimeSeries_infile' must be set"
				" to tell how many sample points a time series has (i.e.,the length of time series.");
		}
		PROTECT(tmpSEXP=coerceVector(tmpSEXP,INTSXP));  nprt++;
		N=INTEGER_VALUE(tmpSEXP);
		tmpPointer=fopen(opt->inputFile,"rb");
		if (isPeriodParGiven==false)
		{
			float *data=(float *)r_malloc(sizeof(float)*N);
			rI32 readBytes=fread(data,sizeof(float),N,(FILE *)tmpPointer);
			opt->period=determine_period(data,N,opt->omissionValue);
			r_free(data);
			r_printf("THE PERIOD PARAMETER IS NOT SET. A BEST GUESS of IT IS%d AND USED IN THE DECOMPOSITION.\n",(int)opt->period);
		}
		fseek((FILE *)tmpPointer,0,SEEK_END);
		int64_t fileSizeinByte=ftell((FILE *)tmpPointer);
		fclose((FILE *)tmpPointer);
		opt->totalPixelNum=(uint32_t)(fileSizeinByte/(4 * N));
		r_printf("File Size:%d ; Number of time series/pixels:%d \n",fileSizeinByte,opt->totalPixelNum);
	}
	else if (TYPEOF(pY)==REALSXP)
	{ 
		opt->inputFromDisk=false;
		opt->yInputData=(float *) REAL(pY);
		if (isMatrix(pY))
		{
			SEXP dims=getAttrib(pY,R_DimSymbol);
			N=(int)INTEGER(dims)[0];
			opt->totalPixelNum=(int)INTEGER(dims)[1];
		}
		else if (isReal(pY))
		{
			N=Rf_length(pY);
			opt->totalPixelNum=1;			
		}
		if (opt->period <0) 
		{
			float *data=(float *)r_malloc(sizeof(float)*N);
			{
				double *tmpDouble=(double *)opt->yInputData;
				for (int i=1; i <=N; i++)
					*data++=(float)(*tmpDouble++);
				data=data - N;
			}
			opt->period=determine_period(data,N,opt->omissionValue);
			r_printf("THE PERIOD PARAMETER IS NOT SET. A BEST GUESS of IT IS%d AND USED IN THE DECOMPOSITION.\n",(int)opt->period);
			r_free(data);
		}
	}
	opt->isSingleyInput=0;
	tmpSEXP=getListElement(pOpt,"startTime");
	opt->startTime=(tmpSEXP !=NULL) ? (PROTECT(tmpSEXP=coerceVector(tmpSEXP,REALSXP)),nprt++,(float)asReal(tmpSEXP)) : 1.f;
	tmpSEXP=getListElement(pOpt,"timeInterval");
	opt->timeInterval=(tmpSEXP !=NULL) ? (PROTECT(tmpSEXP=coerceVector(tmpSEXP,REALSXP)),nprt++,(float)asReal(tmpSEXP)) : 1.f;
	tmpSEXP=getListElement(pOpt,"minTrendOrder");
	opt->minTrendOrder=(tmpSEXP !=NULL) ? (PROTECT(tmpSEXP=coerceVector(tmpSEXP,INTSXP)),nprt++,(unsigned char)asInteger(tmpSEXP)) : 0;
	tmpSEXP=getListElement(pOpt,"maxTrendOrder");
	opt->maxTrendOrder=(tmpSEXP !=NULL) ? (PROTECT(tmpSEXP=coerceVector(tmpSEXP,INTSXP)),nprt++,(unsigned char)asInteger(tmpSEXP)) : 1;
	opt->maxTrendOrder=max(opt->maxTrendOrder,opt->minTrendOrder);
	tmpSEXP=getListElement(pOpt,"minSepDist_Trend");
	opt->minSepDist_Trend=(tmpSEXP !=NULL) ? (PROTECT(tmpSEXP=coerceVector(tmpSEXP,INTSXP)),nprt++,(uint16_t)asInteger(tmpSEXP)) : (uint16_t)(opt->period/2);
	tmpSEXP=getListElement(pOpt,"maxKnotNum_Trend");
	opt->maxKnotNum_Trend=(tmpSEXP !=NULL) ? (PROTECT(tmpSEXP=coerceVector(tmpSEXP,INTSXP)),nprt++,(uint16_t)asInteger(tmpSEXP)) : (uint16_t)floor(N/(opt->minSepDist_Trend+1) - 1);
	opt->maxKnotNum_Trend=(uint16_t)min(floor((float)N/(opt->minSepDist_Trend+1) - 1),opt->maxKnotNum_Trend);
	tmpSEXP=getListElement(pOpt,"maxMoveStepSize");
	opt->maxMoveStepSize=(tmpSEXP !=NULL) ? (PROTECT(tmpSEXP=coerceVector(tmpSEXP,INTSXP)),nprt++,(uint16_t)asInteger(tmpSEXP)) : (uint16_t)(opt->minSepDist_Trend/3+1);
	tmpSEXP=getListElement(pOpt,"samples");
	opt->samples=(tmpSEXP !=NULL) ? (PROTECT(tmpSEXP=coerceVector(tmpSEXP,INTSXP)),nprt++,(mwSize)asInteger(tmpSEXP)) : 10000;
	tmpSEXP=getListElement(pOpt,"thinningFactor");
	opt->thinningFactor=(tmpSEXP !=NULL) ? (PROTECT(tmpSEXP=coerceVector(tmpSEXP,INTSXP)),nprt++,(uint16_t)asInteger(tmpSEXP)) : 1;
	tmpSEXP=getListElement(pOpt,"burnin");
	opt->burnin=(tmpSEXP !=NULL) ? (PROTECT(tmpSEXP=coerceVector(tmpSEXP,INTSXP)),nprt++,(mwSize)asInteger(tmpSEXP)) : 500;
	opt->burnin=max(opt->burnin,500);
	tmpSEXP=getListElement(pOpt,"chainNumber");
	opt->chainNumber=(tmpSEXP !=NULL) ? (PROTECT(tmpSEXP=coerceVector(tmpSEXP,INTSXP)),nprt++,(mwSize)asInteger(tmpSEXP)) : 5;
	opt->chainNumber=max(opt->chainNumber,1);
	tmpSEXP=getListElement(pOpt,"resamplingTrendOrderProb");
	opt->resamplingTrendOrderProb=(tmpSEXP !=NULL) ? (PROTECT(tmpSEXP=coerceVector(tmpSEXP,REALSXP)),nprt++,(float)asReal(tmpSEXP)) : 0.1f;
	tmpSEXP=getListElement(pOpt,"resamplingSeasonOrderProb");
	opt->resamplingSeasonOrderProb=(tmpSEXP !=NULL) ? (PROTECT(tmpSEXP=coerceVector(tmpSEXP,REALSXP)),nprt++,(float)asReal(tmpSEXP)) : 0.17f;
	tmpSEXP=getListElement(pOpt,"seed");
	opt->seed=(tmpSEXP !=NULL) ? (PROTECT(tmpSEXP=coerceVector(tmpSEXP,INTSXP)),nprt++,(uint64_t)asInteger(tmpSEXP)) : 0;
	tmpSEXP=getListElement(pOpt,"outputToDisk");
	opt->outputToDisk=(tmpSEXP !=NULL) ? (PROTECT(tmpSEXP=coerceVector(tmpSEXP,INTSXP)),nprt++,(bool)asInteger(tmpSEXP)) : false;
	if (opt->outputToDisk)
	{
		tmpSEXP=getListElement(pOpt,"outputFolder");
		if (tmpSEXP !=NULL)
		{
			char const *tmpChar=CHAR(asChar(tmpSEXP));
			opt->outputFolder=MEM->alloc(MEM,sizeof(char)*(strlen(tmpChar)+20),0);
			strcpy(opt->outputFolder,tmpChar);
		}
		else
			opt->outputToDisk=false;
	}
	if (!opt->outputToDisk && nlhs==0)
	{
		MEM->free_all(MEM);
		r_error("To run,you must save the output to either a lefthand-side"
			"variable (e.g.,out=func(...) ) or specify the output path \n"
			"in the Option.outputFolder parameter. \n");
		return 0;
	}
	tmpSEXP=getListElement(pOpt,"printToScreen");
	opt->printToScreen=(tmpSEXP !=NULL) ? (PROTECT(tmpSEXP=coerceVector(tmpSEXP,INTSXP)),nprt++,(bool)asInteger(tmpSEXP)) : false;
	if (opt->printToScreen)
	{
		tmpSEXP=getListElement(pOpt,"printCharLen");
		opt->printCharLen=(tmpSEXP !=NULL) ? (PROTECT(tmpSEXP=coerceVector(tmpSEXP,INTSXP)),nprt++,(int16_t)asInteger(tmpSEXP)) : 80;
	}
	tmpSEXP=getListElement(pOpt,"computeCredible");
	opt->computeCredible=(tmpSEXP !=NULL) ? (PROTECT(tmpSEXP=coerceVector(tmpSEXP,INTSXP)),nprt++,(bool)asInteger(tmpSEXP)) : false;
	if (opt->computeCredible)
	{
		tmpSEXP=getListElement(pOpt,"fastCIComputation");
		opt->fastCIComputation=(tmpSEXP !=NULL) ? (PROTECT(tmpSEXP=coerceVector(tmpSEXP,INTSXP)),nprt++,(bool)asInteger(tmpSEXP)) : true;
		tmpSEXP=getListElement(pOpt,"alphaLevel");
		opt->alphaLevel=(tmpSEXP !=NULL) ? (PROTECT(tmpSEXP=coerceVector(tmpSEXP,INTSXP)),nprt++,(float)asReal(tmpSEXP)) : 0.95;
	}
	tmpSEXP=getListElement(pOpt,"computeSlopeSign");
	opt->computeSlopeSign=(tmpSEXP !=NULL) ? (PROTECT(tmpSEXP=coerceVector(tmpSEXP,INTSXP)),nprt++,(bool)asInteger(tmpSEXP)) : false;
	tmpSEXP=getListElement(pOpt,"computeTrendOrder");
	opt->computeTrendOrder=(tmpSEXP !=NULL) ? (PROTECT(tmpSEXP=coerceVector(tmpSEXP,INTSXP)),nprt++,(bool)asInteger(tmpSEXP)) : false;
	tmpSEXP=getListElement(pOpt,"computeChangepoints");
	opt->computeChangepoints=(tmpSEXP !=NULL) ? (PROTECT(tmpSEXP=coerceVector(tmpSEXP,INTSXP)),nprt++,(bool)asInteger(tmpSEXP)) : false;
	tmpSEXP=getListElement(pOpt,"ridgeFactor");
	opt->ridgeFactor=(tmpSEXP !=NULL) ? (PROTECT(tmpSEXP=coerceVector(tmpSEXP,REALSXP)),nprt++,(float)asReal(tmpSEXP)) : 0.0001f;
	UNPROTECT(nprt);
	nprt=0;
	opt->Npad=(int32_t)ceil((float)N/8.0f) * 8;
	opt->Npad16=(int32_t)ceil((float)N/16.0f) * 16;
	opt->N=N;
	opt->separator='$';
	return true;
}
SEXP allocate_output_trend(RESULT * _restrict matOutput,Options * _restrict opt,int * _restrict nprt)
{
	SEXP ANS=NULL;
	SEXP ansListNames=NULL;
	int N=opt->N;
	int totalPixelNum=opt->totalPixelNum;
	SEXP tmpSEXP;
	int idx;
	int32_t elementNum=10+(int32_t)opt->computeCredible * 2+\
		(int32_t)opt->computeSlopeSign+(int32_t)opt->computeTrendOrder+(int32_t)opt->computeChangepoints * 2;
	PROTECT(ANS=allocVector(VECSXP,elementNum));++*nprt;
	PROTECT(ansListNames=allocVector(STRSXP,elementNum));++*nprt;
	idx=0; PROTECT(tmpSEXP=allocVector(REALSXP,N));++*nprt; 
	matOutput->time=(float *)REAL(tmpSEXP); SET_VECTOR_ELT(ANS,idx,tmpSEXP); SET_STRING_ELT(ansListNames,idx,mkChar("time"));
	idx++; PROTECT(tmpSEXP=allocVector(REALSXP,totalPixelNum));++*nprt; 
	matOutput->tN=(float *)REAL(tmpSEXP); SET_VECTOR_ELT(ANS,idx,tmpSEXP); SET_STRING_ELT(ansListNames,idx,mkChar("tN"));
	idx++; PROTECT(tmpSEXP=allocMatrix(REALSXP,(opt->maxKnotNum_Trend+1),totalPixelNum));++*nprt; 
	matOutput->tNProb=(int32_t *)REAL(tmpSEXP); SET_VECTOR_ELT(ANS,idx,tmpSEXP); SET_STRING_ELT(ansListNames,idx,mkChar("tNProb"));
	idx++; PROTECT(tmpSEXP=allocMatrix(REALSXP,N,totalPixelNum));++*nprt;
	matOutput->tProb=(int32_t *)REAL(tmpSEXP); SET_VECTOR_ELT(ANS,idx,tmpSEXP); SET_STRING_ELT(ansListNames,idx,mkChar("tProb"));
	idx++; PROTECT(tmpSEXP=allocMatrix(REALSXP,N,totalPixelNum));++*nprt; 
	matOutput->t=(float *)REAL(tmpSEXP); SET_VECTOR_ELT(ANS,idx,tmpSEXP); SET_STRING_ELT(ansListNames,idx,mkChar("t"));
	if (opt->computeCredible){
		idx++; PROTECT(tmpSEXP=alloc3DArray(REALSXP,N,2,totalPixelNum));++*nprt;
		matOutput->tCI=(float *)REAL(tmpSEXP); SET_VECTOR_ELT(ANS,idx,tmpSEXP); SET_STRING_ELT(ansListNames,idx,mkChar("tCI"));
	}
	idx++; PROTECT(tmpSEXP=allocMatrix(REALSXP,N,totalPixelNum));++*nprt; 
	matOutput->tSD=(float *)REAL(tmpSEXP); SET_VECTOR_ELT(ANS,idx,tmpSEXP); SET_STRING_ELT(ansListNames,idx,mkChar("tSD"));
	idx++; PROTECT(tmpSEXP=allocMatrix(REALSXP,N,totalPixelNum));++*nprt; 
	matOutput->b=(float *)REAL(tmpSEXP); SET_VECTOR_ELT(ANS,idx,tmpSEXP); SET_STRING_ELT(ansListNames,idx,mkChar("b"));
	if (opt->computeCredible){
		idx++; PROTECT(tmpSEXP=alloc3DArray(REALSXP,N,2,totalPixelNum));++*nprt;
		matOutput->bCI=(float *)REAL(tmpSEXP); SET_VECTOR_ELT(ANS,idx,tmpSEXP); SET_STRING_ELT(ansListNames,idx,mkChar("bCI"));
	}
	idx++; PROTECT(tmpSEXP=allocMatrix(REALSXP,N,totalPixelNum));++*nprt; 
	matOutput->bSD=(float *)REAL(tmpSEXP); SET_VECTOR_ELT(ANS,idx,tmpSEXP); SET_STRING_ELT(ansListNames,idx,mkChar("bSD"));
	idx++; PROTECT(tmpSEXP=allocVector(REALSXP,totalPixelNum));++*nprt; 
	matOutput->marg_lik=(float *)REAL(tmpSEXP); SET_VECTOR_ELT(ANS,idx,tmpSEXP); SET_STRING_ELT(ansListNames,idx,mkChar("marg_lik"));
	idx++; PROTECT(tmpSEXP=allocVector(REALSXP,totalPixelNum));++*nprt; 
	matOutput->sig2=(float *)REAL(tmpSEXP); SET_VECTOR_ELT(ANS,idx,tmpSEXP); SET_STRING_ELT(ansListNames,idx,mkChar("sig2"));
	if (opt->computeSlopeSign)
	{
		idx++; PROTECT(tmpSEXP=allocMatrix(REALSXP,N,totalPixelNum));++*nprt; 
		matOutput->bsign=(int32_t *)REAL(tmpSEXP); SET_VECTOR_ELT(ANS,idx,tmpSEXP); SET_STRING_ELT(ansListNames,idx,mkChar("bsign"));
	}
	if (opt->computeTrendOrder)
	{
		idx++; PROTECT(tmpSEXP=allocMatrix(REALSXP,N,totalPixelNum));++*nprt;
		matOutput->torder=(int32_t *)REAL(tmpSEXP); SET_VECTOR_ELT(ANS,idx,tmpSEXP); SET_STRING_ELT(ansListNames,idx,mkChar("toder"));
	}
	if (opt->computeChangepoints)
	{
		idx++; PROTECT(tmpSEXP=allocMatrix(REALSXP,opt->maxKnotNum_Trend,totalPixelNum));++*nprt;  matOutput->tcp=(float *)REAL(tmpSEXP); SET_VECTOR_ELT(ANS,idx,tmpSEXP); SET_STRING_ELT(ansListNames,idx,mkChar("tcp"));
		idx++; PROTECT(tmpSEXP=alloc3DArray(REALSXP,opt->maxKnotNum_Trend,2,totalPixelNum));++*nprt;  matOutput->tcpCI=(float *)REAL(tmpSEXP); SET_VECTOR_ELT(ANS,idx,tmpSEXP); SET_STRING_ELT(ansListNames,idx,mkChar("tcpCI"));
	}
	setAttrib(ANS,R_NamesSymbol,ansListNames);
	PROTECT(tmpSEXP=mkString("beast"));
	setAttrib(ANS,R_ClassSymbol,tmpSEXP);
	UNPROTECT(1);
	PROTECT(tmpSEXP=mkString("beastST"));
	setAttrib(ANS,install("algorithm"),tmpSEXP);
	UNPROTECT(1);
	double *doublePointer=(double*)matOutput->time;
	for (ptrdiff_t i=(ptrdiff_t)(N - 1); i >=0; i--)
		*(doublePointer+i)=(double)*(matOutput->time+i);
	return  ANS;
}
#endif
ENABLE_MANY_WARNINGS
