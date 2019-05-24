#include "abc_000_macro.h"
#include "hemiphoto.h"
#if defined(__GNUC__)||defined(__CLANG__) 
DISABLE_WARNING(unused-variable,unused-variable,NOT_USED)
DISABLE_WARNING(unused-function,unused-function,NOT_USED)
#endif
char fileID='a';
#if R_INTERFACE==1
#include "hemiphoto_lib_R.h"
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
bool get_options(OPT * _restrict opt,SEXP pY1,SEXP pY2,SEXP pOpt,MemPointers * _restrict MEM)
{
	SEXP	tmpSEXP=NULL;
	void	*tmpPointer;
	int nprt=0;
	opt->THETA=REAL(pY1);
	opt->GAP=REAL(pY2);
	if (isMatrix(pY1))
	{
		SEXP dims=getAttrib(pY1,R_DimSymbol);
		opt->NPIXEL=(INTEGER(dims)[0] )  * ( INTEGER(dims)[1]);
	}
	else if (isReal(pY2))
	{
		opt->NPIXEL=Rf_length(pY1);
	}
	else
	{
		MEM->free_all(MEM);
		r_error("Error: The input parameter OPTION should be a List variable!\n");
	}
	if (!isNewList(pOpt))
	{
		MEM->free_all(MEM);
		r_error("Error: The input parameter OPTION should be a List variable!\n");
	}
	opt->ite=-100;
	opt->NBIN=200;
	opt->Nfrac=8;
	opt->GQ_knotNumber=21;
	tmpSEXP=getListElement(pOpt,"ite");
	opt->ite=(tmpSEXP !=NULL) ? (PROTECT(tmpSEXP=coerceVector(tmpSEXP,REALSXP)),nprt++,(int)asReal(tmpSEXP)) : -100;
	tmpSEXP=getListElement(pOpt,"nbin");
	opt->NBIN=(tmpSEXP !=NULL) ? (PROTECT(tmpSEXP=coerceVector(tmpSEXP,REALSXP)),nprt++,(int)asReal(tmpSEXP)) : 200;
	tmpSEXP=getListElement(pOpt,"nfrac");
	opt->Nfrac=(tmpSEXP !=NULL) ? (PROTECT(tmpSEXP=coerceVector(tmpSEXP,INTSXP)),nprt++,(int)asReal(tmpSEXP)) : 8;
	tmpSEXP=getListElement(pOpt,"gq_knotNum");
	opt->GQ_knotNumber=(tmpSEXP !=NULL) ? (PROTECT(tmpSEXP=coerceVector(tmpSEXP,INTSXP)),nprt++,(int)asReal(tmpSEXP)) : 21;
	if (isString(pOpt))
	{
		char  const *tmp=CHAR(STRING_ELT(pOpt,0));
		if (strcicmp(tmp,"mle")==0)
			opt->algorithm=1;
		else
			opt->algorithm=2;
	}
	UNPROTECT(nprt);
	nprt=0;
	return true;
}
SEXP allocate_output( RES * _restrict result,int N,int * _restrict nprt,char **gfuncNames )
{
	SEXP ANS=NULL;
	SEXP ansListNames=NULL; 
	SEXP tmpSEXP;
	int idx;
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
	int32_t elementNum=sizeof(names)/sizeof(names[0]);
	PROTECT(ANS=allocVector(VECSXP,elementNum));++*nprt;
	PROTECT(ansListNames=allocVector(STRSXP,elementNum));++*nprt;
	idx=0; PROTECT(tmpSEXP=allocVector(STRSXP,N));++*nprt;
	for (rINT i=0; i < N; i++)
	{
		SET_STRING_ELT(tmpSEXP,i,mkChar(gfuncNames[i]) );
	}
	result->LAD_TYPE=tmpSEXP; SET_VECTOR_ELT(ANS,idx,tmpSEXP); SET_STRING_ELT(ansListNames,idx,mkChar(names[idx]));
	idx++; PROTECT(tmpSEXP=allocVector(REALSXP,N));++*nprt;
	result->LAI_HA57=REAL(tmpSEXP); SET_VECTOR_ELT(ANS,idx,tmpSEXP); SET_STRING_ELT(ansListNames,idx,mkChar(names[idx]));
	idx++; PROTECT(tmpSEXP=allocVector(REALSXP,N));++*nprt;
	result->LAI_MILLER=REAL(tmpSEXP); SET_VECTOR_ELT(ANS,idx,tmpSEXP); SET_STRING_ELT(ansListNames,idx,mkChar(names[idx]));
	idx++; PROTECT(tmpSEXP=allocVector(REALSXP,N));++*nprt;
	result->LAI_BNR=REAL(tmpSEXP); SET_VECTOR_ELT(ANS,idx,tmpSEXP); SET_STRING_ELT(ansListNames,idx,mkChar(names[idx]));
	idx++; PROTECT(tmpSEXP=allocVector(REALSXP,N));++*nprt; 
	result->AIC_BNR=REAL(tmpSEXP); SET_VECTOR_ELT(ANS,idx,tmpSEXP); SET_STRING_ELT(ansListNames,idx,mkChar(names[idx]));
	idx++; PROTECT(tmpSEXP=allocVector(REALSXP,N));++*nprt;
	result->LAD1_BNR=REAL(tmpSEXP); SET_VECTOR_ELT(ANS,idx,tmpSEXP); SET_STRING_ELT(ansListNames,idx,mkChar(names[idx]));
	idx++; PROTECT(tmpSEXP=allocVector(REALSXP,N));++*nprt;
	result->LAD2_BNR=REAL(tmpSEXP); SET_VECTOR_ELT(ANS,idx,tmpSEXP); SET_STRING_ELT(ansListNames,idx,mkChar(names[idx]));
	idx++; PROTECT(tmpSEXP=allocVector(REALSXP,N));++*nprt;
	result->LAI_L2F=REAL(tmpSEXP); SET_VECTOR_ELT(ANS,idx,tmpSEXP); SET_STRING_ELT(ansListNames,idx,mkChar(names[idx]));
	idx++; PROTECT(tmpSEXP=allocVector(REALSXP,N));++*nprt;
	result->AIC_L2F=REAL(tmpSEXP); SET_VECTOR_ELT(ANS,idx,tmpSEXP); SET_STRING_ELT(ansListNames,idx,mkChar(names[idx]));
	idx++; PROTECT(tmpSEXP=allocVector(REALSXP,N));++*nprt;
	result->LAD1_L2F=REAL(tmpSEXP); SET_VECTOR_ELT(ANS,idx,tmpSEXP); SET_STRING_ELT(ansListNames,idx,mkChar(names[idx]));
	idx++; PROTECT(tmpSEXP=allocVector(REALSXP,N));++*nprt;
	result->LAD2_L2F=REAL(tmpSEXP); SET_VECTOR_ELT(ANS,idx,tmpSEXP); SET_STRING_ELT(ansListNames,idx,mkChar(names[idx]));
	idx++; PROTECT(tmpSEXP=allocVector(REALSXP,N));++*nprt;
	result->LAI_L1F=REAL(tmpSEXP); SET_VECTOR_ELT(ANS,idx,tmpSEXP); SET_STRING_ELT(ansListNames,idx,mkChar(names[idx]));
	idx++; PROTECT(tmpSEXP=allocVector(REALSXP,N));++*nprt;
	result->AIC_L1F=REAL(tmpSEXP); SET_VECTOR_ELT(ANS,idx,tmpSEXP); SET_STRING_ELT(ansListNames,idx,mkChar(names[idx]));
	idx++; PROTECT(tmpSEXP=allocVector(REALSXP,N));++*nprt;
	result->LAD1_L1F=REAL(tmpSEXP); SET_VECTOR_ELT(ANS,idx,tmpSEXP); SET_STRING_ELT(ansListNames,idx,mkChar(names[idx]));
	idx++; PROTECT(tmpSEXP=allocVector(REALSXP,N));++*nprt;
	result->LAD2_L1F=REAL(tmpSEXP); SET_VECTOR_ELT(ANS,idx,tmpSEXP); SET_STRING_ELT(ansListNames,idx,mkChar(names[idx]));
	idx++; PROTECT(tmpSEXP=allocVector(REALSXP,N));++*nprt;
	result->LAI_L2LF=REAL(tmpSEXP); SET_VECTOR_ELT(ANS,idx,tmpSEXP); SET_STRING_ELT(ansListNames,idx,mkChar(names[idx]));
	idx++; PROTECT(tmpSEXP=allocVector(REALSXP,N));++*nprt;
	result->AIC_L2LF=REAL(tmpSEXP); SET_VECTOR_ELT(ANS,idx,tmpSEXP); SET_STRING_ELT(ansListNames,idx,mkChar(names[idx]));
	idx++; PROTECT(tmpSEXP=allocVector(REALSXP,N));++*nprt;
	result->LAD1_L2LF=REAL(tmpSEXP); SET_VECTOR_ELT(ANS,idx,tmpSEXP); SET_STRING_ELT(ansListNames,idx,mkChar(names[idx]));
	idx++; PROTECT(tmpSEXP=allocVector(REALSXP,N));++*nprt;
	result->LAD2_L2LF=REAL(tmpSEXP); SET_VECTOR_ELT(ANS,idx,tmpSEXP); SET_STRING_ELT(ansListNames,idx,mkChar(names[idx]));
	idx++; PROTECT(tmpSEXP=allocVector(REALSXP,N));++*nprt;
	result->LAI_L1LF=REAL(tmpSEXP); SET_VECTOR_ELT(ANS,idx,tmpSEXP); SET_STRING_ELT(ansListNames,idx,mkChar(names[idx]));
	idx++; PROTECT(tmpSEXP=allocVector(REALSXP,N));++*nprt;
	result->AIC_L1LF=REAL(tmpSEXP); SET_VECTOR_ELT(ANS,idx,tmpSEXP); SET_STRING_ELT(ansListNames,idx,mkChar(names[idx]));
	idx++; PROTECT(tmpSEXP=allocVector(REALSXP,N));++*nprt;
	result->LAD1_L1LF=REAL(tmpSEXP); SET_VECTOR_ELT(ANS,idx,tmpSEXP); SET_STRING_ELT(ansListNames,idx,mkChar(names[idx]));
	idx++; PROTECT(tmpSEXP=allocVector(REALSXP,N));++*nprt;
	result->LAD2_L1LF=REAL(tmpSEXP); SET_VECTOR_ELT(ANS,idx,tmpSEXP); SET_STRING_ELT(ansListNames,idx,mkChar(names[idx]));
	idx++; PROTECT(tmpSEXP=allocVector(REALSXP,N));++*nprt;
	result->LAI_L2WL=REAL(tmpSEXP); SET_VECTOR_ELT(ANS,idx,tmpSEXP); SET_STRING_ELT(ansListNames,idx,mkChar(names[idx]));
	idx++; PROTECT(tmpSEXP=allocVector(REALSXP,N));++*nprt;
	result->AIC_L2WL=REAL(tmpSEXP); SET_VECTOR_ELT(ANS,idx,tmpSEXP); SET_STRING_ELT(ansListNames,idx,mkChar(names[idx]));
	idx++; PROTECT(tmpSEXP=allocVector(REALSXP,N));++*nprt;
	result->LAD1_L2WL=REAL(tmpSEXP); SET_VECTOR_ELT(ANS,idx,tmpSEXP); SET_STRING_ELT(ansListNames,idx,mkChar(names[idx]));
	idx++; PROTECT(tmpSEXP=allocVector(REALSXP,N));++*nprt;
	result->LAD2_L2WL=REAL(tmpSEXP); SET_VECTOR_ELT(ANS,idx,tmpSEXP); SET_STRING_ELT(ansListNames,idx,mkChar(names[idx]));
	idx++; PROTECT(tmpSEXP=allocVector(REALSXP,N));++*nprt;
	result->LAI_L1WL=REAL(tmpSEXP); SET_VECTOR_ELT(ANS,idx,tmpSEXP); SET_STRING_ELT(ansListNames,idx,mkChar(names[idx]));
	idx++; PROTECT(tmpSEXP=allocVector(REALSXP,N));++*nprt;
	result->AIC_L1WL=REAL(tmpSEXP); SET_VECTOR_ELT(ANS,idx,tmpSEXP); SET_STRING_ELT(ansListNames,idx,mkChar(names[idx]));
	idx++; PROTECT(tmpSEXP=allocVector(REALSXP,N));++*nprt;
	result->LAD1_L1WL=REAL(tmpSEXP); SET_VECTOR_ELT(ANS,idx,tmpSEXP); SET_STRING_ELT(ansListNames,idx,mkChar(names[idx]));
	idx++; PROTECT(tmpSEXP=allocVector(REALSXP,N));++*nprt;
	result->LAD2_L1WL=REAL(tmpSEXP); SET_VECTOR_ELT(ANS,idx,tmpSEXP); SET_STRING_ELT(ansListNames,idx,mkChar(names[idx]));
	setAttrib(ANS,R_NamesSymbol,ansListNames);
	PROTECT(tmpSEXP=mkString("data.frame"));
	setAttrib(ANS,R_ClassSymbol,ScalarString(mkChar("data.frame")));
	UNPROTECT(1);
	PROTECT(tmpSEXP=mkString("hemiphoto"));
	setAttrib(ANS,install("algorithm"),tmpSEXP);
	UNPROTECT(1);
	PROTECT(tmpSEXP=allocVector(STRSXP,N)); 
	for (rINT i=0; i < N; i++)
	{
		char buffer[20];
		itoa(i+1,buffer,10L);   
		buffer[2]=0;
		buffer[3]=0;
		SET_STRING_ELT(tmpSEXP,i,mkChar(buffer));
	}
	setAttrib(ANS,R_RowNamesSymbol,tmpSEXP);
	UNPROTECT(1);
	return  ANS;
}
#endif
#if defined(__GNUC__)||defined(__CLANG__) 
ENABLE_WARNING(unused-variable,unused-variable,NOT_USED)
ENABLE_WARNING(unused-function,unused-function,NOT_USED)
#endif
