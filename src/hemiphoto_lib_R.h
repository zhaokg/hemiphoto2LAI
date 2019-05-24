#pragma once
#include "abc_000_macro.h"
#if R_INTERFACE==1
	#include "abc_001_config.h"
	#include "abc_mem.h"
	#include "abc_common.h"
	#include "hemiphoto.h"
	extern  int nlhs;
	extern  SEXP	getListElement(SEXP list,const char *str);	
	extern bool		get_options(OPT * _restrict opt,SEXP pY1,SEXP pY2,SEXP pOpt,MemPointers * _restrict MEM);	
	extern SEXP allocate_output(RES * _restrict result,int N,int * _restrict nprt,char **names);
#endif
