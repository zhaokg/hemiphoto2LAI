#pragma once
#include "abc_000_macro.h"
#if R_INTERFACE==1
	#include "abc_001_config.h"
	#include "abc_mem.h"
	#include "abc_common.h"
	#include "beast_common.h"
	extern  int nlhs;
	extern  SEXP	getListElement(SEXP list,const char *str);	
	extern  bool    get_options(Options * _restrict opt,SEXP pY,SEXP pOpt,MemPointers * _restrict MEM);
	extern  SEXP    allocate_output(RESULT * _restrict matOutput,Options * _restrict opt,int * _restrict nprt);
	extern  bool    get_options_trend(Options * _restrict opt,SEXP pY,SEXP pOpt,MemPointers * _restrict MEM);
	extern  SEXP    allocate_output_trend(RESULT * _restrict matOutput,Options * _restrict opt,int * _restrict nprt);
#endif
