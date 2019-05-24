#include "abc_001_config.h"
static char fileID='R';
#if R_INTERFACE==1
#ifdef __GNU__
#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <X11/Xresource.h>
#include <X11/Xlocale.h>
#include <X11/Xatom.h>
#endif
#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <R_ext/Rdynload.h>
#define CALLDEF(name,n) {#name,(DL_FUNC) &name,n}
extern  SEXP DllExport LAI_estimate(SEXP pY1,SEXP pY2,SEXP pOpt);
extern  SEXP DllExport LAI_gfunc(SEXP pY1,SEXP pY2,SEXP pOpt);
extern  SEXP DllExport LAI_simulate(SEXP pY1,SEXP pY2,SEXP pOpt);
#if defined(_WIN32) && defined(_WIN64) && !defined(__i386)  && !defined(__i686) && !defined(i386) && !defined(__i686)
#endif
static const R_CallMethodDef CallEntries[]={
			CALLDEF(LAI_estimate,3),
			CALLDEF(LAI_gfunc,3),
			CALLDEF(LAI_simulate,3),
#if defined(_WIN32) && defined(_WIN64) && !defined(__i386)  && !defined(__i686) && !defined(i386) && !defined(__i686)
#endif
			{ NULL,NULL,0 } 
};
void  R_init_Rbeast(DllInfo *dll)
{
	R_registerRoutines(dll,NULL,CallEntries,NULL,NULL);
	R_useDynamicSymbols(dll,FALSE);
}
#endif
