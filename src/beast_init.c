#include "abc_001_config.h"
DISABLE_MANY_WARNINGS
static char fileID UNUSED_DECORATOR='R';
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
extern SEXP DllExport beastST_multipleChain_fast(SEXP pY,SEXP pOpt);
#if defined(WIN64_OS) 
extern SEXP DllExport WinMainDemoST(SEXP pY,SEXP pOpt);
#endif
static const R_CallMethodDef CallEntries[]={
			CALLDEF(beastST_multipleChain_fast,2),
#if defined(WIN64_OS) 
			CALLDEF(WinMainDemoST,2),
#endif
			{ NULL,NULL,0 } 
};
void  R_init_Rbeast(DllInfo *dll)
{
	R_registerRoutines(dll,NULL,CallEntries,NULL,NULL);
	R_useDynamicSymbols(dll,FALSE);
}
#endif
ENABLE_MANY_WARNINGS
