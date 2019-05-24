
#pragma once
#ifdef _MSC_VER
	#define MSVC 1
	#define GNU 0
#elif defined(__GNUC__)||defined(__clang__)||defined(__APPLE__)||defined(__linux__)||defined(__MINGW32__)||defined(__MINGW64__)||defined(__MACH__)
	#define MSVC 0
	#define GNU 1
#endif
	#define MYMAT_LIBRARY 1
	#define MKL_LIBRARY   0
	#define MATLAB_LIBRARY 0 
	#define R_INTERFACE 1
	#define M_INTERFACE 0
	#define MYRAND_LIBRARY  1
	#define MKLRAND_LIBRARY 0
	#define BASIS_METHODS 2
	#define PTHREAD_INOUT 0
	#define MY_DEBUG 0
	#if MYMAT_LIBRARY==1
		#define MYRAND_LIBRARY 1
		#define MKLRAND_LIBRARY 0
	#endif
	#if MSVC==1
		#define INLINE    __inline
		#define _restrict __restrict
	#elif GNU==1 
		#define INLINE inline
		#define _restrict __restrict__
	#endif
	#define max(a,b)    (((a) > (b)) ? (a) : (b))
	#define min(a,b)    (((a) < (b)) ? (a) : (b))
	#define mv(n,src,dest)	r_cblas_scopy( n,src,1L,dest,1L) 
	#define cp(n,src,dest)    memcpy(dest,src,sizeof(float)*(n))
	#define DIAG_STR(s) #s
	#define DIAG_JOINSTR(x,y) DIAG_STR(x ## y)
	#ifdef _MSC_VER
		#define DIAG_DO_PRAGMA(x) __pragma (#x)
		#define DIAG_PRAGMA(compiler,x) DIAG_DO_PRAGMA(warning(x))
	#else
		#define DIAG_DO_PRAGMA(x) _Pragma (#x)
		#define DIAG_PRAGMA(compiler,x) DIAG_DO_PRAGMA(compiler diagnostic x)
	#endif
	#if defined(__clang__)
		# define DISABLE_WARNING(gcc_unused,clang_option,msvc_unused) DIAG_PRAGMA(clang,push) DIAG_PRAGMA(clang,ignored DIAG_JOINSTR(-W,clang_option))
		# define ENABLE_WARNING(gcc_unused,clang_option,msvc_unused) DIAG_PRAGMA(clang,pop)
	#elif defined(_MSC_VER)
		# define DISABLE_WARNING(gcc_unused,clang_unused,msvc_errorcode) DIAG_PRAGMA(msvc,push) DIAG_DO_PRAGMA(warning(disable:##msvc_errorcode))
		# define ENABLE_WARNING(gcc_unused,clang_unused,msvc_errorcode) DIAG_PRAGMA(msvc,pop)
	#elif defined(__GNUC__)
		#if ((__GNUC__ * 100)+__GNUC_MINOR__) >=406
			# define DISABLE_WARNING(gcc_option,clang_unused,msvc_unused) DIAG_PRAGMA(GCC,push) DIAG_PRAGMA(GCC,ignored DIAG_JOINSTR(-W,gcc_option))
			# define ENABLE_WARNING(gcc_option,clang_unused,msvc_unused) DIAG_PRAGMA(GCC,pop)
		#else
			# define DISABLE_WARNING(gcc_option,clang_unused,msvc_unused) DIAG_PRAGMA(GCC,ignored DIAG_JOINSTR(-W,gcc_option))
			# define ENABLE_WARNING(gcc_option,clang_option,msvc_unused) DIAG_PRAGMA(GCC,warning DIAG_JOINSTR(-W,gcc_option))
		#endif
	#endif
