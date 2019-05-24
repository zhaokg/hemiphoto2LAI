#include "hemiphoto.h"
#include "string.h"
void logitinv(FLOATPTR x,FLOATPTR y,FLOATPTR bnd,int32_t N)
{
	for (rINT i=0; i < N; i++)
	{
		rFLOAT a=*bnd++;
		rFLOAT b=*bnd++;
		rFLOAT xx=*x++;
		*y++=(b - a)*exp(xx)/(1+exp(xx))+a;
	}
}
void logitfwd(FLOATPTR x,FLOATPTR y,FLOATPTR bnd,int32_t N)
{
	for (rINT i=0; i < N; i++)
	{
		rFLOAT a=*bnd++;
		rFLOAT b=*bnd++;
		rFLOAT xx=*x++;
		rFLOAT yy;
		yy=(xx - a)/(b - a);
		yy=log(yy/(1 - yy));
		if (ISNAN(yy)||ISINF(yy))
		{
			yy=log(2 * FLOAT_EPSILON/(1 - 2 * FLOAT_EPSILON));
		}
		*y++=yy;		 
	}
}
FLOAT log_prob_beerlaw(FLOATPTR X,FLOATPTR df,int32_t N,void *par)
{
	FLOATPTR bndL=((PAR_STRUCT*)par)->boundL;
	FLOATPTR bndG=((PAR_STRUCT*)par)->boundG;
	FLOAT  dx=((PAR_STRUCT*)par)->dx;
	FLOATPTR mem=((PAR_STRUCT*)par)->mem;
	FLOAT LAI;
	FLOATPTR gParam=mem+N;   
	logitinv(X,&LAI,bndL,1);
	logitinv(X+1,gParam,bndG,N - 1);
	GFUNC gfun=((PAR_STRUCT*)par)->gfun;
	FLOATPTR G=mem+2*N; 
	FLOATPTR THETA=((PAR_STRUCT*)par)->THETA;
	int32_t  LEN=((PAR_STRUCT*)par)->LEN;
	gfun(THETA,G,LEN,gParam,NULL);
	FLOAT f;
	FLOAT df1;
	f=0;
	df1=0;
	FLOATPTR GAP=((PAR_STRUCT*)par)->GAP;
	for (rINT i=0; i < LEN; i++)
	{
		rFLOAT g=G[i];
		rFLOAT cost=cos(THETA[i]);
		if (GAP[i]==0)
		{
			rFLOAT P0=exp(-LAI*g/cost);
			f+=log(1 -P0);
			df1+=1/(1 - P0)*P0 *g/cost;
		}
		else
		{
			f+=-g * LAI/cost;
			df1+=-g/cost;
		}
	}
	{
		FLOAT LAI1;
		FLOAT X1;
		X1=X[0]+dx;
	    logitinv(&X1,&LAI1,bndL,1);
		df1=df1*(LAI1 - LAI)/dx;
	}
	df[0]=-df1;
	for (rINT i=1; i < N; i++)
	{
		FLOATPTR Xnew=mem;		
		memcpy(Xnew,X,sizeof(FLOAT)*N);
		Xnew[i]=Xnew[i]+dx;
		logitinv(Xnew+1,gParam,bndG,N - 1);
		FLOATPTR Gnew=mem+2*N+LEN;
		gfun(THETA,Gnew,LEN,gParam,NULL);
		rFLOAT DF=0;
		for (rINT j=0; j < LEN; j++)
		{
			rFLOAT g=G[j];
			rFLOAT g1=Gnew[j];
			rFLOAT cost=cos(THETA[j]);
			if (GAP[j]==0)
			{
				rFLOAT P0=exp(-LAI*g/cost);				
				DF+=1/(1 - P0)*P0*LAI*(g1-g)/dx/cost;
			}
			else
			{				
				DF+=-(g1-g)/dx*LAI/cost;
			}
		}
		df[i]=-DF;
	}
	return (FLOAT)(-f);
}
FLOAT alg_BNR(FLOATPTR X,FLOATPTR df,int32_t N,void *par)
{
	FLOATPTR bndL=((PAR_STRUCT*)par)->boundL;
	FLOATPTR bndG=((PAR_STRUCT*)par)->boundG;
	FLOAT  dx=((PAR_STRUCT*)par)->dx;
	FLOATPTR mem=((PAR_STRUCT*)par)->mem;
	FLOAT LAI;
	FLOATPTR gParam=mem+N;   
	logitinv(X,&LAI,bndL,1);
	logitinv(X+1,gParam,bndG,N - 1);
	GFUNC gfun=((PAR_STRUCT*)par)->gfun;
	FLOATPTR G=mem+2 * N; 
	FLOATPTR THETA=((PAR_STRUCT*)par)->THETA01;
	int32_t  LEN=((PAR_STRUCT*)par)->NBIN ;
	gfun(THETA,G,LEN,gParam,NULL);
	FLOAT f;
	FLOAT df1;
	f=0;
	df1=0;
	for (rINT i=0; i < LEN; i++)
	{
		rFLOAT g=G[i];
		rFLOAT cost=cos(THETA[i]);
		{
			rFLOAT P0;
			rFLOAT n=((PAR_STRUCT*)par)->GAP0[i];
			P0=exp(-LAI*g/cost);
			f+=n*log(1 - P0);
			df1+=n*1/(1 - P0)*P0 *g/cost;
		}
		{
			rFLOAT n=((PAR_STRUCT*)par)->GAP1[i];
			f+=n*(-g * LAI/cost);
			df1+=n*(-g/cost);
		}
	}
	{
		FLOAT LAI1;
		FLOAT X1;
		X1=X[0]+dx;
		logitinv(&X1,&LAI1,bndL,1);
		df1=df1*(LAI1 - LAI)/dx;
	}
	df[0]=-df1;
	for (rINT i=1; i < N; i++)
	{
		FLOATPTR Xnew=mem;
		memcpy(Xnew,X,sizeof(FLOAT)*N);
		Xnew[i]=Xnew[i]+dx;
		logitinv(Xnew+1,gParam,bndG,N - 1);
		FLOATPTR Gnew=mem+2 * N+LEN;
		gfun(THETA,Gnew,LEN,gParam,NULL);
		rFLOAT DF=0;
		for (rINT j=0; j < LEN; j++)
		{
			rFLOAT g=G[j];
			rFLOAT g1=Gnew[j];
			rFLOAT cost=cos(THETA[j]);
			{
				rFLOAT n=((PAR_STRUCT*)par)->GAP0[j];
				rFLOAT P0=exp(-LAI*g/cost);				
				DF+=n*( 1/(1 - P0)*P0*LAI*(g1 - g)/dx/cost );
			}
			{
				rFLOAT n=((PAR_STRUCT*)par)->GAP1[j];
				DF+=n*-(g1 - g)/dx*LAI/cost;
			}
		}
		df[i]=-DF;
	}
	return (FLOAT)(-f);
}
FLOAT alg_L2F(FLOATPTR X,FLOATPTR df,int32_t N,void *par)
{
	FLOATPTR bndL=((PAR_STRUCT*)par)->boundL;
	FLOATPTR bndG=((PAR_STRUCT*)par)->boundG;
	FLOAT  dx=((PAR_STRUCT*)par)->dx;
	FLOATPTR mem=((PAR_STRUCT*)par)->mem;
	FLOAT LAI;
	FLOATPTR gParam=mem+N;   
	logitinv(X,&LAI,bndL,1);
	logitinv(X+1,gParam,bndG,N - 1);
	GFUNC gfun=((PAR_STRUCT*)par)->gfun;
	FLOATPTR G=mem+2 * N; 
	FLOATPTR THETA=((PAR_STRUCT*)par)->THETAfrac;
	int32_t  LEN=((PAR_STRUCT*)par)->Nfrac;
	gfun(THETA,G,LEN,gParam,NULL);
	FLOAT f;
	FLOAT df1;
	f=0;
	df1=0;
	for (rINT i=0; i < LEN; i++)
	{
		rFLOAT g=G[i];
		rFLOAT cost=cos(THETA[i]);
		{
			rFLOAT P0;
			rFLOAT frac=((PAR_STRUCT*)par)->FRAC[i];
			P0=exp(-LAI*g/cost);
			f+=(frac - P0)*(frac - P0);
			df1+=-(P0-frac)* P0 *g/cost;
		}
	}
	f=f*0.5;
	{
		FLOAT LAI1;
		FLOAT X1;
		X1=X[0]+dx;
		logitinv(&X1,&LAI1,bndL,1);
		df1=df1*(LAI1 - LAI)/dx;
	}
	df[0]=df1;
	for (rINT i=1; i < N; i++)
	{
		FLOATPTR Xnew=mem;
		memcpy(Xnew,X,sizeof(FLOAT)*N);
		Xnew[i]=Xnew[i]+dx;
		logitinv(Xnew+1,gParam,bndG,N - 1);
		FLOATPTR Gnew=mem+2 * N+LEN;
		gfun(THETA,Gnew,LEN,gParam,NULL);
		rFLOAT DF=0;
		for (rINT j=0; j < LEN; j++)
		{
			rFLOAT g=G[j];
			rFLOAT g1=Gnew[j];
			rFLOAT cost=cos(THETA[j]);
			{
				rFLOAT frac=((PAR_STRUCT*)par)->FRAC[j];
				rFLOAT P0=exp(-LAI*g/cost);
				DF+=-(P0 - frac)*P0*LAI/cost*(g1 - g)/dx;
			}
		}
		df[i]=DF;
	}
	return (FLOAT)(f);
}
FLOAT alg_L1F(FLOATPTR X,FLOATPTR df,int32_t N,void *par)
{
	FLOATPTR bndL=((PAR_STRUCT*)par)->boundL;
	FLOATPTR bndG=((PAR_STRUCT*)par)->boundG;
	FLOAT  dx=((PAR_STRUCT*)par)->dx;
	FLOATPTR mem=((PAR_STRUCT*)par)->mem;
	FLOAT LAI;
	FLOATPTR gParam=mem+N;   
	logitinv(X,&LAI,bndL,1);
	logitinv(X+1,gParam,bndG,N - 1);
	GFUNC gfun=((PAR_STRUCT*)par)->gfun;
	FLOATPTR G=mem+2 * N; 
	FLOATPTR THETA=((PAR_STRUCT*)par)->THETAfrac;
	int32_t  LEN=((PAR_STRUCT*)par)->Nfrac;
	gfun(THETA,G,LEN,gParam,NULL);
	FLOAT f;
	FLOAT df1;
	f=0;
	df1=0;
	for (rINT i=0; i < LEN; i++)
	{
		rFLOAT g=G[i];
		rFLOAT cost=cos(THETA[i]);
		{
			rFLOAT P0;
			rFLOAT frac=((PAR_STRUCT*)par)->FRAC[i];
			P0=exp(-LAI*g/cost);
			rFLOAT w=((PAR_STRUCT*)par)->WEIGHT[i];
			f+=w * (frac - P0)*(frac - P0);
			df1+=w*-(P0 - frac)* P0 *g/cost;
		}
	}
	f=f*0.5;
	{
		FLOAT LAI1;
		FLOAT X1;
		X1=X[0]+dx;
		logitinv(&X1,&LAI1,bndL,1);
		df1=df1*(LAI1 - LAI)/dx;
	}
	df[0]=df1;
	for (rINT i=1; i < N; i++)
	{
		FLOATPTR Xnew=mem;
		memcpy(Xnew,X,sizeof(FLOAT)*N);
		Xnew[i]=Xnew[i]+dx;
		logitinv(Xnew+1,gParam,bndG,N - 1);
		FLOATPTR Gnew=mem+2 * N+LEN;
		gfun(THETA,Gnew,LEN,gParam,NULL);
		rFLOAT DF=0;
		for (rINT j=0; j < LEN; j++)
		{
			rFLOAT g=G[j];
			rFLOAT g1=Gnew[j];
			rFLOAT cost=cos(THETA[j]);
			{
				rFLOAT frac=((PAR_STRUCT*)par)->FRAC[j];
				rFLOAT w=((PAR_STRUCT*)par)->WEIGHT[j];
				rFLOAT P0=exp(-LAI*g/cost);
				DF+=w*-(P0 - frac)*P0*LAI/cost*(g1 - g)/dx;
			}
		}
		df[i]=DF;
	}
	return (FLOAT)(f);
}
FLOAT alg_L1F_weight(FLOATPTR X,int32_t N,void *par)
{
	FLOATPTR bndL=((PAR_STRUCT*)par)->boundL;
	FLOATPTR bndG=((PAR_STRUCT*)par)->boundG;
	FLOAT  dx=((PAR_STRUCT*)par)->dx;
	FLOATPTR mem=((PAR_STRUCT*)par)->mem;
	FLOAT LAI;
	FLOATPTR gParam=mem+N;   
	logitinv(X,&LAI,bndL,1);
	logitinv(X+1,gParam,bndG,N - 1);
	GFUNC gfun=((PAR_STRUCT*)par)->gfun;
	FLOATPTR G=mem+2 * N; 
	FLOATPTR THETA=((PAR_STRUCT*)par)->THETAfrac;
	int32_t  LEN=((PAR_STRUCT*)par)->Nfrac;
	gfun(THETA,G,LEN,gParam,NULL);
	FLOAT f;
 	f=0;
	for (rINT i=0; i < LEN; i++)
	{
		rFLOAT g=G[i];
		rFLOAT cost=cos(THETA[i]);
		{
			rFLOAT P0;
			rFLOAT frac=((PAR_STRUCT*)par)->FRAC[i];
			P0=exp(-LAI*g/cost);
			rFLOAT w=(frac - P0);
			((PAR_STRUCT*)par)->WEIGHT[i]=1/(fabs(w)+FLOAT_EPSILON * 10);
		}
	}
	return (FLOAT)(f);
}
FLOAT alg_L2LF(FLOATPTR X,FLOATPTR df,int32_t N,void *par)
{
	FLOATPTR bndL=((PAR_STRUCT*)par)->boundL;
	FLOATPTR bndG=((PAR_STRUCT*)par)->boundG;
	FLOAT  dx=((PAR_STRUCT*)par)->dx;
	FLOATPTR mem=((PAR_STRUCT*)par)->mem;
	FLOAT LAI;
	FLOATPTR gParam=mem+N;   
	logitinv(X,&LAI,bndL,1);
	logitinv(X+1,gParam,bndG,N - 1);
	GFUNC    gfun=((PAR_STRUCT*)par)->gfun;
	FLOATPTR G=mem+2 * N; 
	FLOATPTR THETA=((PAR_STRUCT*)par)->THETAfrac;
	int32_t  LEN=((PAR_STRUCT*)par)->Nfrac;
	gfun(THETA,G,LEN,gParam,NULL);
	FLOAT f;
	FLOAT df1;
	f=0;
	df1=0;
	for (rINT i=0; i < LEN; i++)
	{
		rFLOAT g=G[i];
		rFLOAT cost=cos(THETA[i]);
		{
			rFLOAT P0;
			rFLOAT frac=((PAR_STRUCT*)par)->FRAC[i];
			rFLOAT logP=-LAI*g/cost;
			f+=(frac - logP)*(frac - logP);
			df1+=-(logP - frac)* g/cost;
		}
	}
	f=f*0.5;
	{
		FLOAT LAI1;
		FLOAT X1;
		X1=X[0]+dx;
		logitinv(&X1,&LAI1,bndL,1);
		df1=df1*(LAI1 - LAI)/dx;
	}
	df[0]=df1;
	for (rINT i=1; i < N; i++)
	{
		FLOATPTR Xnew=mem;
		memcpy(Xnew,X,sizeof(FLOAT)*N);
		Xnew[i]=Xnew[i]+dx;
		logitinv(Xnew+1,gParam,bndG,N - 1);
		FLOATPTR Gnew=mem+2 * N+LEN;
		gfun(THETA,Gnew,LEN,gParam,NULL);
		rFLOAT DF=0;
		for (rINT j=0; j < LEN; j++)
		{
			rFLOAT g=G[j];
			rFLOAT g1=Gnew[j];
			rFLOAT cost=cos(THETA[j]);
			{
				rFLOAT frac=((PAR_STRUCT*)par)->FRAC[j];
				rFLOAT logP=-LAI*g/cost;
				DF+=-(logP - frac)*LAI/cost*(g1 - g)/dx;
			}
		}
		df[i]=DF;
	}
	return (FLOAT)(f);
}
FLOAT alg_L1LF(FLOATPTR X,FLOATPTR df,int32_t N,void *par)
{
	FLOATPTR bndL=((PAR_STRUCT*)par)->boundL;
	FLOATPTR bndG=((PAR_STRUCT*)par)->boundG;
	FLOAT  dx=((PAR_STRUCT*)par)->dx;
	FLOATPTR mem=((PAR_STRUCT*)par)->mem;
	FLOAT LAI;
	FLOATPTR gParam=mem+N;   
	logitinv(X,&LAI,bndL,1);
	logitinv(X+1,gParam,bndG,N - 1);
	GFUNC gfun=((PAR_STRUCT*)par)->gfun;
	FLOATPTR G=mem+2 * N; 
	FLOATPTR THETA=((PAR_STRUCT*)par)->THETAfrac;
	int32_t  LEN=((PAR_STRUCT*)par)->Nfrac;
	gfun(THETA,G,LEN,gParam,NULL);
	FLOAT f;
	FLOAT df1;
	f=0;
	df1=0;
	for (rINT i=0; i < LEN; i++)
	{
		rFLOAT g=G[i];
		rFLOAT cost=cos(THETA[i]);
		{
			rFLOAT P0;
			rFLOAT frac=((PAR_STRUCT*)par)->FRAC[i];
			rFLOAT logP=-LAI*g/cost;
			rFLOAT w=((PAR_STRUCT*)par)->WEIGHT[i];
			f+=w*(frac - logP)*(frac - logP);
			df1+=w*-(logP - frac)* g/cost;
		}
	}
	f=f*0.5;
	{
		FLOAT LAI1;
		FLOAT X1;
		X1=X[0]+dx;
		logitinv(&X1,&LAI1,bndL,1);
		df1=df1*(LAI1 - LAI)/dx;
	}
	df[0]=df1;
	for (rINT i=1; i < N; i++)
	{
		FLOATPTR Xnew=mem;
		memcpy(Xnew,X,sizeof(FLOAT)*N);
		Xnew[i]=Xnew[i]+dx;
		logitinv(Xnew+1,gParam,bndG,N - 1);
		FLOATPTR Gnew=mem+2 * N+LEN;
		gfun(THETA,Gnew,LEN,gParam,NULL);
		rFLOAT DF=0;
		for (rINT j=0; j < LEN; j++)
		{
			rFLOAT g=G[j];
			rFLOAT g1=Gnew[j];
			rFLOAT cost=cos(THETA[j]);
			{
				rFLOAT frac=((PAR_STRUCT*)par)->FRAC[j];
				rFLOAT logP=-LAI*g/cost;
				rFLOAT w=((PAR_STRUCT*)par)->WEIGHT[j];
				DF+=w*-(logP - frac)*LAI/cost*(g1 - g)/dx;
			}
		}
		df[i]=DF;
	}
	return (FLOAT)(f);
}
FLOAT alg_L1LF_weight(FLOATPTR X,int32_t N,void *par)
{
	FLOATPTR bndL=((PAR_STRUCT*)par)->boundL;
	FLOATPTR bndG=((PAR_STRUCT*)par)->boundG;
	FLOAT  dx=((PAR_STRUCT*)par)->dx;
	FLOATPTR mem=((PAR_STRUCT*)par)->mem;
	FLOAT LAI;
	FLOATPTR gParam=mem+N;   
	logitinv(X,&LAI,bndL,1);
	logitinv(X+1,gParam,bndG,N - 1);
	GFUNC gfun=((PAR_STRUCT*)par)->gfun;
	FLOATPTR G=mem+2 * N; 
	FLOATPTR THETA=((PAR_STRUCT*)par)->THETAfrac;
	int32_t  LEN=((PAR_STRUCT*)par)->Nfrac;
	gfun(THETA,G,LEN,gParam,NULL);
	FLOAT f;
	FLOAT df1;
	f=0;
	df1=0;
	for (rINT i=0; i < LEN; i++)
	{
		rFLOAT g=G[i];
		rFLOAT cost=cos(THETA[i]);
		{
			rFLOAT P0;
			rFLOAT frac=((PAR_STRUCT*)par)->FRAC[i];
			rFLOAT logP=-LAI*g/cost;
			rFLOAT w=frac - logP;
			((PAR_STRUCT*)par)->WEIGHT[i]=1/(fabs(w)+FLOAT_EPSILON * 10);
		}
	}
	f=f*0.5;
	return (FLOAT)(f);
}
FLOAT alg_L2WL(FLOATPTR X,FLOATPTR df,int32_t N,void *par)
{
	FLOATPTR bndL=((PAR_STRUCT*)par)->boundL;
	FLOATPTR bndG=((PAR_STRUCT*)par)->boundG;
	FLOAT  dx=((PAR_STRUCT*)par)->dx;
	FLOATPTR mem=((PAR_STRUCT*)par)->mem;
	FLOAT LAI;
	FLOATPTR gParam=mem+N;   
	logitinv(X,&LAI,bndL,1);
	logitinv(X+1,gParam,bndG,N - 1);
	GFUNC gfun=((PAR_STRUCT*)par)->gfun;
	FLOATPTR G=mem+2 * N; 
	FLOATPTR THETA=((PAR_STRUCT*)par)->THETAfrac;
	int32_t  LEN=((PAR_STRUCT*)par)->Nfrac;
	gfun(THETA,G,LEN,gParam,NULL);
	FLOAT f;
	FLOAT df1;
	f=0;
	df1=0;
	for (rINT i=0; i < LEN; i++)
	{
		rFLOAT g=G[i];
		rFLOAT cost=cos(THETA[i]);
		{
			rFLOAT P0;
			rFLOAT frac=((PAR_STRUCT*)par)->FRAC[i];
			rFLOAT clogP=-LAI*g;
			f+=(frac - clogP)*(frac - clogP);
			df1+=-(clogP - frac)* g;
		}
	}
	f=f*0.5;
	{
		FLOAT LAI1;
		FLOAT X1;
		X1=X[0]+dx;
		logitinv(&X1,&LAI1,bndL,1);
		df1=df1*(LAI1 - LAI)/dx;
	}
	df[0]=df1;
	for (rINT i=1; i < N; i++)
	{
		FLOATPTR Xnew=mem;
		memcpy(Xnew,X,sizeof(FLOAT)*N);
		Xnew[i]=Xnew[i]+dx;
		logitinv(Xnew+1,gParam,bndG,N - 1);
		FLOATPTR Gnew=mem+2 * N+LEN;
		gfun(THETA,Gnew,LEN,gParam,NULL);
		rFLOAT DF=0;
		for (rINT j=0; j < LEN; j++)
		{
			rFLOAT g=G[j];
			rFLOAT g1=Gnew[j];
			rFLOAT cost=cos(THETA[j]);
			{
				rFLOAT frac=((PAR_STRUCT*)par)->FRAC[j];
				rFLOAT clogP=-LAI*g;
				DF+=-(clogP - frac)*LAI*(g1 - g)/dx;
			}
		}
		df[i]=DF;
	}
	return (FLOAT)(f);
}
FLOAT alg_L1WL(FLOATPTR X,FLOATPTR df,int32_t N,void *par)
{
	FLOATPTR bndL=((PAR_STRUCT*)par)->boundL;
	FLOATPTR bndG=((PAR_STRUCT*)par)->boundG;
	FLOAT  dx=((PAR_STRUCT*)par)->dx;
	FLOATPTR mem=((PAR_STRUCT*)par)->mem;
	FLOAT LAI;
	FLOATPTR gParam=mem+N;   
	logitinv(X,&LAI,bndL,1);
	logitinv(X+1,gParam,bndG,N - 1);
	GFUNC gfun=((PAR_STRUCT*)par)->gfun;
	FLOATPTR G=mem+2 * N; 
	FLOATPTR THETA=((PAR_STRUCT*)par)->THETAfrac;
	int32_t  LEN=((PAR_STRUCT*)par)->Nfrac;
	gfun(THETA,G,LEN,gParam,NULL);
	FLOAT f;
	FLOAT df1;
	f=0;
	df1=0;
	for (rINT i=0; i < LEN; i++)
	{
		rFLOAT g=G[i];
		rFLOAT cost=cos(THETA[i]);
		{
			rFLOAT P0;
			rFLOAT frac=((PAR_STRUCT*)par)->FRAC[i];
			rFLOAT clogP=-LAI*g/cost;
			rFLOAT w=((PAR_STRUCT*)par)->WEIGHT[i];
			f+=w*(frac - clogP)*(frac - clogP);
			df1+=w*-(clogP - frac)* g;
		}
	}
	f=f*0.5;
	{
		FLOAT LAI1;
		FLOAT X1;
		X1=X[0]+dx;
		logitinv(&X1,&LAI1,bndL,1);
		df1=df1*(LAI1 - LAI)/dx;
	}
	df[0]=df1;
	for (rINT i=1; i < N; i++)
	{
		FLOATPTR Xnew=mem;
		memcpy(Xnew,X,sizeof(FLOAT)*N);
		Xnew[i]=Xnew[i]+dx;
		logitinv(Xnew+1,gParam,bndG,N - 1);
		FLOATPTR Gnew=mem+2 * N+LEN;
		gfun(THETA,Gnew,LEN,gParam,NULL);
		rFLOAT DF=0;
		for (rINT j=0; j < LEN; j++)
		{
			rFLOAT g=G[j];
			rFLOAT g1=Gnew[j];
			rFLOAT cost=cos(THETA[j]);
			{
				rFLOAT frac=((PAR_STRUCT*)par)->FRAC[j];
				rFLOAT clogP=-LAI*g/cost;
				rFLOAT w=((PAR_STRUCT*)par)->WEIGHT[j];
				DF+=w*-(clogP - frac)*LAI *(g1 - g)/dx;
			}
		}
		df[i]=DF;
	}
	return (FLOAT)(f);
}
FLOAT alg_L1WL_weight(FLOATPTR X,int32_t N,void *par)
{
	FLOATPTR bndL=((PAR_STRUCT*)par)->boundL;
	FLOATPTR bndG=((PAR_STRUCT*)par)->boundG;
	FLOAT  dx=((PAR_STRUCT*)par)->dx;
	FLOATPTR mem=((PAR_STRUCT*)par)->mem;
	FLOAT LAI;
	FLOATPTR gParam=mem+N;   
	logitinv(X,&LAI,bndL,1);
	logitinv(X+1,gParam,bndG,N - 1);
	GFUNC gfun=((PAR_STRUCT*)par)->gfun;
	FLOATPTR G=mem+2 * N; 
	FLOATPTR THETA=((PAR_STRUCT*)par)->THETAfrac;
	int32_t  LEN=((PAR_STRUCT*)par)->Nfrac;
	gfun(THETA,G,LEN,gParam,NULL);
	FLOAT f;
	FLOAT df1;
	f=0;
	df1=0;
	for (rINT i=0; i < LEN; i++)
	{
		rFLOAT g=G[i];
		rFLOAT cost=cos(THETA[i]);
		{
			rFLOAT P0;
			rFLOAT frac=((PAR_STRUCT*)par)->FRAC[i];
			rFLOAT clogP=-LAI*g/cost;
			rFLOAT w=frac - clogP;
			((PAR_STRUCT*)par)->WEIGHT[i]=1/(fabs(w)+FLOAT_EPSILON * 10);
		}
	}
	f=f*0.5;
	return (FLOAT)(f);
}
