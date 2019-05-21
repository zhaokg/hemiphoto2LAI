\name{beast}
\alias{beast}
\title{Bayesian estimation of abrupt changepoint, nonlinear trend, and periodicity}
\description{
  Apply a Bayesian model averaging algorithm to decompose time series data into three contrasting components: Abrupt changes, trends, and cyclic/seasonal variations.
  }
\usage{
    beast(data, option=list(),demoGUI=FALSE,...)
}
\arguments{
  \item{data}{a vector or matrix of input data. Missing values are allowed and can be indicated by NA, NaN, or a pre-specified value using the option parameter (e.g., option$omittedValue=-9999). If \code{data} is a vector of dimension Nx1 or 1XN, it is just a single time series of length N. If it is a matrix (e.g., NxM), it contains multiple time series. Its row dimension "N" is the time series length; the col dimension "M" is the number of time series.}
  \item{option}{(optional). If absent, BEAST will use default model parameters. Foremost, the period of the cyclic component of time series will be best guessed via auto-correlation before running BEAST. "Period" must be an integer specifying the number of samples/values per cycle (e.g, a monthly-sampled time series with an annual period has a period of 12). If \code{option} is present, it can be either an integer specifying the period or a list specifying various paramaters for the BEAST decomposition algorithm. Possible parameters are demonstrated below in Example 3 of the Examples Secction.}
  \item{demoGUI}{a boolean indicator. If set to \code{TRUE}, BEAST will be run in a GUI demostration mode, with a GUI window to show an animation of the MCMC sampling in the model space step by step. Note that "\code{demoGUI=TRUE}"  works only for Windows x64 systems not Windows 32 or Linux/Mac systems.}
  \item{\dots}{additional parameters, not used currently but reserved for future extension}
}
\value{
  The output is an object of class "beast". It is a list, consistinig of the following components. Suppose that the input data is of size \code{NxM}: N is the time series length and M is the number of time series (\code{M>=1}): 
  
  \item{time}{a vector of size \code{1xN}: the times at the N sampled locatons. By default, it is simply set to 1:N}
  \item{sN}{a vector of size \code{1xM}. \code{sN} gives the mean number of seasonal changepoints for each of the M time series. If \code{data} is a single time series (i.e., M=1), sN will be a scalar.}
  \item{tN}{a vector of size \code{1xM}. \code{tN} gives the mean number of trend changepoints for each of the M time series. If \code{data} is a single time series, tN will be a scalar.}
  \item{sNProb}{a matrix of size (\code{option$maxKnotNum_Season+1)xM}. For the i-th time series (i.e.,\code{1=<i<=M}), sNProb[,i] gives a probability distribution of having a certain number of changepoints over the range of [0,option$maxKnotNum_Season]; for example, sNProb[1,i] is the probability of having no seasonal changepoint in the ith time series; sNProb[3,i] is the probability of having two changepoints.}
  \item{tNProb}{a matrix of size (option$maxKnotNum_Trend+1)xM. For the i-th time series (i.e.,\code{1=<i<=M}), tNProb[,i] gives a probability distribution of having a certain number of trend changepoints over the range of [0,option$maxKnotNum_Trend]; for example, tNprob[1,i] is the probability of having no trend changepoint in the ith time series; sTprob[4,i] is the probability of having 3 (i.e, 4-1) changepoints.}
  \item{sProb}{a matrix of size \code{NxM}. For the i-th time series (i.e., \code{1=<i<=M}), sProb[,i] gives a probability distribution of having a seasonal changepoint at a certain time over the time range of [1,N]. Plotting sProb[,i] will depict a curve of probability-of-being-changepoint over the time for the i-th time series.}
  \item{tProb}{a matrix of size \code{NxM}. For the i-th time series (i.e., \code{1=<i<=M}), tProb[,i] gives a probability distribution of having a trend changepoint at a certain time over the time range of [1,N]. Plotting tProb[,i] will depict a curve of probability-of-being-trend-changepoint over the time for the i-th time series.}
  \item{s}{a matrix of size \code{NxM}. For the i-th time series (i.e., \code{1=<i<=M}), s[,i] gives the best fitted seasonal component.}
  \item{t}{a matrix of size \code{NxM}. For the i-th time series (i.e., \code{1=<i<=M}), t[,i] gives the best fitted trend component.}
  \item{b}{a matrix of size \code{NxM}. For the i-th time series (i.e., \code{1=<i<=M}), b[,i] gives the estimated slope in the fitted trend component over time.}
  \item{sSD}{a matrix of size \code{NxM}. For the i-th time series (i.e., \code{1=<i<=M}), sSD[,i] gives the standard deviations of the fitted SEASONAL component.}
  \item{tSD}{a matrix of size \code{NxM}. For the i-th time series (i.e., \code{1=<i<=M}), tSD[,i] gives the standard deviations of the fitted TREND component.}
  \item{bSD}{a matrix of size \code{NxM}. For the i-th time series (i.e., \code{1=<i<=M}), bSD[,i] gives the standard deviations of estimated SLOPE.}  
  \item{sCI}{a matrix of size \code{NxM}. For the i-th time series (i.e., \code{1=<i<=M}), sCI[,i] gives the 95\% credible intervals of the fitted SEASONAL component.}  
  \item{tCI}{a matrix of size \code{NxM}. For the i-th time series (i.e., \code{1=<i<=M}), tCI[,i] gives the 95\% credible intervals of the fitted TREND component.} 
  \item{bCI}{a matrix of size \code{NxM}. For the i-th time series (i.e., \code{1=<i<=M}), bCI[,i] gives the the 95\% credible intervals of the estimated SLOPE.}
  \item{horder}{a matrix of size \code{NxM}. For the i-th time series (i.e., \code{1=<i<=M}), horder[,i] gives the estimated harmonic order used to approximate the SEASONAL component over time.}
  \item{torder}{a matrix of size \code{NxM}. For the i-th time series (i.e., \code{1=<i<=M}), torder[,i] gives the estimated polynomial order used to approximate the TREND component over time.}
  \item{scp}{a matrix of size \code{(option$maxKnotNum_Season)xM}. For the i-th time series (i.e., \code{1=<i<=M}), scp[,i] gives the most possible locations of changepoints in the  the SEASONAL component.}
  \item{tscp}{a matrix of size \code{(option$maxKnotNum_Trend)xM}. For the i-th time series (i.e., \code{1=<i<=M}), tcp[,i] gives gives the most possible locations of changepoints in the TREND component.}
  

}

\examples{

library(Rbeast)

#A MODIS time series of NDVI for a forest plot in Ohio. It has
# 23 samples per year (i.e., period=23).

data(modis_ohio)
plot(modis_ohio)

\donttest{
#-------------------Example 1--------------------------------#
#No option supplied below so default parameters are used. The period 
#(i.e., 23) will be estimated via auto-correlation. This is not
# always reliable,so it is always suggested to directly supply 
# the period as in Example 2 and Example 3.

out=beast(modis_ohio)
plot(out)


#-------------------Example 2--------------------------------#
#option is set to 23, specicfying the period of modis_ohio as 23

out=beast(modis_ohio,23)
plot(out)

plot(out$s)     #the same as plot(out$s[,1]): plot the seasonal curve
plot(out$sProb) #plot the probability of observing seasonal changepoints
plot(out$t)     #the same as plot(out$t[,1]): plot the trend
plot(out$sProb) #plot the probability of observing trend changepoints
}


#-------------------Example 3--------------------------------#
#Specify the option parameters explicilty 

opt=list()              #create an empty list to append individual parameters

opt$period=23           #period of the cyclic/seasonal component of the time series
opt$minSeasonOrder=4    #min harmonic order allowed in fitting season component
opt$maxSeasonOrder=8    #max harmonic order allowed in fititing season component
opt$minTrendOrder=0     #min polynomial order allowed to fit trend (0 for constant)
opt$maxTrendOrder=1     #max polynomial order allowed to fit trend (1 for linear term)
opt$minSepDist_Season=20#min seperation distance btw neighboring season changepoints(must be >=0)
opt$minSepDist_Trend=20 #min seperation distance btw neighboring trend changepoints(must be >=0)
opt$maxKnotNum_Season=4 #max number of season changepoints allowed 
opt$maxKnotNum_Trend=10 #max number of trend changepoints allowed  

#A customized value to indicate bad/missing values in the time series, in additon to 
#those NA or NaN values.
opt$omittedValue=-999   


#The following parameters used to configure the reverisible-jump MCMC (RJMCC) sampler
opt$chainNumber=2     #number of parallel MCMC chains 
opt$sample=1000       #number of samples to be collected per chain
opt$thinningFactor=3  #a factor to thin the chains (e.g., samples taken every 3 iterations)
opt$burnin=500        #number of burn-in samples discarded at the start of each chain

#For the move proposal, the max window allowed in jumping from a current changepoint
opt$maxMoveStepSize=30
#The probability of selecting a re-sampling proposal (e.g., resample seasonal harmonic order)
opt$resamplingSeasonOrderProb=0.2  
#The probability of selecting a re-sampling proposal (e.g., resample trend polynomial order)
opt$resamplingTrendOrderProb=0.2  

#A seed for the random generator: If seed=0,random numbers differ for different BEAST runs. 
#Setting seed to a chosen non-zero integer allows reproducing the same result for different
#BEAST runs.
opt$seed=65654   

#If set to 1, compute 95% credible intervals: The results will be saved as 
# sCI, tCI, and bCI in the output variable.
opt$computeCredible=0   

#If set to 1, employ a fast algorithm to compute credible intervals
opt$fastCIComputation=0 

#If set to 1, compute the probability of having a postive slope in the estimated trend.
# The result will be saved as bsign in the output variable.
opt$computeSlopeSign=1  

#If set to 1, compute the mean harmonic order of the fitted seasonal component.
# The result will be saved as horder in the output variable.
opt$computeHarmonicOrder=1  

#If set to 1, compute the mean polynomial order of the fitted trend component.
# The result will be saved as torder in the output variable.
opt$computeTrendOrder=1  


#opt$outputToDisk=0            #(if set to 1, results will be written to files in a folder)
#opt$outputFolder ='c:/out'    #Specify the output folder when outputToDisk=1     
#opt$lengthPerTimeSeries_infile=300#the time series length if input data come from a binary file

opt$printToScreen=1            #if set to 1, display some progress status while running
opt$printCharLen=150           #the length of chars in each status line when printToScreen=1

#Use "opt" defined above in the beast function. Note that to run beast(), not all the individual 
#parameters in option need to be explicitly specified. If an parameter is not given in option, its 
#default value will be used.

\donttest{
out=beast(modis_ohio, opt)
plot(out)
}

\donttest{
#-------------------Example 4--------------------------------#
#Run an interactive GUI to visualize how BEAST is samplinig from
# the possible model spaces in terms of the numbers and timings of 
# seasonal and trend changepoints. The GUI inferface allows changing
# the option parameters interactively.
 
beast(modis_ohio, 23, demoGUI=TRUE)
 

#-------------------Example 5--------------------------------#
#'simdata' is a 300x3 matrix, consisting of three time series

data(simdata)

#Plot individual time series. As a toy example, all the three time series
#are the same.
plot(simdata[,1]) 
plot(simdata[,2]) 

#Below, the option is defined in the command line as a temporary list.
out=beast( simdata, list(period=24, chainNumber=3, sample=1000, burnin=200) )

#"out" contains results for the three time series. Plot the result for the second one
plot(out,2)

}

}
\references{
Zhao K, Wulder M, Hu, T, Bright R, Wu Q, Qin H, Li Y, Toman E, Mallick B, Zhang X, and Brown M (2019).  Detecting change-point, trend, and seasonality in satellite time series data to track abrupt changes and nonlinear dynamics: A Bayesian ensemble algorithm. Remote Sensing of Environment (accepted manuscript) 
}
\keyword{misc}
