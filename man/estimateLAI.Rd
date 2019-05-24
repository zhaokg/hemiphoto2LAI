\name{estimateLAI}
\alias{estimateLAI }
\title{Estimation of leaf area index and leaf angle distribution from hemiphoto}
\description{
  Apply the gap-probability theory to invert LAI and LAD from hemiphotos using more than 150 methods.
  }
\usage{
    estimateLAI(THETA, GAP, option=list(),...)
}
\arguments{
  \item{THETA}{a vector or matrix of input data representing pixel-level zenith angles in radians. Missing values are allowed and can be indicated by NA or NaN.}
  \item{GAP}{a vector or matrix of input data representing pixel-level gap status: 0 for non-gap and 1 for gap. Missing values are allowed and can be indicated by NA or NaN.}
  \item{option}{(optional). If absent, estimateLAI will use default model parameters. If \code{option} is present, it can be a list variable specifying various paramaters for the LAI estimation algorithm. Possible parameters are demonstrated below in Example 3 of the Examples Secction.}  
  \item{\dots}{additional parameters, not used currently but reserved for future extension}
}
\value{
  The output is an object of class "hemiphoto". It is a data frame with 19 rows; each row corresponds to one LAD type. The columns include. 
  
  \item{LAD_TYPE}{three-letter acronyms for the 19 common leaf angle distributon (LAD) models, also known as Ross-Nilson functions.}
  \item{LAI_HA57}{an LAI estimate based on the simple hingle-angle algorithm from the gap fraction observed around the zenith of 57 degree. This algorithim is indepedent of LAD types and the full column is filled with the single LAI estimate throughout.}
  \item{LAI_Miller}{an LAI estimate using the Miller algorithm from gap fractions estimated from the input data. Similar to the hinge angle aglortihim, the Miller algorithim is indepedent of LAD types and the full column is also filled with the single LAI estimate throughout}
  \item{LAI_BNR}{LAI estimates using the binary nonliner regression method proposed in Zhao et al. (2019). These estiamtes are LAD-depedenent. Each row corresponds to one of the 19 LAD types.}
  \item{AIC_BNR}{AIC values calculated by the binary nonliner regression method for the 19 LAD model types. Smaller AIC values indicate beter LAI estimates and better LAD models. The LAD with the smallest AIC value is supposed to give the best LAI estimate.}
  \item{LAD1_BNR}{the LAD model parameters estimated by the binary nonlinear regression method. Of the 19 LAD models, some are parameter-free, as indicated by NAs in the LAD1_BNR or LAD2_BNR columns; some have only one parameter (i.e,, LAI2_BNR will be NAs); others have two parameters.}
  \item{LAD2_BNR}{the LAD model parameters estimated by the binary nonlinear regression method. This is the second LAD parameter for those models with two parameters}
  
  \item{LAI_BNR}{LAI estimates using the binary nonliner regression method proposed in Zhao et al. (2019). These estiamtes are LAD-depedenent. Each row corresponds to one of the 19 LAD types.}
  \item{AIC_BNR}{AIC values calculated by the binary nonliner regression method for the 19 LAD model types. Smaller AIC values indicate beter LAI estimates and better LAD models. The LAD with the smallest AIC value is supposed to give the best LAI estimate.}
  \item{LAD1_BNR}{the LAD model parameters estimated by the binary nonlinear regression method. Of the 19 LAD models, some are parameter-free, as indicated by NAs in the LAD1_BNR or LAD2_BNR columns; some have only one parameter (i.e,, LAI2_BNR will be NAs); others have two parameters.}
  \item{LAD2_BNR}{the LAD model parameters estimated by the binary nonlinear regression method. This is the second LAD parameter for those models with two parameters}
  

}

\examples{

library(hemiphoto2LAI)

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
Zhao et al. (2019). How to better estimate leaf area index and leaf angle distribution from digital hemispherical photography? Switching to a binary nonlinear regression paradigm (under review) 
}
\keyword{misc}
