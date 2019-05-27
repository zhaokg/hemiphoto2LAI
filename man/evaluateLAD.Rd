\name{evaluateLAD}
\docType{methods}
\alias{evaluateLAD}
\title{Evaluate the function values of a leaf area distribution}
\description{Compute the curve or function values of a chosen LAD }
\usage{   evaluateLAD(ladType, ladParameter, option=list(),...) }
\arguments{
  \item{ladType}{a three-leter string to specify one of the 19 LAD models supported (e.g., "\bold{es1}","\bold{uni}",and "\bold{lan}").  Type \code{print(LAD_list)} to see more details about all the 19 possible LAD models.}
  \item{ladParameter}{a vector of two elements specifying the parameter values of the chosen LAD type \code{ladType}. NAs should be filled if the LAD model is parameter-free or has just one parameter. Type \code{print(LAD_list)} to get more information about the valid parameter ranges for the 19 LAD models.}
  \item{option}{(optional). If absent, \code{evaluateLAD} will use default model parameters. If \code{option} is present, it must be a list variable consiting of two elements: "option$scaledParameter" as a binary variable being either 1 or 0 to indicate how to interpret the second argument \code{ladParameter}; "option$N" to specify how many annuli will the zenith range be divided into in evaluation of the LAD.}  
  \item{\dots}{additional parameters, not used currently but reserved for future extension}   

}
\value{
  The output is a LIST object,comprising the following elements:
  
  \item{THETA}{a vector of length \code{option$N x 1}, which is the locations of zenith angles at which the LAD will be evalauted.}
  \item{G}{a vector of length \code{option$N x 1}, giving the evaluted LAD values at the angles of THETA.}
  \item{PAR}{a copy of the input LAD parameters.}
  
}

\examples{

library(hemiphoto2LAI)

#--------------------------------Example 1--------------------------------#
opt=list()             #Create an empty list to append individual parameters
opt$N=100              #The full zenith range (i.e, [0, pi/2]) divided into 
                       #100 intervals
opt$scaledParameter=0  #If set to 0, the ladParamter input should be valid values  
                       #falling within the expected range of the chosen LAD. If set
                       # to 1, ladParamter is two numbers within [0, 1] to 
                       # be scaled to the valid ranges of the chosen LAD model.
result=simHemiphoto("es1", c(0.1, NA), opt)
#*****************************End of Example 1****************************#
 

}
\references{
Zhao et al. (2019). How to better estimate leaf area index and leaf angle distribution from digital hemispherical photography? Switching to a binary nonlinear regression paradigm (under review) 
}
\keyword{misc}
