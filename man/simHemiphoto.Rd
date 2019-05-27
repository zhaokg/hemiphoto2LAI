\name{simHemiphoto}
\docType{methods}
\alias{simHemiphoto}
\title{Simulation of a hemiphoto}
\description{
   Generate a random hemiphoto of gap vs non-gap for a chosen leaf angle distribution model.
  }
\usage{
    simHemiphoto(ladType, ladParameter, option=list(),...)
}
\arguments{
  \item{ladType}{a three-leter string to specify one of the 19 LAD models supported (e.g., "\bold{es1}","\bold{uni}",and "\bold{lan}").  Type \code{print(LAD_list)} to see more details about all the 19 possible LAD models.}
  \item{ladParameter}{a vector of two elements specifying the parameter values of the chosen LAD type \code{ladType}. NAs should be filled if the LAD model is parameter-free or has just one parameter. Type \code{print(LAD_list)} to get more information about the valid parameter ranges for the 19 LAD models.}
  \item{option}{(optional). If absent, \code{simHemiphoto} will use default model parameters. If \code{option} is present, it must be a list variable specifying three input paramaters for the simulation algorithm. Possible parameters are option$LAI--LAI, option$imageSize--the image size, and option$scaledParameter-- a binary variable being either 1 or 0 to indicate how to interpret the second argument \code{ladParameter}.}  
  \item{\dots}{additional parameters, not used currently but reserved for future extension}   

}
\value{
  The output is a LIST object,comprising the following elements:
  
  \item{imageTHETA}{a square matrix of dimension \code{option$imageSize x option$imageSize }. It is the simulated image with pixel values indiciating zenith angles.}
  \item{imageGAP}{a square matrix of dimension \code{option$imageSize x option$imageSize }. It is the simulated binary image with pixel values indiciating gap (i.e., 1's) or non-gaps (i.e., 0's).}
  \item{LAI}{the LAI value used in the simulation. It should be the same as the inputted LAI as specified by option$LAI.}
  \item{PAR}{a copy of the input LAD parameters.}
  \item{THETA}{a vector consisting of 200 zenith angles ranging from 0 to pi/2. It is used to characterize the LAD curve used in simulation. }
  \item{G}{a vector consisting of 200 LAD value evalauted at the zenith agnles of \code{THETA}.}
  
}

\examples{

library(hemiphoto2LAI)

#--------------------------------Example 1--------------------------------#
opt=list()             #Create an empty list to append individual parameters
opt$LAI =3.0           #LAI of the simulated canopy
opt$imageSize=400      #The size of the resulting hemiphoto
opt$scaledParameter=0  #If set to 0, the ladParamter input should be valid values  
                       #falling within the expected range of the chosen LAD. If set
                       # to 1, ladParamter is two numbers within [0, 1] to 
                       # be scaled to the valid ranges of the chosen LAD model.
result=simHemiphoto("es1", c(0.1, NA), opt)
#*****************************End of Example 1****************************#
 

#\donttest
#{
LAI=estimateLAI(result$imageTHETA,result$imageGAP)
#}

}
\references{
Zhao et al. (2019). How to better estimate leaf area index and leaf angle distribution from digital hemispherical photography? Switching to a binary nonlinear regression paradigm (under review) 
}
\keyword{misc}
