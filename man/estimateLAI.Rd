\name{estimateLAI}
\alias{estimateLAI }
\title{Estimation of leaf area index and leaf angle distribution from hemiphoto}
\description{
   A compliation of 152 LAI estimation methods to invert LAI and LAD from hemiphotoss based on the gap probability theory.
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
  The output is an object of class "hemiphoto". It is a data frame with 19 rows; each row corresponds to one LAD type. The column variables include. 
  
  \item{LAD_TYPE}{three-letter acronyms for the 19 common leaf angle distributon (LAD) models, also known as Ross-Nilson functions. Type \code{print(LAD_list)} to see more details about the 19 LAD models. }
  \item{LAI_HA57}{an LAI estimate based on the simple hingle-angle algorithm from the gap fraction observed around the zenith of 57 degree. This algorithim is indepedent of LAD types and the full column is filled with the single LAI estimate throughout.}
  \item{LAI_Miller}{an LAI estimate using the Miller algorithm from gap fractions estimated from the input data. Similar to the hinge angle aglortihim, the Miller algorithim is indepedent of LAD types and the full column is also filled with the single LAI estimate throughout}
  \item{LAI_BNR}{LAI estimates using the binary nonliner regression method proposed in Zhao et al. (2019). These estiamtes are LAD-depedenent. Each row corresponds to one of the 19 LAD types.}
  \item{AIC_BNR}{AIC values calculated by the binary nonliner regression method for the 19 LAD model types. Smaller AIC values indicate beter LAI estimates and better LAD models. The LAD with the smallest AIC value is supposed to give the best LAI estimate.}
  \item{LAD1_BNR}{the LAD model parameters estimated by the binary nonlinear regression method. Of the 19 LAD models, some are parameter-free, as indicated by NAs in the LAD1_BNR or LAD2_BNR columns; some have only one parameter (i.e,, LAI2_BNR will be NAs); others have two parameters.}
  \item{LAD2_BNR}{the LAD model parameters estimated by the binary nonlinear regression method. This is the second LAD parameter for those models with two parameters}
  
  \item{LAI_L2F}{the next four columns are the same as the above, except that the aglorithms used are L2F--a set of classical algorithms to estimate LAI or LAD by mininizing the squared differences (i.e., L2 norm) between observed and modelled gap fractions. More details about the L2F algorithims are given in Zhao et al. (2019)}
  \item{AIC_L2F}{AIC values calculated by the L2F algorithm for individual LAD model types. Smaller AIC values indicate beter LAI estimates and better LAD models. }
  \item{LAD1_L2F}{the LAD model parameters estimated by the L2F algorithms. Of the 19 LAD models, some has zero parameters but others have up to two parameteres. LAD1_L2F and LAD2_L2F report the estimated parameter values. NAs indicate that the associated LAD model either has no or just one parameter.}
  \item{LAD2_L2F}{the LAD model parameters estimated by L2F. This is the second LAD parameter for those models with two parameters}
  
 \item{LAI_L1F}{the next four columns are the same as the above, except that the aglorithms used are L1F--a set of classical algorithms to estimate LAI or LAD by mininizing the absolute differences (i.e., L1 norm) between observed and modelled gap fractions. More details about the L1F algorithims are given in Zhao et al. (2019)}
  \item{AIC_L1F}{AIC values calculated by the L1F algorithm for individual LAD model types. Smaller AIC values indicate beter LAI estimates and better LAD models. }
  \item{LAD1_L1F}{the LAD model parameters estimated by the L1F algorithms. Of the 19 LAD models, some has zero parameters but others have up to two parameteres. LAD1_L1F and LAD2_L1F report the estimated parameter values. NAs indicate that the associated LAD model either has no or just one parameter.}
  \item{LAD2_L1F}{the LAD model parameters estimated by L1F. This is the second LAD parameter for those models with two parameters}
  
   \item{LAI_L2LF}{the next four columns are the same as the above, except that the aglorithms used are L2LF--a set of classical algorithms to estimate LAI or LAD by mininizing the squared differences (i.e., L2 norm) between observed and modelled log-gap fractions. More details about the L2F algorithims are given in Zhao et al. (2019)}
  \item{AIC_L2LF}{AIC values calculated by the L2LF algorithm for individual LAD model types. Smaller AIC values indicate beter LAI estimates and better LAD models. }
  \item{LAD1_L2LF}{the LAD model parameters estimated by the L2LF algorithms. Of the 19 LAD models, some has zero parameters but others have up to two parameteres. LAD1_L2LF and LAD2_L2LF report the estimated parameter values. NAs indicate that the associated LAD model either has no or just one parameter.}
  \item{LAD2_L2LF}{the LAD model parameters estimated by L2LF. This is the second LAD parameter for those models with two parameters}
  
 \item{LAI_L1LF}{the next four columns are the same as the above, except that the aglorithms used are L1LF--a set of classical algorithms to estimate LAI or LAD by mininizing the absolute differences (i.e., L1 norm) between observed and modelled log-gap fractions. More details about the L1LF algorithims are given in Zhao et al. (2019)}
  \item{AIC_L1LF}{AIC values calculated by the L1LF algorithm for individual LAD model types. Smaller AIC values indicate beter LAI estimates and better LAD models. }
  \item{LAD1_L1LF}{the LAD model parameters estimated by the L1LF algorithms. Of the 19 LAD models, some has zero parameters but others have up to two parameteres. LAD1_L2WL and LAD2_L2WL report the estimated parameter values. NAs indicate that the associated LAD model either has no or just one parameter.}
  \item{LAD2_L1LF}{the LAD model parameters estimated by L1LF. This is the second LAD parameter for those models with two parameters}
  
  
     \item{LAI_L2WL}{the next four columns are the same as the above, except that the aglorithms used are L2WL--a set of classical algorithms to estimate LAI or LAD by mininizing the squared differences (i.e., L2 norm) between observed and modelled "weighted log-gap fractions". More details about the L2WL algorithims are given in Zhao et al. (2019)}
  \item{AIC_L2WL}{AIC values calculated by the L2WL algorithm for individual LAD model types. Smaller AIC values indicate beter LAI estimates and better LAD models. }
  \item{LAD1_L2WL}{the LAD model parameters estimated by the L2WL algorithms. Of the 19 LAD models, some has zero parameters but others have up to two parameteres. LAD1_L2LF and LAD2_L2LF report the estimated parameter values. NAs indicate that the associated LAD model either has no or just one parameter.}
  \item{LAD2_L2WL}{the LAD model parameters estimated by L2WL. This is the second LAD parameter for those models with two parameters}
  
 \item{LAI_L1WL}{the next four columns are the same as the above, except that the aglorithms used are L1WL--a set of classical algorithms to estimate LAI or LAD by mininizing the absolute differences (i.e., L1 norm) between observed and modelled "weighted log-gap fractions". More details about the L1WL algorithims are given in Zhao et al. (2019)}
  \item{AIC_L1WL}{AIC values calculated by the L1WL algorithm for individual LAD model types. Smaller AIC values indicate beter LAI estimates and better LAD models. }
  \item{LAD1_L1WL}{the LAD model parameters estimated by the L1WL algorithms. Of the 19 LAD models, some has zero parameters but others have up to two parameteres. LAD1_L1WL and LAD2_L1WL report the estimated parameter values. NAs indicate that the associated LAD model either has no or just one parameter.}
  \item{LAD2_L1WL}{the LAD model parameters estimated by L1WL. This is the second LAD parameter for those models with two parameters}
}

\examples{

library(hemiphoto2LAI)

#--------------------------------Example 1--------------------------------#
data(sampleGapData) #will load two variables, THETA and GAP, into the R environment
plot(THETA, GAP)
result=estimateLAI(THETA, GAP)
#*****************************End of Example 1****************************#
 

}

}
\references{
Zhao et al. (2019). How to better estimate leaf area index and leaf angle distribution from digital hemispherical photography? Switching to a binary nonlinear regression paradigm (under review) 
}
\keyword{misc}
