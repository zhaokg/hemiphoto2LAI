\name{sampleGapData}
\docType{data}
\alias{sampleGapData}
\title{A smaple gap dataset}
\description{
 A list object for a sample gap data
}
\usage{
    data(sampleGapData)
}
\value{
   
  \item{THETA}{a vector of zentihs. }
  \item{GAP}{a vector of binary values indicating gap or non-gap.}
  
  
}
\source{
hemiphoto2LAI v0.1
}
 
\examples{

\donttest{
data(sampleGapData) #will load two variables, THETA and GAP, into the R environment
plot(sampleGapData$THETA, sampleGapData$GAP)
result=estimateLAI(sampleGapData$THETA, sampleGapData$GAP)
}
 
}
\keyword{misc}
