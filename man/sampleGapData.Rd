\name{sampleGapData}
\docType{data}
\alias{sampleGapData}
\title{A smaple gap dataset}
\description{
  \code{sampleGapData} includes two variables--THETA and GAP--the zeniths and gap indicators.
}
\usage{
    data(sampleGapData)
}
\source{
hemiphoto2LAI v0.1
}
 
\examples{

data(sampleGapData) #will load two variables, THETA and GAP, into the R environment
plot(THETA, GAP)
result=estimateLAI(THETA, GAP)
}
}
\keyword{misc}
