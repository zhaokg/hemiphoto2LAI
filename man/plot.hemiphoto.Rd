\name{plot.hemiphoto}
\alias{plot.hemiphoto}
\title{Changepoint Detection}
\description{
  Plot the result obtained from the beast function.
  }
\usage{
  \method{plot}{hemiphoto}(x, index, ...)
}
\arguments{
  \item{x}{x must be an object of class "beast". It is the returned result from the beast function.}
  \item{index}{If x contains results for multiple time series. \code{indx} specifies which of them is plotted. }
  \item{...}{further arguments passed to the \code{\link[graphics]{plot}} function.}
}
\value{
   This function creates various plots to demonstrate the results of a beast decomposition. . 
   }
 
 
\examples{
\donttest{
library(hemiphoto2LAI)
 
}
}
\keyword{misc}
