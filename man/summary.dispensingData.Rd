% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/summary.dispensingData.R
\name{summary.dispensingData}
\alias{summary.dispensingData}
\title{Summarizing a dispensing data set}
\usage{
\method{summary}{dispensingData}(object, ...)
}
\arguments{
\item{object}{an object of class 'dispensingData'.}

\item{...}{further arguments passed to or from other methods.}
}
\value{
Returns summary statistics in form of a matrix for a given dispensing data set.
}
\description{
The function \emph{\code{summary}} generates a concise summary for a dispensing data set. The function provides insight into the scale of an experimental set-up.
}
\details{
This function extends the generic function \code{\link{summary}} for objects of class 'dispensingData'.
It returns the number of unique drug treatments, such as drug combinations (if, a drug combination screen is being carried out) and single drug treatments, as well as the number of drugs and controls used. 
It also provides the number of doses, and in cases where the number of doses differs between single drug treatments and drug combinations, each number is shown respectively (sing/comb). 
A more important number in the planning of a drug screen is the amount of plates an experimental set-up is leading up to, which will also be estimated by this function.
}
\examples{
\donttest{\dontrun{
summary(dispensingData)
}}

}
\seealso{
\code{\link{generateDispensingData}}
}
\keyword{combination}
\keyword{dispensing}
\keyword{drug}
\keyword{screen}
\keyword{summary}
