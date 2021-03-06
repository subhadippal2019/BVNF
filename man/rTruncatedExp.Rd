% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rMPG_nu_Half.R
\name{rTruncatedExp}
\alias{rTruncatedExp}
\title{Generates a sample from the truncated Exponential random variable
 rTruncatedExp(n=1, rate=1,lowerLim=0, upperLim=Inf, eps=1e-100 )}
\usage{
rTruncatedExp(n = 1, rate = 1, lowerLim = 0, upperLim = Inf, eps = 1e-100)
}
\arguments{
\item{n}{number of samples to be generated.}

\item{rate}{rate of the Exponential Distribution.}

\item{lowerLim}{Lower limit of the support. A number greater than 0.}

\item{upperLim}{UpperLimit of the Support. A number greater then \code{lowerLim}}

\item{eps=1e-100}{Incase a discrete sampler is to be deployed, \code{eps} pertains to the bin size.}
}
\value{
Samples from the truncated Exponential distribution
}
\description{
Generates a sample from the truncated Exponential random variable
 rTruncatedExp(n=1, rate=1,lowerLim=0, upperLim=Inf, eps=1e-100 )
}
\examples{
rTruncatedExp(n=1, rate=2)
rTruncatedExp(n=1, rate=200, lowerLim=20, upperLim=21)
}
