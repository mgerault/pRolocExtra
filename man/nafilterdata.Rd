% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nafilterdata.R
\name{nafilterdata}
\alias{nafilterdata}
\title{Filter your data}
\usage{
nafilterdata(object)
}
\arguments{
\item{object}{A \code{\link{MSnSet}} object}
}
\value{
The object without the proteins with NA or infinite value
}
\description{
\code{nafilterdata} take off all proteins from your data set with NA or infinite value.
}
\examples{

library(pRolocExtra)
nafilterdata(tan2009r1)
}
