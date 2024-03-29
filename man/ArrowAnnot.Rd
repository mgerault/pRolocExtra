% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ArrowAnnot.R
\name{ArrowAnnot}
\alias{ArrowAnnot}
\title{Plot arrow on plotly graph}
\usage{
ArrowAnnot(data, d = 10)
}
\arguments{
\item{data}{is a data.frame of four columns, containing the coordinates
of the depart point (xd, yd) and the arrival point (xe, ye)
of each vector you want to print (at least one)}

\item{d}{A numeric parameter corresponding to the Startstandoff parameter;
allow to display the arrow with a specific length}
}
\value{
a list of list containing the parameter for the layout function
}
\description{
\code{ArrowAnnot} create a list of list of parameter which will be useful
in the layout function for a plotly graph.
The aim is to print an arrow on the plotly graph, because we loose the arrow of geom_segment
when we use plotly. This create a parameter list for each vector of your data,
whether your have one vector or any number.
}
\seealso{
\code{\link{ggplotly}} and \code{\link{layout}}
}
