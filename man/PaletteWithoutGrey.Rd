% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PaletteWithoutGrey.R
\name{PaletteWithoutGrey}
\alias{PaletteWithoutGrey}
\title{Create a color palette}
\usage{
PaletteWithoutGrey(markers, l = 2)
}
\arguments{
\item{markers}{A list containing all the markers of your data but it can be any list of character, factors, numeric, etc.}

\item{l}{l is a numeric argument.
if it is equal to 1 : same number of colors than the number of unique features of the markers list.
if it is equal to 2 : same number of colors than the number of unique features of the markers list minus 1.
(for example, if we want to set manually a color for one feature)}
}
\value{
A list of color which can be used for a plot
}
\description{
\code{PaletteWithoutGrey} create a color palette which not contains grey colors
 of the same length as the number of different markers of your \code{\link{MSnSet}} object.
}
\examples{
library(pRolocdata)
data(tan2009r1)
mypalette(fData(tan2009r1)$markers)
}
