#' @title Fit the pRolocdata format
#'
#' @description
#' \code{nainunknown} transform all NA values in the marker column of the fData of the object in "unknown".
#'
#' @param object A \code{\link{MSnSet}} object
#' @param fcol The name of the column where are the markers (organelles names)
#'
#' @return The object without NA values in the marker column of its  fData
#'
#' @export
#'
#' @examples
#'
#' library(pRolocExtra)
#' nainunknown(tan2009r1)

nainunknown <- function(object, fcol = "markers"){
  naind <- which(is.na(MSnbase::fData(object)[,fcol]))
  MSnbase::fData(object)[,fcol][naind] <- "unknown"
  return(object)
}
