#' @title Filter your data
#'
#' @description
#' \code{nafilterdata} take off all proteins from your data set with NA or infinite value.
#'
#' @param object A \code{\link{MSnSet}} object
#'
#' @return The object without the proteins with NA or infinite value
#'
#' @export
#'
#' @examples
#' library(pRolocdata)
#' data(tan2009r1)
#' nafilterdata(tan2009r1)
nafilterdata <- function(object) {
  mat = MSnbase::exprs(object)                                     #get the experiment data

  prot.has.na = apply(mat, 1, function(x) sum(is.na(x)))           #find the proteins with NA values
  prot.has.inf = apply(mat, 1, function(x) sum(is.infinite(x)))    #find the proteins with infinite values

  return(object[prot.has.na == 0 & prot.has.inf == 0])
}
