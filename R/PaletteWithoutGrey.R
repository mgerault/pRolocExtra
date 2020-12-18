#' @title Create a color palette
#'
#' @description
#' \code{PaletteWithoutGrey} create a color palette which not contains grey colors
#'  of the same length as the number of different markers of your \code{\link{MSnSet}} object.
#'
#' @param markers  A list containing all the markers of your data but it can be any list of character, factors, numeric, etc.
#' @param l l is a numeric argument.
#'  if it is equal to 1 : same number of colors than the number of unique features of the markers list.
#'  if it is equal to 2 : same number of colors than the number of unique features of the markers list minus 1.
#'  (for example, if we want to set manually a color for one feature)
#'
#' @return A list of color which can be used for a plot
#'
#' @export
#'
#' @examples
#' library(pRolocdata)
#' data(tan2009r1)
#' mypalette(fData(tan2009r1)$markers)
PaletteWithoutGrey <- function(markers, l = 2){

  n = length(unique(markers))
  x <- grDevices::colors(distinct = TRUE)                           #all the color from R
  mycol <- x[which(is.na(stringr::str_extract(x, "gr(e|a)y")))]   #keep only colors that are not grey

  listcolor <- c()
  for (i in 0:(n-l))
    listcolor <- append(listcolor, mycol[i*20 + 9])      #save a color from the list (the number 20 and 9 were chosen in order to have distincts colors, this is empirical, can be changed)

  return(listcolor)
}
