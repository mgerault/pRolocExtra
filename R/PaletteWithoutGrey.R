#' @title Create a color palette
#'
#' @description
#' \code{PaletteWithoutGrey} create a color palette which not contains grey colors
#'  of the same length as the number of different markers of your \code{\link{MSnSet}} object.
#'
#' @param markers  A list containing all the markers of your data but it can be any list of character, factors, numeric, etc.
#'
#' @return A list of color which can be used for a plot
#'
#' @export

PaletteWithoutGrey <- function(markers){
  listcolor <- list()

  mrk <- as.character(unique(markers))
  mrk <- mrk[order(mrk)]
  if("unknown" %in% mrk){
    mrk <- mrk[!(mrk %in% "unknown")]
    listcolor[["unknown"]] <- "darkgrey"
  }

  x <- grDevices::colors(distinct = TRUE)                         #all the color from R
  mycol <- x[-stringr::str_which(x, "gr(e|a)y")]   #keep only colors that are not grey

  #save a color from the list (the number 20 and 9 were chosen in order to have distincts colors, this is empirical, can be changed)
  for (i in 1:length(mrk)){
    listcolor[[mrk[i]]] <- mycol[i*20 + 9]
  }

  return(unlist(listcolor))
}
