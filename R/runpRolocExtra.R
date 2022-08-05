#' @title Run the app to visualize all data from pRolocdata
#'
#' @description
#'  \code{runpRolocExtra} works exactly the same as runExample from \code{\link{shiny}} package.
#'
#' @export
runpRolocExtra <- function() {
 appDir <- system.file("shiny-examples", "myapp", package = "pRolocExtra")
  if (appDir == "") {
    stop("Could not find example directory. Try re-installing `pRolocExtra`.", call. = FALSE)
  }
  shiny::runApp(appDir, display.mode = "normal")
}
