#' @title Plot arrow on plotly graph
#'
#' @description
#' \code{ArrowAnnot} create a list of list of parameter which will be useful
#' in the layout function for a plotly graph.
#' The aim is to print an arrow on the plotly graph, because we loose the arrow of geom_segment
#' when we use plotly. This create a parameter list for each vector of your data,
#' whether your have one vector or any number.
#'
#' @param data is a data.frame of four columns, containing the coordinates
#'  of the depart point (xd, yd) and the arrival point (xe, ye)
#'  of each vector you want to print (at least one)
#' @param d A numeric parameter corresponding to the Startstandoff parameter;
#'  allow to display the arrow with a specific length
#'
#' @return a list of list containing the parameter for the layout function
#'
#' @seealso \code{\link{ggplotly}} and \code{\link{layout}}
#'
#' @export

ArrowAnnot <- function(data, d = 10){
  n <- nrow(data)  #the number of vector
  whole_list <- list()
  for (i in 1:n){
    #e is for end and d for depart
    Start <- sqrt(((data$xe[i] - data$xd[i])^2) + ((data$ye[i] - data$yd[i])^2)) #length of the vector
    whole_list[[length(whole_list)+1]] <- list(x=data$xe[i],
                                               y=data$ye[i],
                                               showarrow=TRUE,
                                               xref = "x",
                                               yref = "y",
                                               arrowsize = 2,
                                               arrowhead =4,
                                               startstandoff = d*Start,
                                               ax=data$xd[i],
                                               ay=data$yd[i],
                                               axref="x",
                                               ayref="y")
  }
  return(whole_list)
}
