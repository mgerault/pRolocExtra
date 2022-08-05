#' @title Box plot profile for dynamic organellar map
#'
#' @description
#' \code{BoxProfile} plot a box plot of specific proteins between two condition (control vs treated).
#' Perform also an anova analysis on the protein selected.
#' The box are only useful when there are replicates.
#' This allow to visualize quickly if the protein profile changed between the two conditions.
#' Use only when you have dataset with two different conditions.
#'
#' @param obj A \code{\link{MSnSet}} object
#' @param prot A character vector which contains the protein(s) that you want to see the profiles between the two conditions
#' @param Condition The name of the column of your data which contains the condition for each protein (control or treated)
#' @param mytit Logical to tell if you want to print your own first element of the title
#' @param tit Character corresponding to the first element of the title, use when mytit = TRUE
#' @param cmet Character correspoonding to the column name containing the localization of each proteins
#'
#' @return a list containing the plot and the summary of the anova analysis
#'
#' @export

BoxProfile <- function(obj, prot = "GRB2", Condition = "cond", mytit = TRUE, tit = "Box profile", cmet = "markers"){

  prot_data <- as.data.frame(MSnbase::exprs(obj)[which(MSnbase::fData(obj)$gene_name == prot),])  #take the data from the protein chosen
  Loca <- unique(MSnbase::fData(obj)[[cmet]][which(MSnbase::fData(obj)$gene_name == prot)])       #take the marker from the protein
  if (cmet != "markers"){
    Loca <- paste("Clustered in", Loca, "for condition", MSnbase::fData(obj)[[Condition]][which(MSnbase::fData(obj)$gene_name == prot)], collapse = ", ")
  }
  else{
    Loca <- paste("Clustered in", Loca, "for all condition")
  }

  if (mytit){
    tit_ <- tit
  }
  else{
    tit_ <- deparse(substitute(obj))
  }

  ncond <- length(unique(MSnbase::fData(obj)[[Condition]])) #take the number different value in condition --> always two conditions but there can be replicates
  Condi <- c()

  for (i in 1:(ncond/2)){ #assume that the first condition condition are control, and the last are treated --> can be change
    Condi <- append(Condi, "C")
  }
  for (i in 1:(ncond/2)){
    Condi <- append(Condi, "T")
  }

  prot_data$cond <- factor(Condi)  #now we got only two conditions (no mattter the number of replicates) : C and T

  prot_data <- prot_data %>%       #change the organisation of the data in order to plot the two profiles
    dplyr::select(cond, colnames(prot_data)[1:5]) %>%
    tidyr::gather(key = "fraction", value = "value", -cond)

  #the plot is a box plot --> show the interest of the replicates if there is
  profile <- ggplot2::ggplot(prot_data, aes(x = fraction, y = value, colour = cond, fill = cond))+
    ggplot2::geom_point(position = position_jitterdodge(dodge.width = 0), size = 2)+
    ggplot2::geom_boxplot(alpha = 0.4, position = position_dodge(width=0), fatten = NULL)+
    ggplot2::stat_summary(fun = mean, geom = "errorbar", aes(ymax = ..y.., ymin = ..y..), #put the error bar in the box, position on the same fraction axe
                 width = 0.65, size = 1.5, linetype = "solid", position = position_dodge(width=0))+
    ggplot2::stat_summary(fun = mean, geom = "line", aes(group = cond), position = position_dodge(width=0), size = 1.5)+
    ggplot2::ggtitle(paste(tit_, ",", prot, Loca)) + #the localization in the title : the first is for the control condition, the second is for the treated (can change between the replicates if the localisation is calculated by clustering (error))
    ggplot2::theme_classic()

  fit <- stats::aov(value ~ cond*fraction, prot_data) #perform an anova on the data --> if the p-value is significant, the profiles are different between the conditions --> the protein have probably moved

  return(list(graph = profile, anova = summary(fit))) #return the plot and the anova results
}
