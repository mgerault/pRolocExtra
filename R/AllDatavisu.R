#' @title Complete data visualization
#'
#' @description
#' \code{AllDatavisu} allows you to visualize your data with many options.
#' You can choose the reduction method (PCA, t-SNE, ...), the clustering method, export the plot as a png or pdf file.
#' The function adapts according to the data, so if you have replicates or two different experiment conditions you can visualize it.
#' This is very useful when you work with dynamic organellar map.
#' With this type of data you can visualize the movement of protein between the two conditions with a vector.
#'
#' @param object2 A \code{\link{MSnSet}} object
#' @param redmet The reduction method from the pRoloc package : "PCA", "t-SNE", "MDS", "nipals", "lda", "kpca".
#' @param cmet the name of the column of your data which contains the markers, if it contains "unknown" assignment, please precise unknow = TRUE
#' @param ax A numeric vector of length two, you can choose on which axes you want to see the plot.
#' (depend of the number of fraction of the data)
#' @param Title The title of the plot, use only when vect = TRUE
#' @param Condition A character argument, it is the name of the column of your data which contains the condition of the experiment (control or treated).
#'  You can have replicates.
#'  If you have no dynamic experiment, don't bother with this parameter.
#' @param highpr A logical argument to highlight or not specific protein
#' @param proteins A character vector containing the proteins you want to highlight
#' @param vect A logical argument to use when you have dynamic experiment.
#'  If TRUE, allows you to see the proteins movement between the two conditions with vectors.
#'  You also have to choose the proteins you want to see with the proteins argument.
#' @param same_size A logical argument to tell if you want all the proteins having the same size on the plot.
#'  if FALSE, the size is the score of the clustering algorithm (also if unknow = FALSE).
#' @param sz An integer for the size you want (if same_size = TRUE)
#' @param unknow A logical argument to tell if the marker column contain "unknown" features
#' @param expor.png A logical argument to export the figure in a png file
#' @param expor.pdf A logical argument to export the figure in a pdf file
#' @param yourseed An integer for the t-SNE algorithm in order to having same plot if there is several
#'
#' @return A figure showing your data, depending of your chosen parameters
#'
#' @seealso \code{\link{plot2D}} and \code{\link{svmOptimisation}} from pRoloc package for more details
#'
#' @export
#'
#' @examples
#' library(pRolocdata)
#' data(tan2009r1)
#' AllDatavisu(datavisupca(tan2009r1, method ="knn", sh.gr = FALSE), redmet = "t-SNE", cmet = "knn")
AllDatavisu <- function(object2, redmet = "PCA", cmet = "svm", ax = c(1,2),
                        Title = "Cellular map", Condition = "cond",
                        highpr = FALSE, proteins = c("PKN2", "GRB2", "SHC1", "EGFR"),
                        vect = FALSE, same_size = FALSE, sz = 1, unknow = FALSE,
                        expor.png = FALSE, expor.pdf = FALSE, yourseed = 500){

  set.seed(yourseed) #set the seed for the t-SNE

  print("1. Create the data plot")
  p2 <- pRoloc::plot2D(object2, fcol = cmet, method = redmet, plot = FALSE,  #get the data from the plot2D function of pRoloc package (see documentation of pRoloc)
               dims = ax)
  p2 <- as.data.frame(p2)
  axname <- colnames(p2)  #save the names of the axes


  print("2. Data process")

  if (redmet != "t-SNE") #t-SNE only work with unique features --> filter first
    p2$markers <- MSnbase::fData(object2)[[cmet]]

  else {
    if (length(unique(object2@assayData$exprs)[,1]) == length(MSnbase::fData(object2)[[cmet]]))
      p2$markers <- MSnbase::fData(object2)[[cmet]]

    else
      p2$markers <- MSnbase::fData(object2)[[cmet]][-which(duplicated(object2@assayData$exprs) == TRUE)]
  }

  p2$markers <- as.factor(p2$markers)

  if (!unknow){ #set the size of the points on the plot according their clustering score or not
    rep.loss <- round(mean(MSnbase::fData(object2)[[paste0(cmet, ".scores")]]),3)
    capt <- paste("The score of", cmet, "clustering method is", rep.loss)
    if (!same_size){
      if (redmet != "t-SNE")
        p2$score <- as.factor(exp(MSnbase::fData(object2)[[paste0(cmet, ".scores")]]) - 1)
      else {
        if (length(unique(object2@assayData$exprs)[,1]) == length(MSnbase::fData(object2)[[cmet]]))
          p2$score <- as.factor(exp(MSnbase::fData(object2)[[paste0(cmet, ".scores")]]) - 1)

        else
          p2$score <-  as.factor(exp(MSnbase::fData(object2)[[paste0(cmet, ".scores")]][-which(duplicated(object2@assayData$exprs) == TRUE)]) - 1)
      }
    }

    if (same_size){
      p2$score <- as.factor(matlab::linspace(sz,sz, nrow(p2)))
    }
  }

  else { #if unknow = TRUE, no clustering algorithm was performed, so to change the size of the unknown features on the plot we have to define it
    p2$score <- p2$markers
    capt <- NULL
  }

  colnames(p2) <- c("Dim1", "Dim2", "markers", "score")    #reorganize the data, adding new columns (save the proteins names)
  p2$prot <- rownames(p2)
  p2$unprot <- MSnbase::fData(object2)$gene_name


  print("3. Building the base plot")

  grap2 <- ggplot2::ggplot(data = p2, aes(x=Dim1, y=Dim2, fill=markers,
                                 size = score)) +
    ggplot2::geom_point(alpha = 0.7, shape = 21) + ggplot2::guides(alpha = FALSE, size = FALSE)



  if (is.null(MSnbase::fData(object2)[[Condition]])){  #if no condition in the data (no dynamic)
    MSnbase::fData(object2)[[Condition]] <- matlab::linspace(1,1,nrow(fData(object2))) #create one in order to have always the same data organization for after
    Main <- c("Map") #the plot title
    k <- length(unique(MSnbase::fData(object2)[[Condition]]))
    all_cond <- unique(MSnbase::fData(object2)[[Condition]])
    if (redmet != "t-SNE")
      p2$cond <- MSnbase::fData(object2)[[Condition]]
    else{
      if (length(unique(object2@assayData$exprs)[,1]) == length(MSnbase::fData(object2)[[cmet]]))
        p2$cond <- MSnbase::fData(object2)[[Condition]]
      else
        p2$cond <- MSnbase::fData(object2)[[Condition]][-which(duplicated(object2@assayData$exprs) == TRUE)]
    }

  }

  else{
    k <- length(unique(MSnbase::fData(object2)[[Condition]]))  #the number of condition (always two but we can have replicates)
    all_cond <- unique(MSnbase::fData(object2)[[Condition]])   #save all the conditions names of the data
    if (redmet != "t-SNE")
      p2$cond <- MSnbase::fData(object2)[[Condition]]
    else{
      if (length(unique(object2@assayData$exprs)[,1]) == length(MSnbase::fData(object2)[[cmet]]))
        p2$cond <- MSnbase::fData(object2)[[Condition]]
      else
        p2$cond <- MSnbase::fData(object2)[[Condition]][-which(duplicated(object2@assayData$exprs) == TRUE)]
    }

    Main <- c("Control", "Treated")                   #creation of the plot title according to the number of condition
    if (k > 2){
      for (i in 1:(k/2)){
        Main[i] <- paste0("Control_", i)
        Main[k/2 +i] <- paste0("Treated_", i)
      }
    }
  }


  print("4. Building the plot(s)")

  if (!highpr & !vect  & !unknow){ #creating the plot according to the number of condition, then highlight the condition
    Gr <- list()
    for (i in 1:k){
      print(i)
      graph1a <- ggpubr::ggpar(grap2, main = Main[i], subtitle = paste("data :", deparse(substitute(object2)), redmet),
                       legend.title = "Organelles", legend.position = "top",
                       ggtheme = ggplot2::theme_minimal(), xlab = axname[1], ylab = axname[2],
                       caption = capt) +
        gghighlight::gghighlight(p2$cond == all_cond[i], use_direct_label = FALSE,
                    unhighlighted_params = list(alpha = 0.1))  +
        ggplot2::guides(fill = ggplot2::guide_legend(override.aes =  list(size = 6)), size = FALSE) +
        ggplot2::scale_fill_manual(values = c(PaletteWithoutGrey(p2$markers, 1)))


      Gr[[paste(i)]] <- graph1a #saving in a list the plot created in order to all plot them on one figure
    }
  }

  if (!highpr & !vect & unknow){ #same thing as before but with unknown assignment in the markers column
    Gr <- list()
    for (i in 1:k){
      print(i)

      graph1a <- ggpubr::ggpar(grap2, main = Main[i], subtitle = paste("data :", deparse(substitute(object2)), redmet),
                       legend.title = "Organelles", legend.position = "top",
                       ggtheme = ggplot2::theme_minimal(), xlab = axname[1], ylab = axname[2]) +
        gghighlight::gghighlight(p2$cond == all_cond[i], use_direct_label = FALSE,
                    unhighlighted_params = list(alpha = 0.1))  +
        ggplot2::guides(fill = ggplot2::guide_legend(override.aes =  list(size = 6)), size = FALSE) +
        ggplot2::scale_fill_manual(values = c(PaletteWithoutGrey(p2$markers), "unknown" = "darkgrey")) +
        ggplot2::scale_size_manual(values = c(matlab::linspace(5,5,length(unique(p2$markers)) -1), "unknown" = 2))


      Gr[[paste(i)]] <- graph1a
    }

  }


  if (highpr & !vect){ #same thing as before but with highlighting sepecifics proteins
    Gr <- list()
    for (i in 1:k){
      print(i)

      graph1a <- ggpubr::ggpar(grap2, main = Main[i], subtitle = paste("data :", deparse(substitute(object2)), redmet),
                       legend.title = "Organelles", legend.position = "top",
                       ggtheme = ggplot2::theme_minimal(), xlab = axname[1], ylab = axname[2],
                       caption = capt) +
        gghighlight::gghighlight(p2$cond == all_cond[i] & is.na(match(p2$unprot, proteins)) == FALSE,
                    unhighlighted_params = list(alpha = 0.1),
                    label_key = prot,
                    label_params = list(size = 4, fill = "lightblue", color = "red"))  +
        ggplot2::scale_fill_manual(values = c(PaletteWithoutGrey(p2$markers, 1))) +
        ggplot2::geom_blank(aes(fill = markers), p2) +  #keep same color
        ggplot2::guides(fill = ggplot2::guide_legend(override.aes = list(size = 5))) #keep the legend

      Gr[[paste(i)]] <- graph1a
    }
  }



  if (vect){ #see the moving proteins with vector
    print("Vect = TRUE, data process")

    p2_f <- p2 %>% #filter the proteins we selected
      group_by(unprot) %>%
      filter(!is.na(match(unprot, proteins))) %>%
      ungroup()

    #creation of a data.frame which contains the coordinates from each proteins
    p2_f1 <- p2_f[,c("Dim1", "Dim2", "cond")]
    p2_f1$cond[1:(length(p2_f$cond)/2)] <- 1
    p2_f1$cond[(length(p2_f$cond)/2 + 1):length(p2_f$cond)] <- 2
    p2_f1 <- data.frame(p2_f1[which(p2_f1$cond == 1), "Dim1"],  p2_f1[which(p2_f1$cond == 2), "Dim1"],
                        p2_f1[which(p2_f1$cond == 1), "Dim2"], p2_f1[which(p2_f1$cond == 2), "Dim2"])
    colnames(p2_f1) <- c("xd", "xe", "yd", "ye")
    p2_f1$prot <- unique(p2_f$unprot)

    print("2. Building the plot")

    #creation of the plot : highlight the specifics proteins we chosed + the vectors showing their shifting
    graph_v <- ggpubr::ggpar(grap2, main = Title, subtitle = paste("data :", deparse(substitute(object2)), redmet),
                     legend.title = "Organelles", legend.position = "top",
                     ggtheme = ggplot2::theme_minimal(), xlab = paste("Axe", ax[1]), ylab = paste("Axe", ax[2]),
                     caption = capt) +
      gghighlight::gghighlight(!is.na(match(p2$unprot, proteins)),                                         #highlight proteins
                  unhighlighted_params = list(alpha = 0.1),
                  label_key = prot,
                  label_params = list(size = 4, fill = "lightblue", color = "red"))  +
      ggplot2::geom_segment(aes(x = xd, y = yd, xend=xe, yend=ye, color = prot), data = p2_f1,                        #plot the vectors
                   inherit.aes = FALSE, arrow = arrow(length = unit(0.05, "npc")),
                   size = 1.2, alpha = 0.6) +
      ggplot2::scale_fill_manual(values = c(PaletteWithoutGrey(p2$markers, 1))) +
      ggplot2::geom_blank(aes(fill = markers), p2) +  #keep same color
      ggplot2::guides(fill = ggplot2::guide_legend(override.aes = list(size = 5))) #keep the legend

    print("3. End of plot building")
  }

  if (vect){
    if (expor.png){ #if you want to export the plot created in png of pdf format
      ggplot2::ggsave(paste0(deparse(substitute(object2)), "_", redmet, "_vect.png"), plot = graph_v, width = 2000/300, height = 1000/300)
    }
    if (expor.pdf){
      ggplot2::ggsave(paste0(deparse(substitute(object2)), "_", redmet, "_vect.pdf"), plot = graph_v, width = 2000/300, height = 1000/300)
    }
    return(graph_v) #return th plot with the vectors
  }

  else { #when vect = FALSE, we can have several plot --> several condition
    #here we want to show all the plots on one figure in order to compare them (if there is several)
    if (length(Gr) == 1){ #this lines are for define the number of columns and rows in ggarrange according the number of plot = the number of conditions
      nR <- 1
      nC <- 1
    }
    else{
      if (length(Gr) == 2){
        nR <- 1
        nC <- 2
      }
      else{
        nR <- 2
        nC <- length(Gr)/2
      }
    }

    print("Building the final plot")
    #put all the plots created in one figure
    final_graph <- ggpubr::ggarrange(plotlist = Gr, nrow = nR, ncol = nC,
                             common.legend = TRUE, legend = "right")


    if (expor.png){ #export the figure in png or pdf format
      ggplot2::ggsave(paste0(deparse(substitute(object2)), "_", redmet, "_map.png"), plot = final_graph, width = 2000/72, height = 1000/72)
    }
    if (expor.pdf){
      ggplot2::ggsave(paste0(deparse(substitute(object2)), "_", redmet, "_map.pdf"), plot = final_graph, width = 2000/72, height = 1000/72)
    }

    return(final_graph)   #return the figure
  }
}
