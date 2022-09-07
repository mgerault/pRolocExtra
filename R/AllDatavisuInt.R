#' @title Interactive (or not) data visualization
#'
#' @description
#' \code{AllDatavisuInt} allows you to visualize your data with many options.
#' You can choose the reduction method (PCA, t-SNE, ...), the clustering method, export the plot as a png or pdf file.
#' The function adapts according to the data, so if you have replicates or two different experiment conditions you can visualize it.
#' This is very useful when you work with dynamic organellar map.
#' With this type of data you can visualize the movement of protein between the two conditions with a vector.
#' The two many differences between \code{AllDatavisuInt} and \code{\link{AllDatavisuInt}}, are that you can visualize your
#' interactively with plotly, and it shows only one graph per figure.
#' This function was mainly created for shiny purpose. (for example run \code{rnShinyVisualization})
#'
#' @param object2 A \code{\link{MSnSet}} object
#' @param redmet The reduction method from the pRoloc package : "PCA", "t-SNE", "MDS", "nipals", "lda", "kpca", and "umap".
#' @param cmet the name of the column of your data which contains the markers, if it contains "unknown" assignment, please precise unknow = TRUE
#' @param ax A numeric vector of length two, you can choose on which axes you want to see the plot.
#' (depend of the number of fraction of the data)
#' @param Interact A logical argument to tell if you want interactive plot
#' @param Mean_point A logical argument to tell if you want to print the mean point of each cluster on the plot
#' @param Title The title of the plot, use only when vect = TRUE
#' @param name_cond A character argument, it is the name of the column of your data which contains the condition of the experiment (control or treated).
#'  You can have replicates.
#'  If you have no dynamic experiment, don't bother with this parameter.
#' @param highpr A logical argument to highlight or not specific protein
#' @param proteins A character vector containing the proteins you want to highlight
#' @param vect A logical argument to use when you have dynamic experiment.
#'  If TRUE, allows you to see the proteins movement between the two conditions with vectors.
#'  You also have to choose the proteins you want to see with the proteins argument.
#' @param highcond A logical argument if you want to highlight specific condition (replicates / control and treated)
#' @param Condition A numeric argument to tell which condition you want to highlight
#' @param expor.png A logical argument to export the figure in a png file
#' @param expor.pdf A logical argument to export the figure in a pdf file
#' @param yourseed An integer for the t-SNE algorithm in order to having same plot if there is several
#' @param Source A character specifying the source of the plotly graph (use when Interact = TRUE)
#' @param mysubtitle A logical argument to tell if you want to put your own subtitle
#' @param subtitle A character vector which is your subtitle
#'
#' @return A figure showing your data, depending of your chosen parameters
#'
#' @seealso \code{\link{plot2D}}, \code{\link{ggplotly}} and \code{\link{svmOptimisation}} from pRoloc package for more details
#'
#' @export
#'
#' @examples
#'
#' library(pRolocExtra)
#' tan2009r1_clustered <- datavisupca(tan2009r1, method ="knn", sh.gr = FALSE)
#' AllDatavisuInt(tan2009r1_clustered, cmet = "knn", Interact = TRUE)

AllDatavisuInt <- function(object2, redmet = "PCA", cmet = "svm", ax = c(1,2), Interact = FALSE,
                        Title = "Cellular map", name_cond = "cond", Mean_point = FALSE,
                        highpr = FALSE, proteins = c("PKN2", "GRB2", "SHC1", "EGFR"),
                        vect = FALSE, highcond = FALSE, Condition = 1,
                        expor.png = FALSE, expor.pdf = FALSE, yourseed = 500, Source = "AA",
                        mysubtitle = FALSE, subtitle = ""){

  set.seed(yourseed) #set the seed for the t-SNE

  message("1. Create the data plot")
  if (redmet == "umap"){
    p2 <- uwot::umap(exprs(object2))
    p2 <- as.data.frame(p2)
    rownames(p2) <- rownames(exprs(object2))
    colnames(p2) <- c("Dimension 1", "Dimension 2")
  }
  else{
    p2 <- pRoloc::plot2D(object2, fcol = cmet, method = redmet,
                         plot = FALSE,  #get the data from the plot2D function of pRoloc package (see documentation of pRoloc)
                         dims = ax)
    p2 <- as.data.frame(p2)
  }

  axname <- colnames(p2)  #save the names of the axes


  message("2. Data process")
  if (redmet != "t-SNE"){#t-SNE only work with unique features --> filter first
    p2$markers <- MSnbase::fData(object2)[[cmet]]
  }
  else {
    if (length(unique(object2@assayData$exprs)[,1]) == length(MSnbase::fData(object2)[[cmet]]))
      p2$markers <- MSnbase::fData(object2)[[cmet]]

    else
      p2$markers <- MSnbase::fData(object2)[[cmet]][-which(duplicated(object2@assayData$exprs) == TRUE)]
  }

  p2$markers <- as.factor(p2$markers)
  unknow <- "unknown" %in% p2$markers

  if (!unknow){ #set the size of the points on the plot according their clustering score or not
    rep.loss <- round(mean(MSnbase::fData(object2)[[paste0(cmet, ".scores")]]),3)
    capt <- paste("The score of", cmet, "clustering method is", rep.loss)

      if (redmet != "t-SNE")
        p2$score <- as.factor(exp(MSnbase::fData(object2)[[paste0(cmet, ".scores")]]) - 1)
      else {
        if (length(unique(object2@assayData$exprs)[,1]) == length(MSnbase::fData(object2)[[cmet]]))
          p2$score <- as.factor(exp(MSnbase::fData(object2)[[paste0(cmet, ".scores")]]) - 1)

        else
          p2$score <-  as.factor(exp(MSnbase::fData(object2)[[paste0(cmet, ".scores")]][-which(duplicated(object2@assayData$exprs) == TRUE)]) - 1)
      }
  }
  else { #if unknow = TRUE, no clustering algorithm was performed, so to change the size of the unknown features on the plot we have to define it
    p2$score <- p2$markers
    capt <- NULL
  }

  colnames(p2) <- c("Dim1", "Dim2", "markers", "score")    #reorganize the data, adding new columns (save the proteins names)
  p2$prot <- rownames(p2)
  if (is.null(MSnbase::fData(object2)$gene_name)){
    p2$unprot <- p2$prot
  }
  else {
    p2$unprot <- MSnbase::fData(object2)$gene_name
  }

  if (!is.null(MSnbase::fData(object2)[[name_cond]])){
    if (redmet != "t-SNE")
      p2$cond <- as.factor(MSnbase::fData(object2)[[name_cond]])
    else {
      if (length(unique(object2@assayData$exprs)[,1]) == length(MSnbase::fData(object2)[[cmet]]))
        p2$cond <- as.factor(MSnbase::fData(object2)[[name_cond]])

      else
        p2$cond <-  as.factor(MSnbase::fData(object2)[[name_cond]][-which(duplicated(object2@assayData$exprs) == TRUE)])
    }
  }

  if(Mean_point){
    mrkdep2 <- as.character(unique(p2$markers))

    mean_p2_1 <- c()
    mean_p2_2 <- c()

    for (i in 1:length(mrkdep2)){
      mrk_p2 <- p2[which(p2$markers == mrkdep2[i]),]
      mean_p2_1 <- append(mean_p2_1, mean(mrk_p2$Dim1))
      mean_p2_2 <- append(mean_p2_2, mean(mrk_p2$Dim2))
    }

    mean_p2 <- data.frame(Dim1 = mean_p2_1, Dim2 = mean_p2_2, markers = mrkdep2)
  }
  message("3. Building the base plot")

  grap2 <- ggplot2::ggplot(data = p2, aes(x=Dim1, y=Dim2, fill=markers,
                                         size = score)) +
    ggplot2::geom_point(alpha = 0.7, shape = 21, aes(key = prot)) +
    ggplot2::guides(alpha = "none", size = "none")

  if(mysubtitle) {
    sub <- subtitle
  }
  else{
    sub <- deparse(substitute(object2))
  }

  message("4. Building the plot(s)")
  if (!Interact) {
    if (!highpr & !vect & !highcond){ #creating the plot according to the number of condition, then highlight the condition
      if (!unknow){
        graph1a <- grap2 +
          ggplot2::labs(title = Title, subtitle = paste("data :", sub, redmet),
                        x = axname[1], y = axname[2], caption = capt, fill = "Organelles") +
          ggplot2::theme_bw() +
          ggplot2::guides(fill = ggplot2::guide_legend(override.aes =  list(size = 6)), size = "none") +
          ggplot2::scale_fill_manual(values = PaletteWithoutGrey(p2$markers))
      }
      else {
        graph1a <- grap2 +
          ggplot2::labs(title = Title, subtitle = paste("data :", sub, redmet),
                        x = axname[1], y = axname[2], caption = capt, fill = "Organelles") +
          ggplot2::theme_bw() +
          ggplot2::guides(fill = ggplot2::guide_legend(override.aes =  list(size = 6)), size = "none") +
          ggplot2::scale_fill_manual(values = PaletteWithoutGrey(p2$markers)) +
          ggplot2::scale_size_manual(values = unlist(lapply(unique(p2$markers), function(x){
                                         if(x == "unknown"){y <- 2}
                                         else{y <- 5}
                                         names(y) <- x;
                                         y
                                         })
                                         )
                                     )
      }
    }

    if (!highpr & !vect & highcond){ #creating the plot according to the number of condition, then highlight the condition
      if (!unknow) {
        graph1a <- grap2 +
          ggplot2::labs(title = Title, subtitle = paste("data :", sub, redmet),
                        x = axname[1], y = axname[2], caption = capt, fill = "Organelles") +
          ggplot2::theme_bw() +
          gghighlight::gghighlight(p2$cond == Condition, use_direct_label = FALSE,
                                   unhighlighted_params = list(alpha = 0.1))  +
          ggplot2::guides(fill = ggplot2::guide_legend(override.aes =  list(size = 6)), size = "none") +
          ggplot2::scale_fill_manual(values = PaletteWithoutGrey(p2$markers))
      }

      else {
        graph1a <- grap2 +
          ggplot2::labs(title = Title, subtitle = paste("data :", sub, redmet),
                        x = axname[1], y = axname[2], caption = capt, fill = "Organelles") +
          ggplot2::theme_bw() +
          gghighlight::gghighlight(p2$cond == Condition, use_direct_label = FALSE,
                                   unhighlighted_params = list(alpha = 0.1))  +
          ggplot2::guides(fill = ggplot2::guide_legend(override.aes =  list(size = 6)), size = "none") +
          ggplot2::scale_fill_manual(values = PaletteWithoutGrey(p2$markers)) +
          ggplot2::scale_size_manual(values = unlist(lapply(unique(p2$markers), function(x){
                                                     if(x == "unknown"){y <- 2}
                                                     else{y <- 5}
                                                     names(y) <- x;
                                                     y
                                                     })
                                                     )
                                     )
      }

    }
    if (highpr & !vect & !highcond){ #same thing as before but with highlighting sepecifics proteins
      if (!unknow) {
        graph1a <- grap2 +
          ggplot2::labs(title = Title, subtitle = paste("data :", sub, redmet),
                        x = axname[1], y = axname[2], caption = capt, fill = "Organelles") +
          ggplot2::theme_bw() +
          gghighlight::gghighlight(is.na(match(p2$prot, proteins)) == FALSE,
                                   unhighlighted_params = list(alpha = 0.1),
                                   label_key = prot,
                                   label_params = list(size = 4, fill = "lightblue", color = "red"))  +
          ggplot2::scale_fill_manual(values = PaletteWithoutGrey(p2$markers)) +
          ggplot2::geom_blank(aes(fill = markers), p2) +  #keep same color
          ggplot2::guides(fill = ggplot2::guide_legend(override.aes = list(size = 5))) #keep the legend
      }
      else {
        graph1a <- grap2 +
          ggplot2::labs(title = Title, subtitle = paste("data :", sub, redmet),
                        x = axname[1], y = axname[2], caption = capt, fill = "Organelles") +
          ggplot2::theme_bw() +
          gghighlight::gghighlight(is.na(match(p2$prot, proteins)) == FALSE,
                                   unhighlighted_params = list(alpha = 0.1),
                                   label_key = prot,
                                   label_params = list(size = 4, fill = "lightblue", color = "red"))  +
          ggplot2::scale_fill_manual(values = PaletteWithoutGrey(p2$markers)) +
          ggplot2::scale_size_manual(values = unlist(lapply(unique(p2$markers), function(x){
                                                            if(x == "unknown"){y <- 2}
                                                            else{y <- 5}
                                                            names(y) <- x;
                                                            y
                                                            })
                                                     )
                                     ) +
          ggplot2::geom_blank(aes(fill = markers), p2) +  #keep same color
          ggplot2::guides(fill = ggplot2::guide_legend(override.aes = list(size = 5))) #keep the legend
      }
    }

    if (vect){ #see the moving proteins with vector
      message("Vect = TRUE, data process")

      try(if(is.null(MSnbase::fData(object2)[[name_cond]])) stop("Your data doesn't contain a condition columns"))

      p2_f <- p2[which(!is.na(match(p2$unprot, proteins))),]

      #creation of a data.frame which contains the coordinates from each proteins
      p2_f1 <- p2_f[,c("Dim1", "Dim2", "cond")]
      p2_f1$cond[1:(length(p2_f$cond)/2)] <- 1
      p2_f1$cond[(length(p2_f$cond)/2 + 1):length(p2_f$cond)] <- 2
      p2_f1 <- data.frame(p2_f1[which(p2_f1$cond == 1), "Dim1"],  p2_f1[which(p2_f1$cond == 2), "Dim1"],
                          p2_f1[which(p2_f1$cond == 1), "Dim2"], p2_f1[which(p2_f1$cond == 2), "Dim2"])
      colnames(p2_f1) <- c("xd", "xe", "yd", "ye")
      p2_f1$prot <- unique(p2_f$unprot)

      message("2. Building the plot")

      #creation of the plot : highlight the specifics proteins we chosen + the vectors showing their shifting
      if (!unknow) {
        graph1a <- grap2 +
          ggplot2::labs(title = Title, subtitle = paste("data :", sub, redmet),
                        x = axname[1], y = axname[2], caption = capt, fill = "Organelles") +
          ggplot2::theme_bw() +
          gghighlight::gghighlight(!is.na(match(p2$unprot, proteins)),                                         #highlight proteins
                                   unhighlighted_params = list(alpha = 0.1),
                                   label_key = prot,
                                   label_params = list(size = 4, fill = "lightblue", color = "red"))  +
          ggplot2::geom_segment(aes(x = xd, y = yd, xend=xe, yend=ye, color = prot), data = p2_f1,             #plot the vectors
                                inherit.aes = FALSE, arrow = arrow(length = unit(0.05, "npc")),
                                size = 1.2, alpha = 0.6) +
          ggplot2::scale_fill_manual(values = PaletteWithoutGrey(p2$markers)) +
          ggplot2::geom_blank(aes(fill = markers), p2) +  #keep same color
          ggplot2::guides(fill = ggplot2::guide_legend(override.aes = list(size = 5))) #keep the legend

      }
      else {
        graph1a <- grap2 +
          ggplot2::labs(title = Title, subtitle = paste("data :", sub, redmet),
                        x = axname[1], y = axname[2], caption = capt, fill = "Organelles") +
          ggplot2::theme_bw() +
          gghighlight::gghighlight(!is.na(match(p2$unprot, proteins)),                                         #highlight proteins
                                   unhighlighted_params = list(alpha = 0.1),
                                   label_key = prot,
                                   label_params = list(size = 4, fill = "lightblue", color = "red"))  +
          ggplot2::geom_segment(aes(x = xd, y = yd, xend=xe, yend=ye, color = prot), data = p2_f1,                        #plot the vectors
                                inherit.aes = FALSE, arrow = arrow(length = unit(0.05, "npc")),
                                size = 1.2, alpha = 0.6) +
          ggplot2::scale_fill_manual(values = PaletteWithoutGrey(p2$markers)) +
          ggplot2::scale_size_manual(values = unlist(lapply(unique(p2$markers), function(x){
                                                            if(x == "unknown"){y <- 2}
                                                            else{y <- 5}
                                                            names(y) <- x;
                                                            y   })
                                                     )
                                     ) +
          ggplot2::guides(fill = ggplot2::guide_legend(override.aes = list(size = 5))) #keep the legend

      }
    }
  }

  #Interactive graph
  else {
    if (!highpr & !vect & !highcond){ #creating the plot according to the number of condition, then highlight the condition
        graph1a <- grap2
    }
    if (!highpr & !vect & highcond){ #creating the plot according to the number of condition, then highlight the condition
      p2_filtered <- p2[which(p2$cond == Condition),]

        graph1a <- ggplot2::ggplot(p2_filtered) +
          ggplot2::geom_point(aes(x=Dim1, y=Dim2, fill = markers,
                         size = markers, key = prot), data = p2, colour = alpha("lightblue", 0.7), size = 0.9) +
          ggplot2::geom_point(aes(x=Dim1, y=Dim2, fill = markers,
                         size = markers, key = prot), shape= 21) +
          ggplot2::guides(alpha = "none", size = "none")
    }
    if (highpr & !vect & !highcond){ #same thing as before but with highlighting sepecifics proteins
      p2_filtered <- p2[proteins, ]

        graph1a <- ggplot2::ggplot(p2_filtered) +
          ggplot2::geom_point(aes(x=Dim1, y=Dim2, fill = markers,
                         size = markers, key = prot), data = p2, colour = alpha("lightblue", 0.7), size = 0.9) +
          ggplot2::geom_point(aes(x=Dim1, y=Dim2, fill = markers,
                         size = markers, key = prot), shape= 21) +
          ggplot2::guides(alpha = "none", size = "none")


    }
    if (vect){ #see the moving proteins with vector
      message("Vect = TRUE, data process")

      try(if(is.null(MSnbase::fData(object2)[[name_cond]])) stop("Your data doesn't contain a condition columns"))

      p2_f <- p2[which(!is.na(match(p2$unprot, proteins))),]

      #creation of a data.frame which contains the coordinates from each proteins
      p2_f1 <- p2_f[,c("Dim1", "Dim2", "cond")]
      p2_f1$cond[1:(length(p2_f$cond)/2)] <- 1
      p2_f1$cond[(length(p2_f$cond)/2 + 1):length(p2_f$cond)] <- 2
      p2_f1 <- data.frame(p2_f1[which(p2_f1$cond == 1), "Dim1"],  p2_f1[which(p2_f1$cond == 2), "Dim1"],
                          p2_f1[which(p2_f1$cond == 1), "Dim2"], p2_f1[which(p2_f1$cond == 2), "Dim2"])
      colnames(p2_f1) <- c("xd", "xe", "yd", "ye")
      p2_f1$prot <- unique(p2_f$unprot)

      message("2. Building the plot")

      #creation of the plot : highlight the specifics proteins we chosen + the vectors showing their shifting
        grap2 <- ggplot2::ggplot(p2_f) +
          ggplot2::geom_point(aes(x=Dim1, y=Dim2, fill = markers,
                         size = markers, key = prot), data = p2, colour = alpha("lightblue", 0.7), size = 0.9) +
          ggplot2::geom_point(aes(x=Dim1, y=Dim2, fill = markers,
                         size = markers, key = prot), shape= 21) + guides(alpha = "none", size = "none") +
          ggplot2::geom_segment(aes(x = xd, y = yd, xend=xe, yend=ye, color = prot), data = p2_f1,
                       inherit.aes = FALSE, arrow = arrow(length = unit(0.05, "npc")),
                       size = 1.2, alpha = 0.6)

    }
    if (unknow) {
      graph1a <- grap2 +
        ggplot2::labs(title = Title, subtitle = paste("data :", sub, redmet),
                      x = axname[1], y = axname[2], caption = capt, fill = "Organelles") +
        ggplot2::theme_bw() +
        ggplot2::scale_fill_manual(values = PaletteWithoutGrey(p2$markers)) +
        ggplot2::scale_size_manual(values = unlist(lapply(unique(p2$markers), function(x){
                                    if(x == "unknown"){y <- 2}
                                    else{y <- 5}
                                    names(y) <- x;
                                    y })
                                    )
                                   ) +
        ggplot2::guides(fill = ggplot2::guide_legend(override.aes = list(size = 5))) #keep the legend
    }
    else {
      graph1a <- grap2 +
        ggplot2::labs(title = Title, subtitle = paste("data :", sub, redmet),
                      x = axname[1], y = axname[2], caption = capt, fill = "Organelles") +
        ggplot2::theme_bw() +
        ggplot2::scale_fill_manual(values = PaletteWithoutGrey(p2$markers)) +  #keep same color
        ggplot2::guides(fill = ggplot2::guide_legend(override.aes = list(size = 5))) #keep the legend
    }

  }

  if (Mean_point) {
    graph1a <- graph1a + ggplot2::geom_point(aes(x=Dim1, y=Dim2, fill = markers),
                                             data = mean_p2, shape = 22, size = 5, show.legend = FALSE) +
      ggplot2::guides(shape = "none", size = "none")
  }
  if (Interact){
    Subt <- sub
    if (vect){
      if (redmet != "t-SNE"){
        graph1a <- plotly::ggplotly(graph1a, tooltip = c("fill", "key"), dynamicTicks = TRUE, source = Source) %>%
          plotly::layout(hovermode = "closest", clickmode = "event+select", dragmode = "select",
                         title = list(text = paste0(Title,
                                                    '<br>',
                                                    '<sup>',
                                                    paste("data :", Subt, redmet),
                                                    '<br>',
                                                    capt, '</sup>')),
                         annotations = ArrowAnnot(p2_f1)
          )
      }
      else{
        graph1a <- plotly::ggplotly(graph1a, tooltip = c("fill", "key"), dynamicTicks = TRUE, source = Source) %>%
          plotly::layout(hovermode = "closest", clickmode = "event+select", dragmode = "select",
                         title = list(text = paste0(Title,
                                                    '<br>',
                                                    '<sup>',
                                                    paste("data :", Subt, redmet),
                                                    '<br>',
                                                    capt, '</sup>')),
                         annotations = ArrowAnnot(p2_f1, d = 1)
          )
      }


      #change the legend generated by plotly : take off the parentheses and keep all legends element
      for (i in 1:length(graph1a$x$data)){
        if (!is.null(graph1a$x$data[[i]]$name)){
          graph1a$x$data[[i]]$name = gsub("\\(","", stringr::str_split(graph1a$x$data[[i]]$name,",")[[1]][1])
        }
        if (graph1a$x$data[[i]]$showlegend == FALSE){
          graph1a$x$data[[i]]$showlegend <- TRUE
        }
      }
    }
    else{
      graph1a <- plotly::ggplotly(graph1a, tooltip = c("fill", "key"), dynamicTicks = TRUE, source = Source) %>%
        plotly::layout(hovermode = "closest", clickmode = "event+select", dragmode = "select",
               title = list(text = paste0(Title,
                                          '<br>',
                                          '<sup>',
                                          paste("data :", Subt, redmet),
                                          '<br>',
                                          capt, '</sup>'))
               )


      #same thing as before for the plotly legend
      for (i in 1:length(graph1a$x$data)){
        if (!is.null(graph1a$x$data[[i]]$name)){
          graph1a$x$data[[i]]$name = gsub("\\(","", stringr::str_split(graph1a$x$data[[i]]$name,",")[[1]][1])
        }
      }
    }

    if (!unknow){
      #problem with the plotly legend : take also the size (size is according to the clustering score)
      #keep only the biggest point in the legend (only solution I found for now)
      mrk <- c()
      for (i in 1:(length(graph1a$x$data))){
        mrk <- append(mrk,graph1a$x$data[[i]]$name)
      }

      idnotl <- c()
      for (i in 1:length(unique(mrk))){
        a <- which(mrk == unique(mrk)[i])
        b <- which(a != max(a))
        if (Mean_point){
          c <- which(b != max(b))
          idnotl <- append(idnotl, a[c])
        }
        else
          idnotl <- append(idnotl, a[b])
      }

      for (i in idnotl){
        graph1a$x$data[[i]]$showlegend <- FALSE
      }
    }
  }

  message("Building the final plot")
  if (expor.png == TRUE){ #export the figure in png or pdf format
    ggplot2::ggsave(paste0(sub, "_", redmet, "_map.png"), plot = graph1a, width = 2000/72, height = 1000/72)
  }
  if (expor.pdf == TRUE){
    ggplot2::ggsave(paste0(sub, "_", redmet, "_map.pdf"), plot =graph1a, width = 2000/72, height = 1000/72)
  }

  return(graph1a)   #return the figure
}

