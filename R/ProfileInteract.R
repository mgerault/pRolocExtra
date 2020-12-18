#' @title Visualize interactive protein profile
#'
#' @description
#' \code{ProfileInteract} allows you to visualize protein profile from the organelle(s) you want
#' interactively or not, and highlight specific proteins.
#' It is based on the plotDist function from pRoloc package.
#'
#' @param obj A \code{\link{MSnSet}} object
#' @param mrk The name of the column which contains the markers from the data
#' @param Organelle A character vector which contains the organelle(s) that you want to see the protein profile
#' @param Interact A logical argument to tell if you want an interactive plot
#' @param mytitle A logical argument to tell if you want to change the first part of the subtitle (name of the data)
#' @param TITLE A character which be your title (use with mytitle = TRUE)
#' @param Clust A logical argument to tell if you use clustered data.
#' (mrk should be equal to "knn" or orther clustering algorithm or any name of your data).
#' @param one_pr A logical argument to tell if you want to highlight a specific protein profile (or several)
#' @param protein A character vector which contains the proteins you want to highlight.
#' if NULL, show all the proteins
#' @param xLab A character which will be the name of the x-axis
#' @param yLab A character which will be the name of the y-axis
#'
#' @return A gpplot or ggplotly object : the protein profiles
#'
#' @seealso \code{\link{ggplotly}}, \code{\link{gghighlight}} and \code{\link{plotDist}} from Proloc package
#'
#' @export
#'
#' @examples
#' library(pRolocdata)
#' data(tan2009r1)
#' Profil_interact(tan2009r1)
ProfileInteract <- function(obj, mrk = "markers", Organelle = "Golgi", Interact = TRUE,
                            mytitle = FALSE, TITLE = "Profile", Clust = FALSE,
                            one_pr = FALSE, protein = NULL, xLab = "Fractions", yLab = "Intensity") {

  #check the fractions names (from plotDist from pRoloc package)
  if (is.character(MSnbase::sampleNames(obj)) & length(MSnbase::sampleNames(obj)) == 1) {
    if (sum(MSnbase::sampleNames(obj) %in% names(MSnbase::pData(obj))) != 1)
      stop("'fractions' must be a single pData name.")
    MSnbase::sampleNames(obj) <- as.character(MSnbase::pData(obj)[, MSnbase::sampleNames(obj)])
  }

  all_mrk <- as.character(unique(MSnbase::fData(obj)[[mrk]]))  #get all the markers from the column selected

  if (Clust){
    if (!one_pr){
      #get the number of proteins of which we know their localization (the training set for each organelle)
      npr_kn <- length(which(MSnbase::fData(obj)$markers == Organelle))
    }
    else {
      npr_kn <- length(which(MSnbase::fData(obj)$markers != "unknown"))      #same thing but for all organelle
    }

    prtx <- paste("proteins", "(based on", npr_kn, "proteins)")     #will be used for the graph title
  }
  else {
    prtx <- "proteins from known location"
  }

  if (!mytitle){ #put your own title
    subt <- deparse(substitute(obj)) #if not, take the name of the data (the obj argument)
  }
  else {
    subt <- TITLE
  }


  if (!one_pr){
    i <- which(!is.na(match(MSnbase::fData(obj)[[mrk]], Organelle)))  #find the proteins that are in the organelle
    npr <- length(i)                                         #number of protein in this organelle
    m_tit <- paste(Organelle, collapse = ", ")               #will be used for the main graph title
  }
  else {
    i <- which(!is.na(match(MSnbase::fData(obj)[[mrk]], all_mrk)))    #same thing but for all organelle
    npr <- length(which(fData(obj)$markers != "unknown"))
    m_tit <- "all organelle"
  }

  .data <- as.data.frame(t(MSnbase::exprs(obj)[i,]))             #the data from each protein we selected
  .data$fraction <- rownames(.data)                     #the row are the fraction, le column the proteins
  .data <- .data %>%                                    #reorganise the data in order to plot them : Intensity versus Fraction
    dplyr::select(fraction, colnames(.data)) %>%
    tidyr::gather(key = "prot", value = "value", -fraction)
  .data$mark <- fData(obj)[mrk][.data$prot,]            #get the organelle for each proteins


  if (Interact | one_pr){
    g <- ggplot2::ggplot(.data, aes(x=fraction, y=value)) +
      ggplot2::geom_line(aes(color = prot, linetype=mark, group = prot)) +  #we can differentiate every protein according to their organelle (only in interactive graph)
      labs(linetype = "Organelle", color = "Proteins")
  }

  if (!Interact & !one_pr){
    g <- ggplot(.data, aes(x=fraction, y=value)) +
      geom_line(aes(color = mark, group = prot)) + scale_color_manual(values = PaletteWithoutGrey(.data$mark,1)) + #differentiate every organelle, set color for each marker with PaletteWithoutGrey
      ggplot2::labs(color = "Organelle")
  }

  #change the title, subtitle and axis names
  g <- g +
    ggplot2::ggtitle(label= paste("Profile of", m_tit),
            subtitle =  paste(subt,",", npr, prtx)) + ggplot2::theme(plot.title = element_text(hjust = 0.5),
                                                            plot.subtitle = element_text(hjust = 0.5)) +
    ggplot2::xlab(xLab) + ggplot2::ylab(yLab)


  if (one_pr & !is.null(protein)){ #highlight specific proteins
    if (!Interact) {
      g <- g + gghighlight::gghighlight(!is.na(match(prot, protein)), #the predicate
                           label_key = prot,
                           max_highlight = 50,
                           label_params = list(size = 4, fill = "lightblue", color = "red", alpha =1),
                           unhighlighted_params = list(size = 0.5)) +
        ggplot2::aes(alpha=0.8, size = 1.4) + ggplot2::guides(size = FALSE, alpha = FALSE)  #bigger size for profile we want to highlight
    }
    if (Interact){
      g <- g + gghighlight::gghighlight(!is.na(match(prot, protein)),
                           label_key = prot,
                           max_highlight = 50,
                           use_direct_label = FALSE,    #ggplotly (for now) can't handle labeling
                           unhighlighted_params = list(size = 0.5))
    }

  }

  if (Interact){
    gI <- plotly::ggplotly(g, tooltip = c("color", "linetype", "group")) %>%  #the tooltip is the legend which is shown by the cursor
      plotly::layout(hovermode = "closest", clickmode = "select",
             margin = list(t = 75),
             title = list(text = paste0(paste("Profile of", m_tit),     #personalize the title and subtitle
                                        '<br>',
                                        '<sup>',
                                        paste(subt,",", npr, prtx),
                                        '</sup>'

             ))
      )


    for (i in 1:length(gI$x$data)){ #change the name of the element of the legend and how they are grouped (here by their organelle)
      if (!is.null(gI$x$data[[i]]$name)){
        gI$x$data[[i]]$name = gsub("\\)", "", gsub("\\(","",gI$x$data[[i]]$name))
        gI$x$data[[i]]$legendgroup = gsub("\\)","", stringr::str_split(gI$x$data[[i]]$legendgroup,",")[[1]][2])
      }
      if (gI$x$data[[i]]$showlegend == FALSE){
        gI$x$data[[i]]$showlegend <- TRUE      #show every legend
      }
    }
  }


  if (Interact){
    return(gI) #the interactive graph
  }
  else {
    return(g)
  }

}
