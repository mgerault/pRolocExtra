#' @title PCA visualization of your \code{\link{MSnSet}} data
#'
#' @description
#' \code{datavisupca} allows you to visualize the PCA plot of your data, clustered and not clustered on the same figure.
#' You can also choose to not see the plot and only get back your clustered data.
#' The clustering method available are from the pRoloc package (Laurent Gatto and al.).
#'
#' @param object A \code{\link{MSnSet}} object
#' @param mfcol The name of the column which contains the markers from your data
#' @param method The clustering method, available : svm, ksvm, knn, perTurbo, nnet (neural network),
#'               rf (random forest), naiveBayes, xgboost, CPA (constrained proportionate assignment),
#'               CNN or SpatialTransformer.
#' @param ax A numeric vector of length 2, the axes on which you want to see the PCA plot (depend of the number of fraction of the data)
#' @param sh.gr A logical argument, to show or not the plot.
#' if FALSE, return only a \code{\link{MSnSet}} object : your data + the clustering information
#' @param tm An integer corresponding to the times parameter of clustering optimization functions from pRoloc package
#' @param cv An integer corresponding to the cross validation parameter of clustering optimization functions from pRoloc package
#'
#' @return A list containing the two PCA plots on the same figure (clustered and not clustered) and the clustered MSnSet object if sh.gr = TRUE
#' else, it return only the clustered \code{\link{MSnSet}} object
#'
#' @seealso \code{\link{svmOptimisation}} from Roloc package and \code{\link{fviz_pca_ind}} from factoextra package
#'
#' @export
#'
#' @examples
#'
#' library(pRolocExtra)
#' datavisupca(tan2009r1)

datavisupca <- function(object, mfcol = "markers", method = "knn",
                        ax = c(1,2), sh.gr = TRUE, tm = 5, cv = 5){

  object1 <- nainunknown(object, fcol = mfcol)                #take off the NA values from the markers column and replace it by "unknown"
  object2 <- nafilterdata(object1)                            #take off proteins with NA or infinite values from the data

  if (sh.gr){
    rep.pca <- FactoMineR::PCA(MSnbase::exprs(object2), graph = FALSE)

    myalpha <- rep(0,length(MSnbase::fData(object2)[[mfcol]]))  #creation of a low alpha for the markers that are unknown (better visualization on the plot)
    id <- which(MSnbase::fData(object2)[[mfcol]] == "unknown")
    myalpha[id] <- 0.2
    myalpha[-id] <- 0.8

    rep.fgr <- factoextra::fviz_pca_ind(rep.pca, geom.ind = "point",                           #build the first PCA graph
                            pointshape = 21,
                            pointsize = MSnbase::fData(object2)[[mfcol]],                   #we define the size in order to change it manually
                            fill.ind = MSnbase::fData(object2)[[mfcol]],                    #we color the protein according their organelle
                            alpha.ind = myalpha, mean.point = FALSE, axes = ax) +
      ggplot2::guides(fill = ggplot2::guide_legend(override.aes = list(size = 6)), size = "none")     #change the size component from the fill legend and not showing the size legend

    graph1 <- rep.fgr +
      ggplot2::labs(title = "PCA", subtitle = paste("data :", deparse(substitute(object))),
                    fill = "Organelles") +
      ggplot2::theme_bw() +
      ggplot2::scale_fill_manual(values = PaletteWithoutGrey(MSnbase::fData(object2)[[mfcol]])) +            #change the color manually with PaletteWithoutGrey function and set the color of protein which are assign to unknown in grey
      ggplot2::scale_size_manual(values = unlist(lapply(unique(MSnbase::fData(object2)[[mfcol]]), function(x){
                                                        if(x == "unknown"){y <- 2}
                                                        else{y <- 5}
                                                        names(y) <- x;
                                                        y
                                                      })
                                                      )
                                 ) #change the size manually : the protein which are assign to unknown are smaller than the others
  }

  #clustering method from the pRoloc package; some have specifics argument --> pRoloc documentation
  #first, optimizes the hyperparameters of the clustering algorithm according to the data
  #then cluster the data, on the starting space not the one from the PCA
  if (method == "knn") {
    paramsrep <- pRoloc::knnOptimisation(object2, k = seq(2, 20, 2),
                                 times = tm, fcol=mfcol, xval = cv)

    object2 <- pRoloc::knnClassification(object2, paramsrep, fcol=mfcol)
  }
  else if (method == "svm") {
    paramsrep <- pRoloc::svmOptimisation(object2, times = tm, fcol=mfcol, xval = cv)

    object2 <- pRoloc::svmClassification(object2, paramsrep, fcol=mfcol)
  }
  else if (method == "ksvm") {
    paramsrep <- pRoloc::ksvmOptimisation(object2, times = tm, fcol=mfcol, xval = cv)

    object2 <- pRoloc::ksvmClassification(object2, paramsrep, fcol=mfcol)
  }
  else if (method == "rf") {
    paramsrep <- pRoloc::rfOptimisation(object2, times = tm, fcol=mfcol, xval = cv)

    object2 <- pRoloc::rfClassification(object2, paramsrep, fcol=mfcol, mtry = c(2,5,10))
  }
  else if (method == "naiveBayes") {
    paramsrep <- pRoloc::nbOptimisation(object2, times = tm, fcol=mfcol, xval = cv)

    object2 <- pRoloc::nbClassification(object2, paramsrep, fcol=mfcol)
  }
  else if (method == "nnet") {
    paramsrep <- pRoloc::nnetOptimisation(object2, times = tm, fcol=mfcol, xval = cv)

    object2 <- pRoloc::nnetClassification(object2, paramsrep, fcol=mfcol)
  }
  else if (method == "perTurbo") {
    paramsrep <- pRoloc::perTurboOptimisation(object2, times = tm, fcol=mfcol, xval = cv,
                                      pRegul = 10^seq(-1,0,0.2), sigma = 10^seq(-1,1,0.5),
                                      inv = "Inversion Cholesky",
                                      reg = "tikhonov")

    object2 <- pRoloc::perTurboClassification(object2, paramsrep, fcol=mfcol)
  }
  else if (method == "xgboost") {
    paramsrep <- xgboostOptimisation(object2, times = tm, fcol=mfcol, xval = cv)

    object2 <- xgboostClassification(object2, paramsrep, fcol=mfcol)
  }
  else if (method == "CPA") {
    object2 <- CPAClassification(object2, fcol=mfcol)
  }
  else if (method == "CNN") {
    object2 <- CNNimgClassification(object2, fcol=mfcol)
  }
  else if (method == "SpatialTransformer") {
    object2 <- CNNimgClassification(object2, fcol=mfcol, model_type = "SpatialTransformer")
  }

  method_name <- method
  if(method_name == "SpatialTransformer"){
    method_name <- "CNN"
  }
  if (sh.gr){
    #a score is calculated for each protein clustered when the clustering algorithm is called
    ptszerep <- exp(MSnbase::fData(object2)[[paste0(method_name, ".scores")]]) - 1        #get back this score and modify it in order to change size of each point on the plot according to their clustering score
    rep.loss <- round(mean(MSnbase::fData(object2)[[paste0(method_name, ".scores")]]),3)  #get the mean of the clustering score to have it on the plot

    #exact same step as before but with the new assignment, the clustering score, etc.
    rep.clus <- factoextra::fviz_pca_ind(rep.pca, geom.ind = "point",
                             pointshape = 21,
                             pointsize = ptszerep,
                             fill.ind = MSnbase::fData(object2)[[method_name]],
                             alpha.ind = 0.6, mean.point = FALSE, axes = ax) +
      ggplot2::guides(fill = ggplot2::guide_legend(override.aes = list(size = 5)), size = "none")


    graph2 <- rep.clus +
      labs(title = method, subtitle = paste("data :", deparse(substitute(object))),
           fill = "Organelles", caption = paste("Scores =", rep.loss)) +
      theme_bw() +
      ggplot2::scale_fill_manual(values = PaletteWithoutGrey(MSnbase::fData(object2)[[method_name]]))

    #return the plot and the clustered object under a list
    return(list("graph" = ggpubr::ggarrange(graph1, graph2, nrow = 1, ncol = 2, common.legend = TRUE), "MSn" = object2))
  }
  else{
    #return the clustered object only
    return(object2)
  }
}
