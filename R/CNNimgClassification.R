#' Classification using convolutionnal network on image profile.
#'
#' @title CNN img classification
#'
#' @param object An instance of class \code{"\linkS4class{MSnSet}"}.
#' @param scores One of \code{"prediction"}, \code{"all"} or
#'     \code{"none"} to report the score for the predicted class
#'     only, for all classes or none.
#' @param fcol The feature meta-data containing marker definitions.
#'     Default is \code{markers}.
#' @param resample_data Logical to tell if you want to rebalance the data (remove proteins with bad profiles and
#'                      add artificial ones to get a more uniform distribution of the number of profiles per class)
#' @param model_type The type of neural network you want (either 'CNN' or 'SpatialTransformer').
#' @param keep_img Logical to tell if you want to keep the images profiles created.
#' @param epochs The number of epochs for training; default is 40.
#' @param batch_size The size of the batch for training; default is 16.
#' @param lr The learning rate of the neural network; default is 0.001.
#' @param show_plot Logical to tell if you want to print the loss and accuracy plots.
#'
#' @details This algorithm will not apply clustering on raw data but on the image of the profiles from each protein.
#'          Then by using a convolutional network, it will learn to classify those images and then make predictions.
#'          Also, since the output is a log softmax, we can interpret these scores as probabilities. So we can handle
#'          proteins with composite profiles.
#'
#'          The first model available is a classical CNN containg two conv2d layers with dropout, Tanh activation and maxpool.
#'          The second one is a spatial transformers using convolutional layers and affine transformation.
#'
#'          This function calls python scripts in order to use torch python library, so it totally depends on \code{reticulate}
#'          r package.
#'
#' @return An instance of class \code{"\linkS4class{MSnSet}"} with
#'     \code{CNNimg} and \code{CNNimg.scores} feature variables storing the
#'     classification results and scores respectively.
#'
#' @export
#'
#' @examples
#' library(pRolocExtra)
#' data(tan2009r1)
#' res <- CNNimgClassification(tan2009r1)
#' getPredictions(res, fcol = "CNNimg")
#' getPredictions(res, fcol = "CNNimg", t = 0.75)
#' plot2D(res, fcol = "CNNimg")

CNNimgClassification <- function (object,
                                  scores = c("prediction", "all", "none"),
                                  fcol = "markers",
                                  resample_data = TRUE,
                                  model_type = "CNN",
                                  keep_img = FALSE,
                                  epochs = 40, batch_size = 32,
                                  lr = 0.001, show_plot = FALSE)
{
  if(sum(model_type %in% c("CNN", "SpatialTransformer")) == 0){
    stop("model_type can only be 'CNN' or 'SpatialTransformer'")
  }
  is_conda <- tryCatch(reticulate::conda_list(),
                            error = function(e)
                              "no_conda")
  if(class(is_conda) == "character"){
    stop("No conda installation was found in your computer.
         \nPlease reload the pRolocExtra package and answer yes to configure miniconda.")
  }
  for(pk in  c("matplotlib", "pandas", "cv2", "pytorch", "torchvision")){
    is_inst <- reticulate::py_module_available(stringr::str_remove(pk, "py"))
    if(!is_inst){
      if(stringr::str_detect(pk, "torch")){
        reticulate::py_install(pk, channel = c("pytorch")) ## torchvision is only accessible in channel pytorch
      }
      else if(pk == "cv2"){
        reticulate::py_install("opencv")
      }
      else{
        reticulate::py_install(pk)
      }
    }
  }

  scores <- match.arg(scores)

  file_name = paste0(format(Sys.time(), "%y%m%d_%H%M_"), deparse(substitute(object)), "_")

  # set train and test set
  trainInd <- which(fData(object)[, fcol] != "unknown")
  testInd <- which(fData(object)[, fcol] == "unknown")

  trainSet <- pRoloc:::subsetAsDataFrame(object, fcol, train = TRUE)
  testSet <- pRoloc:::subsetAsDataFrame(object, fcol, train = FALSE)
  write.csv(trainSet, paste0(file_name, "train.csv"))
  write.csv(testSet, paste0(file_name, "test.csv"))

  .trainSep <- pRoloc:::separateDataSet(trainSet, fcol)
  .testSep <- pRoloc:::separateDataSet(testSet, fcol)


  py_script <- reticulate::import_from_path("image_profile_classification", path = system.file("python", package = "pRolocExtra"))
  #reticulate::source_python(system.file("python", "image_profile_classification.py", package = "pRolocExtra"))
  model = py_script$ImgProfileClassification(paste0(file_name, "train.csv"),
                                             resample = resample_data,
                                             model_type = model_type,
                                             keep_img = keep_img,
                                             epoch = as.integer(epochs),
                                             batch_size = as.integer(batch_size),
                                             lr = lr, show_plot = show_plot)
  his = model$fit()

  #get predictions from CNNimg
  ans = model$pred(paste0(file_name, "test.csv"))
  predictedLabels <- unlist(lapply(ans, function(x) {y = list(x$markers)
                                                     names(y) <- x$name;
                                                     y})
                            )

  fData(object)$CNNimg <-  rep("", length(trainInd) + length(testInd))
  fData(object)$CNNimg[trainInd] <- as.character(.trainSep$theLabels)
  fData(object)[names(predictedLabels),]$CNNimg <- predictedLabels

  #get back scores from CNNimg
  if (scores == "all") {
    nbLabels <- length(levels(.trainSep$theLabels))
    tempScores <- matrix(rep(0, nbLabels * (length(trainInd) +
                                              length(testInd))), ncol = nbLabels)
    colnames(tempScores) <- names(ans[[1]]$score)
    rownames(tempScores) <- rownames(fData(object))

    for(i in 1:length(trainInd)){
      tempScores[trainInd[i], .trainSep$theLabels[i]] <- 1
    }
    i <- 1:length(testInd)
    sc <- lapply(ans, function(x) unlist(x$score))
    names(sc) <- unlist(lapply(ans, function(x) x$name))
    for(n in names(sc)){
      tempScores[n,] <- sc[[n]]
    }
    scoreMat <- tempScores
    colnames(scoreMat) <- paste0(colnames(scoreMat), ".CNNimg.scores")
    fData(object)$CNNimg.all.scores <- scoreMat
  }
  else if (scores == "prediction") {
    fData(object)$CNNimg.scores <- 1
    sc <- unlist(lapply(ans, function(x) max(unlist(x$score))))
    names(sc) <- unlist(lapply(ans, function(x) x$name))
    fData(object)[names(sc),]$CNNimg.scores <- sc
  }
  object@processingData@processing <- c(processingData(object)@processing,
                                        paste0("Performed CNNimg prediction ", date()))

  unlink(file.path(getwd(), paste0(file_name, "train.csv")))
  unlink(file.path(getwd(), paste0(file_name, "test.csv")))

  if (validObject(object))
    return(object)
}


