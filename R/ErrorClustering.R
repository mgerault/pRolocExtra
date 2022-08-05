#' @title Get the approximate error of a clustering algorithm
#'
#' @description
#' \code{ErrorClustering} perform n loop to calculate the loss of the clustering algorithm from the pRoloc package.
#' Take only the protein from known locations, separate randomly in 80/20 train set / test set.
#' Then count the number of mistake made by the algorithm.
#'
#' @param data A \code{\link{MSnSet}} object
#' @param n An integer which correspond to the number of loop
#' @param cmet A character which correspond to the clusterig method of the pRoloc package
#' knn, svm, naiveBayes, perTurbo, nnet (neural network) and rf (random forest).
#' @param train_size A numeric between 0 and 1 corresponding to train size in percentage
#'
#' @return a list containing the loss and the clustering score of the proteins which were miss located at each loop;
#' The mean loss and the mean clustering score of the proteins which were miss located
#'
#' @seealso  \code{\link{svmOptimisation}} from pRoloc package and  \code{\link{datavisupca}} for more details
#'
#' @export
#'
#' @examples
#'
#' library(pRolocExtra)
#' err_func(tan2009r1, 5)

ErrorClustering <- function(data, n, cmet = "knn", train_size = 0.8){
  all_err <- c()
  all_score_err <- c()

  unk_id <- which(MSnbase::fData(data)$markers == "unknown")  #get only the protein which we know their location : the original training set from the data
  .data <- data[-unk_id,]

  for (i in 1:n){ #the number of loop
    train_id <- sample(1:nrow(MSnbase::exprs(.data)), train_size*nrow(MSnbase::exprs(.data))) #building the training and the test set
    test_id <- setdiff(1:nrow(MSnbase::exprs(.data)), train_id)

    .data_test <- MSnbase::fData(.data)[test_id,]
    MSnbase::fData(.data)$markers[test_id] <- "unknown"  #change 20% of the marker training set, will be used for testing

    .data_clus <- datavisupca(.data, mfcol = "markers", method = cmet, sh.gr = FALSE)  #cluster the data with the method of our choice thanks datavisupca function

    err_id <- which(as.character(.data_test$markers) != as.character(MSnbase::fData(.data_clus)[test_id,][[cmet]])) #find the mistakes made by the algorithm

    err <- length(err_id)/length(test_id)                                             #the number of mistakes divided by the number of proteins contained by the test set
    score_err <- mean(MSnbase::fData(.data_clus[test_id,])[[paste0(cmet, ".scores")]][err_id]) #the mean of clustering score from the proteins that were incorrectly assigned

    all_err <- append(all_err, err)                   #save the values in a list
    all_score_err <- append(all_score_err, score_err)
  }

  return(list("Clus_loss" = all_err, "Clus_score" = all_score_err,               #return the list created and the mean of this list
              "Mean_loss" = mean(all_err), "Mean_score" = mean(all_score_err)))
  #a problem persist : perTurbo algorithm lead to more than 90% of error
}
