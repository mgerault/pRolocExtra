#' Classification parameter optimisation for the xgboost algorithm.
#'
#' Note that when performance scores precision, recall and (macro) F1
#' are calculated, any NA values are replaced by 0. This decision is
#' motivated by the fact that any class that would have either a NA
#' precision or recall would result in an NA F1 score and,
#' eventually, a NA macro F1 (i.e. mean(F1)). Replacing NAs by 0s
#' leads to F1 values of 0 and a reduced yet defined final macro F1
#' score.
#'
#' @title xgboost parameter optimisation
#'
#' @param object An instance of class \code{"\linkS4class{MSnSet}"}.
#' @param fcol The feature meta-data containing marker definitions.
#' Default is \code{markers}.
#' @param max_depth The hyper-parameter; the max depth of the tree. Default values are \code{4:7}.
#' @param gamma The hyper-parameter; Minimum loss reduction required to make a further partition on
#'              a leaf node of the tree. The larger gamma is, the more conservative the algorithm will be.
#'              Default values are \code{seq(2.6,3.2, 0.3)}.
#' @param nrounds The max number of boosting iterations.
#' @param times The number of times internal cross-validation is performed. Default is 100.
#' @param test.size The size of test data. Default is 0.2 (20 percent).
#' @param xval The \code{n}-cross validation. Default is 5.
#' @param fun The function used to summarise the \code{xval} macro F1 matrices.
#' @param seed The optional random number generator seed.
#' @param verbose A \code{logical} defining whether a progress bar is displayed.
#' @param ... Additional parameters passed to \code{\link{xgb.train}} from package \code{xgboost}.
#'
#' @return An instance of class \code{"\linkS4class{GenRegRes}"}.
#'
#' @seealso \code{\link{xgboostClassification}} and example therein.
#'
#' @aliases  xgboostOptimization
#'
#' @export

xgboostOptimisation <- function (object, fcol = "markers",
                                 max_depth = 4:7, gamma = seq(2.6,3.2, 0.3),
                                 nrounds = 100,
                                 times = 100, test.size = 0.2, xval = 5,
                                 fun = mean, seed, verbose = TRUE,
                                 ...)
{
  nparams <- 2
  mydata <- pRoloc:::subsetAsDataFrame(object, fcol, train = TRUE)
  if (missing(seed)) {
    seed <- sample(.Machine$integer.max, 1)
  }
  .seed <- as.integer(seed)
  set.seed(.seed)
  .warnings <- NULL
  .f1Matrices <- vector("list", length = times)
  .testPartitions <- .cmMatrices <- vector("list", length = times)
  .results <- matrix(NA, nrow = times, ncol = nparams + 1)
  colnames(.results) <- c("F1", "max_depth", "gamma")
  if (verbose) {
    pb <- txtProgressBar(min = 0, max = xval * times, style = 3)
    ._k <- 0
  }
  #searching for best parameters
  for (.times in 1:times) {
    .size <- ceiling(table(mydata$markers) * test.size)
    .size <- .size[unique(mydata$markers)]
    test.idx <- sampling::strata(mydata, "markers", size = .size,
                       method = "srswor")$ID_unit
    .testPartitions[[.times]] <- test.idx
    .test1 <- mydata[test.idx, ]
    .train1 <- mydata[-test.idx, ]
    xfolds <- caret::createFolds(.train1$markers, xval, returnTrain = TRUE)
    .matrixF1L <- vector("list", length = xval)
    for (.xval in 1:xval) {
      if (verbose) {
        setTxtProgressBar(pb, ._k)
        ._k <- ._k + 1
      }
      .train2 <- .train1[xfolds[[.xval]], ]
      .test2 <- .train1[-xfolds[[.xval]], ]
      .train2_xg <- xgboost::xgb.DMatrix(data = as.matrix(.train2[,-ncol(.train2)]),
                                         label = as.integer(.train2[,ncol(.train2)]) - 1)
      .test2_xg <- xgboost::xgb.DMatrix(data = as.matrix(.test2[,-ncol(.test2)]),
                                        label = as.integer(.test2[,ncol(.test2)]) - 1)
      .matrixF1 <- pRoloc:::makeF1matrix(list(max_depth = max_depth, gamma = gamma))
      for (.max_depth in max_depth) {
        for (.gamma in gamma) {
          params_xg <- list(
            booster = "gbtree",
            eta = 0.001,
            max_depth = .max_depth,
            gamma = .gamma,
            subsample = 0.75,
            colsample_bytree = 1,
            lambda = 1, alpha = 0,
            objective = "multi:softprob",
            eval_metric = "mlogloss",
            num_class = length(levels(.train2$markers))
          )
          model <- xgboost::xgb.train(params = params_xg,
                             data = .train2_xg,
                             nrounds = nrounds,
                             early_stopping_rounds = 10,
                             watchlist = list(train = .train2_xg),
                             verbose = 0,
                             ...)
          ans <- as.data.frame(predict(model, .test2_xg, reshape = TRUE))
          colnames(ans) <- levels(.train2$markers)
          ans <- factor(apply(ans, 1, function(x) colnames(ans)[which.max(x)]), levels = levels(.test2$markers))
          conf <- caret::confusionMatrix(ans, .test2$markers)$table
          .p <- pRoloc:::checkNumbers(MLInterfaces:::.precision(conf))
          .r <- pRoloc:::checkNumbers(MLInterfaces:::.recall(conf))
          .f1 <- MLInterfaces:::.macroF1(.p, .r, naAs0 = TRUE)
          .matrixF1[as.character(.gamma), as.character(.max_depth)] <- .f1
        }
      }
      .matrixF1L[[.xval]] <- .matrixF1
    }
    .summaryF1 <- pRoloc:::summariseMatList(.matrixF1L, fun)
    .f1Matrices[[.times]] <- .summaryF1
    .bestParams <- pRoloc:::getBestParams(.summaryF1)[1:nparams,
                                             1]
    .clcol <- which(names(.train1) == "markers")
    params_xg <- list(
      booster = "gbtree",
      eta = 0.001,
      max_depth = .bestParams["max_depth"],
      gamma = .bestParams["gamma"],
      subsample = 0.75,
      colsample_bytree = 1,
      lambda = 1, alpha = 0,
      objective = "multi:softprob",
      eval_metric = "mlogloss",
      num_class = length(levels(.train2$markers))
    )
    .train1_xg <- xgboost::xgb.DMatrix(data = as.matrix(.train1[,-ncol(.train1)]),
                                       label = as.integer(.train1[,ncol(.train1)]) - 1)
    .test1_xg <- xgboost::xgb.DMatrix(data = as.matrix(.test1[,-ncol(.test1)]),
                                       label = as.integer(.test1[,ncol(.test1)]) - 1)
    ans <- xgboost::xgb.train(params = params_xg,
                       data = .train1_xg,
                       nrounds = nrounds,
                       early_stopping_rounds = 10,
                       watchlist = list(train = .train1_xg),
                       verbose = 0,
                       ...)
    ans <- as.data.frame(predict(ans, .test1_xg, reshape = TRUE))
    colnames(ans) <- levels(.train1$markers)
    ans <- factor(apply(ans, 1, function(x) colnames(ans)[which.max(x)]), levels = levels(.test1$markers))
    .cmMatrices[[.times]] <- conf <- caret::confusionMatrix(ans,
                                                            .test1$markers)$table
    p <- pRoloc:::checkNumbers(MLInterfaces:::.precision(conf), tag = "precision",
                      params = .bestParams)
    r <- pRoloc:::checkNumbers(MLInterfaces:::.recall(conf), tag = "recall",
                      params = .bestParams)
    f1 <- MLInterfaces:::.macroF1(p, r, naAs0 = TRUE)
    .results[.times, ] <- c(f1, .bestParams["max_depth"], .bestParams["gamma"])
  }
  if (verbose) {
    setTxtProgressBar(pb, ._k)
    close(pb)
  }
  # saving results
  .hyperparams <- list(max_depth = max_depth, gamma = gamma)
  .design <- c(xval = xval, test.size = test.size, times = times)
  ans <- new("GenRegRes", algorithm = "xgboost", seed = .seed,
             hyperparameters = .hyperparams, design = .design, results = .results,
             f1Matrices = .f1Matrices, cmMatrices = .cmMatrices,
             testPartitions = .testPartitions, datasize = list(data = dim(mydata),
                                                               data.markers = table(mydata[, "markers"]), train1 = dim(.train1),
                                                               test1 = dim(.test1), train1.markers = table(.train1[,
                                                                                                                   "markers"]), train2 = dim(.train2), test2 = dim(.test2),
                                                               train2.markers = table(.train2[, "markers"])))
  if (!is.null(.warnings)) {
    ans@log <- list(warnings = .warnings)
    sapply(.warnings, warning)
  }
  return(ans)
}

xgboostOptimization <- xgboostOptimisation
