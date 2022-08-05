#' Classification using constrained proportionate assignment method
#'
#' @title CPA classification
#'
#' @param object An instance of class \code{"\linkS4class{MSnSet}"}.
#' @param scores One of \code{"prediction"}, \code{"all"} or
#'     \code{"none"} to report the score for the predicted class
#'     only, for all classes or none.
#' @param fcol The feature meta-data containing marker definitions.
#'     Default is \code{markers}.
#' @param method An integer (1, 2 or 3) specifying wich Barzilai-Borwein steplength to use; default is 3.
#'               For more details see documentation from \code{\link{spg}} function from package \code{BB}.
#' @param ... Additional parameters passed to \code{\link{spg}} from
#'     package \code{BB}.
#'
#' @details CPA was originally described by Michel Jadot et al. in 'Accounting for protein subcellular localization:
#'          A compartmental map of the rat liver proteome.', 2017. The goal of this method is to assign
#'          probabilities to each protein of belonging to an organelle. This way, we can potentially assign multiple
#'          organelles to proteins (a protain may have multiple locations).
#'          For this, it uses the mean profile of each organelle of the train dataset and then find the best coefficients
#'          to each profile  to obtain the profile from a given protein. The coefficients are constrained to sum to one and
#'          bounded between 0 and 1 (so we can interpret them as probabilities). To solve this constrained optimization
#'          problem, we use \code{spg} function from \code{BB} package, as in the article.
#'          However, in this function, any treatments on the mean profiles are done, unlike in the article. So be sure
#'          of your reference proteins and use normalized data.
#'
#' @return An instance of class \code{"\linkS4class{MSnSet}"} with
#'     \code{CPA} and \code{CPA.scores} feature variables storing the
#'     classification results and scores respectively.
#'
#' @export
#'
#' @examples
#' library(pRolocExtra)
#' data(tan2009r1)
#' res <- CPAClassification(tan2009r1)
#' getPredictions(res, fcol = "CPA")
#' getPredictions(res, fcol = "CPA", t = 0.75)
#' plot2D(res, fcol = "CPA")

CPAClassification <- function (object, scores = c("prediction", "all", "none"),
                               fcol = "markers", method = 3, ...)
{
  scores <- match.arg(scores)

  # set train and test set
  trainInd <- which(fData(object)[, fcol] != "unknown")
  testInd <- which(fData(object)[, fcol] == "unknown")
  trainSet <- pRoloc:::subsetAsDataFrame(object, fcol, train = TRUE)
  testSet <- pRoloc:::subsetAsDataFrame(object, fcol, train = FALSE)
  .trainSep <- pRoloc:::separateDataSet(trainSet, fcol)
  .testSep <- pRoloc:::separateDataSet(testSet, fcol)

  # get men profiles
  mean_prof = trainSet %>%
    dplyr::group_by_at(fcol) %>%
    dplyr::summarise(across(where(is.numeric), mean))

  mrk = mean_prof[[fcol]]
  m = length(mrk)
  mean_prof[[fcol]] <- NULL

  obj.fn <- function(x, y) sum((y - x %*% as.matrix(mean_prof))**2)
  # gradient from obj fun
  grobj.fn <- function(x, y){
    m <- as.matrix(mean_prof)
    gr <- c()
    for(i in 1:length(x)){
      r <- c()
      for(f in 1:ncol(mean_prof)){
        r1 <- m[i,f] * (y[f] - x %*% m[,f])
        r <- append(r, as.numeric(r1))
      }
      r <- -2*sum(r)
      gr <- append(gr, r)
    }
    return(gr)
  }

  #bounds
  lb <- rep(0, m)
  ub <- rep(1, m)

  # constraint
  Amat = matrix(rep(1, m), nrow = 1, ncol = m, byrow = TRUE)
  b <- 1
  meq <- 1 # first condition is an equality

  #initials
  x0 = runif(m)

  ans <- apply(.testSep$theData, 1, function(x){
    res <- BB::spg(par = x0,
                  fn = obj.fn,
                  gr = grobj.fn,
                  y = x,
                  project = "projectLinear",
                  projectArgs = list(A = Amat, b = b, meq = meq),
                  lower = lb, upper = ub,
                  method = method,                        # may meth 1
                  control = list("trace" = FALSE), ...);
    res$par
  }, simplify = FALSE)

  ans <- t(as.data.frame(ans))
  colnames(ans) <- levels(mrk)
  ans <- as.data.frame(ans)

  predictedLabels <- factor(apply(ans, 1, function(x) colnames(ans)[which.max(x)]), levels = levels(mrk))
  temp <- rep("", length(trainInd) + length(testInd))
  i <- 1:length(trainInd)
  temp[trainInd[i]] <- as.character(.trainSep$theLabels[i])
  i <- 1:length(testInd)
  temp[testInd[i]] <- as.character(predictedLabels[i])
  fData(object)$CPA <- temp

  # get back probabilities from CPA
  if (scores == "all") {
    nbLabels <- length(levels(mrk))
    tempScores <- matrix(rep(0, nbLabels * (length(trainInd) +
                                              length(testInd))), ncol = nbLabels)
    colnames(tempScores) <- levels(.trainSep$theLabels)
    for(i in 1:length(trainInd)){
      tempScores[trainInd[i], .trainSep$theLabels[i]] <- 1
    }
    ans <- as.matrix(ans)
    i <- 1:length(testInd)
    tempScores[testInd[i], ] <- ans[i, ]
    scoreMat <- tempScores
    colnames(scoreMat) <- paste0(colnames(scoreMat), ".CPA.scores")
    fData(object)$CPA.all.scores <- scoreMat
  }
  else if (scores == "prediction") {
    nbLabels <- length(levels(mrk))
    tempScores <- rep(0, length(trainInd) + length(testInd))
    i <- 1:length(trainInd)
    tempScores[trainInd[i]] <- 1
    i <- 1:length(testInd)
    tempScores[testInd[i]] <- apply(ans, 1, max)
    fData(object)$CPA.scores <- tempScores
  }
  object@processingData@processing <- c(processingData(object)@processing,
                                        paste0("Performed CPA prediction", date()))
  if (validObject(object))
    return(object)
}

