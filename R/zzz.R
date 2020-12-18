#load all the data needed for using the package
#especially to avoid error when using runShinyVisualization
.onLoad <- function(...){
  print("All the data needed for this package are loading.")
  library(pRolocdata)
  alldatapRoloc <- c()
  for (i in 1:(length(pRolocdata::pRolocdata()[[3]])/4)){
    alldatapRoloc <- append(alldatapRoloc, pRolocdata::pRolocdata()[[3]][i,][[3]])
  }

  alldatapRoloc <- append(alldatapRoloc, c("FdynCon_n_Msn", "SdynCon_n_Msn", "TdynCon_n_Msn",
                                           "FdynEGF_n_Msn", "SdynEGF_n_Msn", "TdynEGF_n_Msn",
                                           "alldyn", "alldyn_mean", "alldyn_two"))
  alldatapRoloc <<- alldatapRoloc
  d <- c()
  for (i in 1:length(alldatapRoloc)){
    if (exists(alldatapRoloc[i]) == FALSE){
      d <- append(d, alldatapRoloc[i])
    }

  }

  if (!is.null(d)){
    data(list = d)
  }
}
