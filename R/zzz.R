#load all the data needed for using the package
#especially to avoid error when using runpRolocExtra
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
    if (exists(alldatapRoloc[i]) == FALSE){  # if some data already loaded, don't do it twice
      d <- append(d, alldatapRoloc[i])
    }

  }

  if (!is.null(d)){
    data(list = d)
  }


  inst_conda <- tolower(readline(prompt = "Do you want to configure miniconda now with reticulate package \nin order to use CNNimgClassification ? (yes/no)"
                                 )
                        )
  if(inst_conda == "yes" | inst_conda == "y"){
    x <- tryCatch(reticulate::conda_list(),
                 error = function(e)
                   "no_conda")
    y <- x
    if(class(x) != "character"){
      n_conda <- nrow(x)
      vers_conda <- x$name
      todo <- readline(prompt = paste(n_conda, "versions of conda was found in your computer. \nWhich one do you want to use ? Or do you want to install one with reticulate anyway ? \nChoices :", paste(vers_conda, collapse = ", "), "or install.")
                               )

      if(todo == "install"){
        y <- "inst"
      }
      else if(sum(todo %in% c("default", vers_conda)) == 0){
        print("Invalid choice")
      }
      else{
        reticulate::use_condaenv(condaenv = todo)
      }
    }
    if(class(y) == "character"){
      message("No conda installation was found in your computer. \nreticulate will install it; the environment name will be 'r-reticulate' \nAlso, all the python librairies needed will be installed.")

      reticulate::install_miniconda()

      reticulate::use_condaenv(condaenv = "r-reticulate")

      message("Installing librairies...")
      reticulate::py_install(packages = c("matplotlib",
                                          "pandas",
                                          "opencv")
                             )

      reticulate::py_install(c("pytorch", "torchvision"),
                             channel = c("pytorch")
                             )  ## torchvision is only accessible in channel pytorch
    }
  }
  else if(inst_conda == "no" | inst_conda == "n"){
    message("If you want to use CNNimgClassification, be sure that conda is in your computer.")
  }
  else{
    print("Invalid choice")
  }
}
