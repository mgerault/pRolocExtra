# pRolocExtra
pRolocExtra is a package to complete the pRoloc package from L. Gatto and al. (http://bioconductor.org/packages/release/bioc/html/pRoloc.html) 

With this package you will be able to have better visualization of protein cellular location, with interactive plot or not. 
It also include new data from the paper from [Borner et al. 2016](https://elifesciences.org/articles/16950), which include a dynamic component. 
You will be able to compare protein cellular location between two conditions, see the protein movement with vectors, etc.

The package contains also an app to have an interactive visualization of all the data from pRolocdata and the data from Borner and al.

Moreover, pRolocExtra contains three new clustering methods. First is xgboost from the r package of the same name. Second is Constrained Proportionate Assignment (CPA) described by [Jadot et al.](https://pubmed.ncbi.nlm.nih.gov/27923875/) in 2017. Quickly, this method allow to assign probabilities of belonging to each organelle to each protein using a constrained optimization to know the contribution of each mean profile of each organelle to proteins profiles.
Third uses a CNN on the image protein profiles. This function calls a python script that plot the image profile of each protein and then train a CNN on these.
The functions are coded in the same shape as in pRoloc. Xgboost has an optimization step so, as in pRoloc, there are two functions which are called xgboostOptimization and xgboostClassification. The two others are named CPAClassification and CNNimgClassification.


# How to install and use pRolocExtra ?
First, go to Rstudio. Before installing pRolocExtra, you will need to install [pRoloc](http://bioconductor.org/packages/release/bioc/html/pRoloc.html), [pRolocdata](http://bioconductor.org/packages/release/data/experiment/html/pRolocdata.html) and [MSnbase](https://bioconductor.org/packages/release/bioc/html/MSnbase.html) packages from bioconductor in order to use all functionnalities from the app.
Run this commands :

```c
if(!requireNamespace("BiocManager", quietly = TRUE)){
   install.packages("BiocManager")  
}
BiocManager::install(c("pRoloc", "pRolocdata", "MSnbase"))  
```

You can now install pRolocExtra from github : 

```c
if(!requireNamespace("devtools", quietly = TRUE)){
   install.packages("devtools") 
}
devtools::install_github("mgerault/pRolocExtra")
```

You can now load it and run the app with this commands : 

```c
library(pRolocExtra)
runpRolocExtra()
```
