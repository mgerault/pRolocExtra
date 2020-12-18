library(rio)
library(tidyverse)
library(pRoloc)
library(matlab)


Borner2016 <- import_list("data-raw/data_SILAC_dynamic.xlsx")
Dyn_data <- list("Con" = as.data.frame(Borner2016$`SILAC Dynamic Con`), "EGF" = as.data.frame(Borner2016$`SILAC Dynamic +EGF`))

#Take off proteins which are duplicated in order to use rownames
Dupl_Con <- which(duplicated(Dyn_data$Con$`Lead gene name`) == TRUE)
Dupl_EGF <- which(duplicated(Dyn_data$EGF$`Lead gene name`) == TRUE)
Dyn_data$Con <- Dyn_data$Con[-Dupl_Con,]
Dyn_data$EGF <- Dyn_data$EGF[-Dupl_EGF,]

rownames(Dyn_data$Con) <- Dyn_data$Con$`Lead gene name`
rownames(Dyn_data$EGF) <- Dyn_data$EGF$`Lead gene name`

#Find the common proteins between the two experiment conditions
common_prot <- intersect(Dyn_data$Con$`Lead gene name`, Dyn_data$EGF$`Lead gene name`)

Dyn_data$Con <- Dyn_data$Con[common_prot,]
Dyn_data$EGF <- Dyn_data$EGF[common_prot,]



#Three replicates Con and three replicates EGF -> MSnSet object
FdynCon <- readMSnSet2(Dyn_data$Con,
                       ecol = c(7:11), fnames = "Lead gene name")
SdynCon <- readMSnSet2(Dyn_data$Con,
                       ecol = c(12:16), fnames = "Lead gene name")
TdynCon <- readMSnSet2(Dyn_data$Con,
                       ecol = c(17:21), fnames = "Lead gene name")
FdynEGF <- readMSnSet2(Dyn_data$EGF,
                       ecol = c(7:11), fnames = "Lead gene name")
SdynEGF <- readMSnSet2(Dyn_data$EGF,
                       ecol = c(12:16), fnames = "Lead gene name")
TdynEGF <- readMSnSet2(Dyn_data$EGF,
                       ecol = c(17:21), fnames = "Lead gene name")

#change NA values in unknown in the organellar column
FdynCon = nainunknown(FdynCon, fcol = "Organellar_markers")
SdynCon = nainunknown(SdynCon, fcol = "Organellar_markers")
TdynCon = nainunknown(TdynCon, fcol = "Organellar_markers")
FdynEGF = nainunknown(FdynEGF, fcol = "Organellar_markers")
SdynEGF = nainunknown(SdynEGF, fcol = "Organellar_markers")
TdynEGF = nainunknown(TdynEGF, fcol = "Organellar_markers")

#take off protein with NA or infinite values
FdynCon = nafilterdata(FdynCon)
SdynCon = nafilterdata(SdynCon)
TdynCon = nafilterdata(TdynCon)
FdynEGF = nafilterdata(FdynEGF)
SdynEGF = nafilterdata(SdynEGF)
TdynEGF = nafilterdata(TdynEGF)

#keep only the common proteins between the six different data
all_common_prot <- Reduce(intersect, list(featureNames(FdynCon), featureNames(SdynCon), featureNames(TdynCon),
                                          featureNames(FdynEGF), featureNames(SdynEGF), featureNames(TdynEGF)))



#normalization of the as Borner and al. 2016
#The ratios are inverted in order to have light on heavy (the data in the excel are log2 transformed)
FdynCon_n <- 1/2^exprs(FdynCon)[all_common_prot,]
SdynCon_n <- 1/2^exprs(SdynCon)[all_common_prot,]
TdynCon_n <- 1/2^exprs(TdynCon)[all_common_prot,]
FdynEGF_n <- 1/2^exprs(FdynEGF)[all_common_prot,]
SdynEGF_n <- 1/2^exprs(SdynEGF)[all_common_prot,]
TdynEGF_n <- 1/2^exprs(TdynEGF)[all_common_prot,]

#each fractions is multiplied by its average yields
#(find approximately on a figure of the supplementary file of Borner and al. 2016)
#1 : 27, 2 : 27.8, 3 : 10, 4 : 24, 5 : 11.2
weight_yield <- c(0.27,0.278,0.1,0.24,0.112)
for (i in 1:5){
  FdynCon_n[,i] <- FdynCon_n[,i]*weight_yield[i]
  SdynCon_n[,i] <- SdynCon_n[,i]*weight_yield[i]
  TdynCon_n[,i] <- TdynCon_n[,i]*weight_yield[i]
  FdynEGF_n[,i] <- FdynEGF_n[,i]*weight_yield[i]
  SdynEGF_n[,i] <- SdynEGF_n[,i]*weight_yield[i]
  TdynEGF_n[,i] <- TdynEGF_n[,i]*weight_yield[i]
}

#then, each ratio L/H for each protein is divided by the total sum of the five ratio L/H
Fs <- c()
for (i in 1:nrow(FdynCon_n)){
  Fs <- append(Fs, sum(FdynCon_n[i,]))
}
FdynCon_n <- FdynCon_n/Fs

Ss <- c()
for (i in 1:nrow(FdynCon_n)){
  Ss <- append(Ss, sum(SdynCon_n[i,]))
}
SdynCon_n <- SdynCon_n/Ss

Ts <- c()
for (i in 1:nrow(FdynCon_n)){
  Ts <- append(Ts, sum(TdynCon_n[i,]))
}
TdynCon_n <- TdynCon_n/Ts

EFs <- c()
for (i in 1:nrow(FdynCon_n)){
  EFs <- append(EFs, sum(FdynEGF_n[i,]))
}
FdynEGF_n <- FdynEGF_n/EFs

ESs <- c()
for (i in 1:nrow(FdynCon_n)){
  ESs <- append(ESs, sum(SdynEGF_n[i,]))
}
SdynEGF_n <- SdynEGF_n/ESs

ETs <- c()
for (i in 1:nrow(FdynCon_n)){
  ETs <- append(ETs, sum(TdynEGF_n[i,]))
}
TdynEGF_n <- TdynEGF_n/ETs

#the data are now normalized

#adding the gene_name and the markers to the data frame
FdynCon_n <- as.data.frame(FdynCon_n)
FdynCon_n$gene_name <- rownames(FdynCon_n)
FdynCon_n$markers <- as.factor(fData(FdynCon)[all_common_prot,]$Organellar_markers)

SdynCon_n <- as.data.frame(SdynCon_n)
SdynCon_n$gene_name <- rownames(SdynCon_n)
SdynCon_n$markers <- fData(SdynCon)[all_common_prot,]$Organellar_markers

TdynCon_n <- as.data.frame(TdynCon_n)
TdynCon_n$gene_name <- rownames(TdynCon_n)
TdynCon_n$markers <- fData(TdynCon)[all_common_prot,]$Organellar_markers

FdynEGF_n <- as.data.frame(FdynEGF_n)
FdynEGF_n$gene_name <- rownames(FdynEGF_n)
FdynEGF_n$markers <- fData(FdynEGF)[all_common_prot,]$Organellar_markers

SdynEGF_n <- as.data.frame(SdynEGF_n)
SdynEGF_n$gene_name <- rownames(SdynEGF_n)
SdynEGF_n$markers <- fData(SdynEGF)[all_common_prot,]$Organellar_markers

TdynEGF_n <- as.data.frame(TdynEGF_n)
TdynEGF_n$gene_name <- rownames(TdynEGF_n)
TdynEGF_n$markers <- fData(TdynEGF)[all_common_prot,]$Organellar_markers


colnames(FdynCon_n)[1:5] <- c("F1", "F2", "F3", "F4", "F5")
colnames(SdynCon_n)[1:5] <- c("F1", "F2", "F3", "F4", "F5")
colnames(TdynCon_n)[1:5] <- c("F1", "F2", "F3", "F4", "F5")
colnames(FdynEGF_n)[1:5] <- c("F1", "F2", "F3", "F4", "F5")
colnames(SdynEGF_n)[1:5] <- c("F1", "F2", "F3", "F4", "F5")
colnames(TdynEGF_n)[1:5] <- c("F1", "F2", "F3", "F4", "F5")

#convert in MSnSet object
FdynCon_n_Msn <- readMSnSet2(FdynCon_n,ecol = c(1:5), fnames = "gene_name")
SdynCon_n_Msn <- readMSnSet2(SdynCon_n,ecol = c(1:5), fnames = "gene_name")
TdynCon_n_Msn <- readMSnSet2(TdynCon_n,ecol = c(1:5), fnames = "gene_name")
FdynEGF_n_Msn <- readMSnSet2(FdynEGF_n,ecol = c(1:5), fnames = "gene_name")
SdynEGF_n_Msn <- readMSnSet2(SdynEGF_n,ecol = c(1:5), fnames = "gene_name")
TdynEGF_n_Msn <- readMSnSet2(TdynEGF_n,ecol = c(1:5), fnames = "gene_name")



#new form for further analysis : the three replicates times the two conditions put in one tab (5 columns, 2235*6 rows)
alldyn <- rbind(FdynCon_n, SdynCon_n, TdynCon_n,
                FdynEGF_n, SdynEGF_n, TdynEGF_n)

#differentiate proteins betwween replicates and conditions
rownames(alldyn)[1:nrow(FdynCon_n)] <- paste0(FdynCon_n$gene_name, "_C1")
rownames(alldyn)[(nrow(FdynCon_n)+1):(2*nrow(FdynCon_n))] <- paste0(FdynCon_n$gene_name, "_C2")
rownames(alldyn)[(2*nrow(FdynCon_n)+1):(3*nrow(FdynCon_n))] <- paste0(FdynCon_n$gene_name, "_C3")
rownames(alldyn)[(3*nrow(FdynCon_n)+1):(4*nrow(FdynCon_n))] <- paste0(FdynCon_n$gene_name, "_E1")
rownames(alldyn)[(4*nrow(FdynCon_n)+1):(5*nrow(FdynCon_n))] <- paste0(FdynCon_n$gene_name, "_E2")
rownames(alldyn)[(5*nrow(FdynCon_n)+1):(6*nrow(FdynCon_n))] <- paste0(FdynCon_n$gene_name, "_E3")

alldyn$protcond_name <- rownames(alldyn)

#convert in MSnSet object, create a condition column
alldyn <- readMSnSet2(alldyn,ecol = c(1:5), fnames = "protcond_name")
fData(alldyn)$cond <- c(linspace(1,1,2235), linspace(2,2,2235), linspace(3,3,2235),
                        linspace(4,4,2235), linspace(5,5,2235), linspace(6,6,2235))

#same process but with another form : three replicates control + three replicates EGF
#a tab with 15 columns and 2235*2 rows
alldyn_two <- rbind(cbind(FdynCon_n, SdynCon_n, TdynCon_n),
                    cbind(FdynEGF_n, SdynEGF_n, TdynEGF_n))

colnames(alldyn_two) <- c("F11", "F12", "F13", "F14", "F15", "gene_name1", "markers1",
                          "F21", "F22", "F23", "F24", "F25", "gene_name2", "markers2",
                          "F31", "F32", "F33", "F34", "F35", "gene_name", "markers")

#same markers and gene name between the replicates  --> take off the duplicates
alldyn_two <- alldyn_two[,-c(6,7,13,14)]

#rename the proteins to differentiate them between the two conditions
rownames(alldyn_two)[1:nrow(FdynCon_n)] <- paste0(FdynCon_n$gene_name, "_C")
rownames(alldyn_two)[(nrow(FdynCon_n)+1):(2*nrow(FdynCon_n))] <- paste0(FdynCon_n$gene_name, "_E")

alldyn_two$protcond_name <- rownames(alldyn_two)


alldyn_two<- readMSnSet2(alldyn_two,ecol = c(1:15), fnames = "protcond_name")
fData(alldyn_two)$cond <- c(linspace(1,1,2235), linspace(2,2,2235))



#same process but averaging the replicates : a tab with 5 columns and 2235*2 rows
dynCon_mean <- FdynCon_n
dynEGF_mean <- FdynEGF_n
for (i in 1:5){
  dynCon_mean[[i]] <- rowMeans(cbind(FdynCon_n[[i]], SdynCon_n[[i]], TdynCon_n[[i]]))
  dynEGF_mean[[i]] <- rowMeans(cbind(FdynEGF_n[[i]], SdynEGF_n[[i]], TdynEGF_n[[i]]))
}


alldyn_mean <- rbind(dynCon_mean, dynEGF_mean)


rownames(alldyn_mean)[1:nrow(FdynCon_n)] <- paste0(FdynCon_n$gene_name, "_C")
rownames(alldyn_mean)[(nrow(FdynCon_n)+1):(2*nrow(FdynCon_n))] <- paste0(FdynCon_n$gene_name, "_E")

alldyn_mean$protcond_name <- rownames(alldyn_mean)

alldyn_mean <- readMSnSet2(alldyn_mean,ecol = c(1:5), fnames = "protcond_name")
fData(alldyn_mean)$cond <- c(linspace(1,1,2235), linspace(2,2,2235))


usethis::use_data(alldyn, overwrite = TRUE)
usethis::use_data(alldyn_mean, overwrite = TRUE)
usethis::use_data(alldyn_two, overwrite = TRUE)
usethis::use_data(Borner2016, overwrite = TRUE)

usethis::use_data(FdynCon_n_Msn, overwrite = TRUE)
usethis::use_data(SdynCon_n_Msn, overwrite = TRUE)
usethis::use_data(TdynCon_n_Msn, overwrite = TRUE)
usethis::use_data(FdynEGF_n_Msn, overwrite = TRUE)
usethis::use_data(SdynEGF_n_Msn, overwrite = TRUE)
usethis::use_data(TdynEGF_n_Msn, overwrite = TRUE)

