#' Data from Borner and al. 2016, preotein cellular localization
#'
#' A large list containing the data from Borner and al. 2016.
#'  It contains two differents experiment, both using SILAC.
#'  The first is "static", with 6 replicates of 5 fractions.
#'  The second is dynamic, in the sense that the samples are separated with control and treated by EGF.
#'  The control is made of three replicates of 5 fractions, as the EGF treated.
#'
#' @docType data
#' @keywords datasets
#' @name Borner2016
#' @usage data(Borner2016)
#' @format A large list of 12 elements
NULL

#' Synthetic data of dynamic organellar map, from Borner and al. 2016
#'
#' A large \code{\link{MSnSet}} containing the data from Borner and al. 2016.
#'  It contains 2235 proteins, with five fractions, three replicates for two conditions.
#'  The data were normalized (see Borner and al. 2016 for more details)
#'  The best way to access its content is to use the functions :
#'  \code{\link{fData}} and \code{\link{exprs}},
#'  For alldyn : exprs show a tab of 5 columns (the five fractions) and 2235*6 rows (three replicates times two conditions).
#'  For alldyn_two : exprs show a tab of 15 columns (3 replicates times five fractions) and 2235*2 rows (two conditions).
#'  For alldyn_mean : exprs show a tab of 5 columns (the three replicates were averaged) and 2235*2 rows (two conditions).
#'  fData always show a tab of 4 columns containing the gene name, the markers of each protein
#'  (their cellular location), a columns named protcond_name to differentiate proteins between the replicates and the conditions
#'  and a column named cond to differentiate replicates and the conditions.
#'  Here \code{\link{pData}} isn't useful
#'
#' @docType data
#' @keywords datasets
#' @name alldyn
#' @usage alldyn
#' @format A large MSnSet

#' @rdname alldyn
#' @usage alldyn_two
#' @usage alldyn_mean
NULL

#' Data of dynamic organellar map, from Borner and al. 2016
#'
#' A \code{\link{MSnSet}} containing the data from Borner and al. 2016.
#'  It contains 2235 proteins, with five fractions.
#'  F is for first replicate, S the second and T the third.
#'  Con is for control and EGF id for treated by EGF.
#'  n_Msn is because the data were normalized (see Borner and al. 2016 for more details)
#'  The best way to access its content is to use the functions :
#'  \code{\link{fData}} and \code{\link{exprs}},
#'  exprs show a tab of 5 columns (the five fractions) and 2235 rows (three replicates times two conditions)
#'  fData show a tab of 2 columns and 2235 rows containing the gene_name and the markers of each protein
#'  (their cellular location).
#'  Here \code{\link{pData}} isn't useful
#'
#' @docType data
#' @keywords datasets
#' @name FdynCon_n_Msn
#' @usage FdynCon_n_Msn
#' @format A MSnSet

#' @rdname FdynCon_n_Msn
#' @usage SdynCon_n_Msn
#' @usage TdynCon_n_Msn
#' @usage FdynEGF_n_Msn
#' @usage SdynEGF_n_Msn
#' @usage TdynEGF_n_Msn
NULL

