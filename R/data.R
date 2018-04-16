#' Cell composition for 650 individuals from GSE40279. Estimated by running an array method on the original (array) data.
#'
#' @format A data frame with 650 rows and 6 variables
#' @source \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE40279}
"baseline_GSE40279"


#' Simulated amount of methylated reads for 650 individuals from GSE40279. The data was sub-sampled to simulate a 30X coverage.
#' Only 241 sites that are known to differ substaintally between cell types are recorded.
#'
#' @format A data frame with 650 rows and 241 variables
#' @source \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE40279}
"methylation_GSE40279"


#' Simulated amount of total reads reads for 650 individuals from GSE40279. The data was sub-sampled to simulate a 30X coverage.
#' Only 241 sites that are known to differ substaintally between cell types are recorded.
#'
#' @format A data frame with 650 rows and 241 variables
#' @source \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE40279}
"total_reads_GSE40279"


#' A reference containing methylation proportions for pure cell types in the 241 chosen sites.
#'
#' @format A data frame with 241 rows and 7 variables
#' @source \url{https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-016-0943-7}
"reference_blood"


#' Recommended values for a Dirichlet prior on cell composition in blood samples.
#' Estimated by fitting Dirichlet distribution to cell counts data.
#' @format A vector of length 6
"alpha_blood"
