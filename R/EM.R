E_1_step <- function(P, Pi) {
    N <- nrow(P)
    J <- nrow(Pi)
    K <- ncol(P)

    E <- array(dim = c(N, J, K))

    constants <- tcrossprod(P, Pi)  #  P %*% t(Pi)
    for (k in 1:K) {
        E[, , k] <- outer(P[, k], Pi[, k])/constants
    }

    return(E)
}

E_0_step <- function(P, Pi) {
    N <- nrow(P)
    J <- nrow(Pi)
    K <- ncol(P)

    E <- array(dim = c(N, J, K))

    constants <- tcrossprod(P, 1 - Pi)  # P %*% (1-t(Pi))
    for (k in 1:K) {
        E[, , k] <- outer(P[, k], 1 - Pi[, k])/constants
    }

    return(E)
}

Maximize_P <- function(E_0, E_1, R, Y, reads_i, alpha) {
    K <- dim(E_0)[3]
    constants <- reads_i + sum(alpha) - K

    N <- nrow(R)

    result <- matrix(nrow = N, ncol = K)

    for (k in 1:K) {
        nominator <- rowSums(E_1[, , k] * Y) + rowSums(E_0[, , k] * (R - Y)) + alpha[k] - 1
        result[, k] <- nominator/constants
    }

    return(result)
}

Maximize_Pi <- function(E_0, E_1, Y, R) {
    J <- ncol(Y)
    K <- dim(E_0)[3]
    result <- matrix(nrow = J, ncol = K)
    for (k in 1:K) {
        nominator <- colSums(Y * E_1[, , k])
        denominator <- colSums(Y * E_1[, , k]) + colSums((R - Y) * E_0[, , k])

        result[, k] <- nominator/denominator
    }

    return(result)
}


## R - a matrix of individuals (rows) on sites (columns), containing the total number of reads for each site, in each individual.
## Y - a matrix of individuals (rows) on sites (columns), containing the number of methylated reads for each site, in each individual.
## Pi - A matrix of sites (rows) on cell types (columns), containing the probability for methylation in each site, in each cell type.  If you need a random Pi matrix use
## 'initialize.Pi()'.  P - A matrix of individuals (rows) on cell types (columns), consisted of the estimated proportion of each cell type in each
## individual.  The top n_known_samples row should be the cell proportions of the known individuals (if any), and the rest should be initialized with
## initial.P().  n_known_samples: an integer, the number of known individuals.  If not known, initialize to a random value between 0 and 1.  alpha - a
## list of length K (amount of cell types), containing the hyper-paramters of the dirichlet prior. If None, a non-informative prior is used.
run_EM <- function(R, Y, Pi, P, alpha, n_known_samples = 0, iterations = 200, estimate_Pi = F, estimate_P = T) {

    minimum_cell_proportion <- 0.001

    # The samples to estimate + the known samples.
    reads_i <- rowSums(R)
    progress_bar <- txtProgressBar(min = 0, max = iterations)
    for (i in 1:iterations) {
        E_0 <- E_0_step(P, Pi)
        E_1 <- E_1_step(P, Pi)

        if (estimate_P) {
            maximize_P_results <- Maximize_P(E_0, E_1, R, Y, reads_i, alpha = alpha)
            if (n_known_samples > 0) {
                P[-(1:n_known_samples), ] <- maximize_P_results[-(1:n_known_samples), ]
            } else {
                P <- maximize_P_results
            }

            ## make sure that P doesn't leave allowed range to prevent numeric errors.
            P[P < 0] <- minimum_cell_proportion
            P[P > 1] <- 1 - minimum_cell_proportion
        }

        if (estimate_Pi) {
            Pi <- Maximize_Pi(E_0, E_1, Y, R)
        }

        setTxtProgressBar(progress_bar, i)
    }

    return(list(P = P, Pi = Pi))
}

initialize_P <- function(n_individuals, n_cell_types, alpha) {
    initial_P <- rep(alpha/sum(alpha), n_individuals)
    initial_P <- matrix(initial_P, nrow = n_individuals, byrow = T)

    return(initial_P)
}

#' @importFrom stats runif
initialize_Pi <- function(n_sites, n_cell_types) {
    Pi <- runif(n_sites * n_cell_types)
    Pi <- matrix(Pi, nrow = n_sites)

    return(Pi)
}


### Exported Functions ---------------------------------------------------
#' Add together two numbers.
#'
#' @param methylation a matrix of individuals (rows) on sites (columns), containing the number of methylated reads for each site, in each individual.
#' @param total_reads a matrix of individuals (rows) on sites (columns), containing the total number of reads for each site, in each individual.
#' @param reference a matrix of sites (rows) on cell types (columns), containing the probability for methylation in each site, in each cell type.
#' @param alpha a vector containing the hyper-parameters for the dirichelt prior. One value for each cell type. If NA, it is initiallized to 1/(number of cell types).
#' @param iterations the number of iterations to use in the EM algorithm.
#' @return A matrix of individuals (rows) on cell types (columns) containing the estimated proportion of each cell type, in each individual.
#' @examples
#' ## Prepare the methylation and total reads matrices
#' methylation <- as.matrix(methylation_GSE40279)
#' total_reads <- as.matrix(total_reads_GSE40279)
#' ## Remove the IDs column from the reference
#' Pi <- as.matrix(reference_blood[,-1])
#'
#' ## Run Bisect. You should use around 200 iterations. I choose than to accelarate the example.
#' results <- bisect_supervised(methylation, total_reads, Pi, alpha_blood, iterations = 10)
#' @importFrom utils txtProgressBar
#' @importFrom utils setTxtProgressBar
#' @export
bisect_supervised <- function(methylation, total_reads, reference, alpha = NA, iterations = 200) {
    N <- nrow(methylation)
    J <- ncol(methylation)
    K <- ncol(reference)

    if (length(alpha) == 1) {
        if (is.na(alpha)) {
            alpha <- rep(1, K)
        }
    }

    initial_P <- initialize_P(N, J, alpha)

    results <- run_EM(total_reads, methylation, reference, initial_P, alpha, iterations = iterations, estimate_P = T)

    return(results$P)
}

#' Add together two numbers.
#'
#' @param methylation_unkown_samples a matrix of individuals (rows) on sites (columns), containing the number of methylated reads for each site, in each individual for the samples with unknown cell composition.
#' @param total_reads_unknown_samples a matrix of individuals (rows) on sites (columns), containing the total number of reads for each site, in each individual for the samples with unknown cell composition.
#' @param methylation_known_samples a matrix of individuals (rows) on sites (columns), containing the number of methylated reads for each site, in each individual for the samples with known cell composition.
#' @param total_reads_known_samples a matrix of individuals (rows) on sites (columns), containing the total number of reads for each site, in each individual for the samples with known cell composition.
#' @param cell_composition_known_samples a matrix of individuals (rows) on cell types (columns), containing the proportion of each cell type, in each known sample.
#' @param alpha a vector containing the hyper-parameters for the dirichelt prior. One value for each cell type. If NA, it is initiallized to 1/(number of cell types).
#' @param iterations the number of iterations to use in the EM algorithm.
#' @return A list containing P, a matrix of estimated cell proportions for the unknown samples, and Pi, an estimated reference (the probability of methylation in each cell type).
#' @examples
#' ## Randomly choose samples to be used as known
#' n_known_samples <- 50
#' known_samples_indices <- sample.int(nrow(baseline_GSE40279), size = n_known_samples)
#' known_samples <- as.matrix(baseline_GSE40279[known_samples_indices, ])
#'
#' ## Fit a dirichlet distribution to known samples to use as prior
#' fit_dirichlet <- sirt::dirichlet.mle(as.matrix(known_samples))
#' alpha <- fit_dirichlet$alpha
#'
#' ## Prepare the 4 needed matrices
#' methylation_known <- methylation_GSE40279[known_samples_indices, ]
#' methylation_unknown <-methylation_GSE40279[-known_samples_indices, ]
#' total_known <- total_reads_GSE40279[known_samples_indices, ]
#' total_unknown <- total_reads_GSE40279[-known_samples_indices, ]
#'
#' ## Run Bisect. You should use around 200 iterations. I choose than to accelarate the example.
#' results <- bisect_semi_suprevised(methylation_unknown, total_unknown,
#'                                   methylation_known, total_known,
#'                                   known_samples, alpha, iterations = 10)
#' @export
bisect_semi_suprevised <- function(methylation_unkown_samples, total_reads_unknown_samples, methylation_known_samples, total_reads_known_samples, cell_composition_known_samples,
    alpha = NA, iterations = 200) {
    N_unknown <- nrow(methylation_unkown_samples)
    J <- ncol(methylation_unkown_samples)
    K <- ncol(cell_composition_known_samples)
    N_known <- nrow(methylation_known_samples)

    if (length(alpha) == 1) {
        if (is.na(alpha)) {
            alpha <- rep(1, K)
        }
    }

    ## Estimate a starting point for Pi using only known samples
    print("estimating reference .........")

    initial_Pi <- initialize_Pi(J, K)
    results_known_samples <- run_EM(R = total_reads_known_samples, Y = methylation_known_samples, Pi = initial_Pi, P = cell_composition_known_samples,
        alpha = alpha, n_known_samples = N_known, iterations = iterations, estimate_Pi = T, estimate_P = F)

    initial_Pi <- results_known_samples$Pi

    ## Estimate P and Pi together, use the inital Pi from before as a starting point
    unknown_P <- initialize_P(N_unknown, J, alpha)
    initial_P <- rbind(cell_composition_known_samples, unknown_P)

    Y <- rbind(methylation_known_samples, methylation_unkown_samples)
    R <- rbind(total_reads_known_samples, total_reads_unknown_samples)

    n_known_samples <- nrow(methylation_known_samples)

    print("estimating cell composition .........")
    results <- run_EM(R, Y, initial_Pi, initial_P, alpha, n_known_samples = n_known_samples, estimate_Pi = T, iterations = iterations)

    ## Only return P for unknown samples
    results$P <- results$P[-(1:N_known), ]
    return(results)
}
