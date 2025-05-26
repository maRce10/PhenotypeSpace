#' @title Calculates rarefacted space overlaps
#'
#' @description \code{rarefact_space_similarity}
#' @inheritParams template_params
#' @param n Integer vector of length 1 indicating the number of samples to be use for rarefaction (i.e. how many samples per group will be gather at each iteration). Default is the minimum sample size across groups.
#' @param replace Logical argument to control if sampling is done with replacement. Default is \code{FALSE}.
#' @param iterations Integer vector of length 1. Controls how the number of times the rarefaction routine is iterated. Default is 30.
#' @param ... Additional arguments to be passed to \code{\link{space_similarity}} for customizing similarity measurements.
#' @return A data frame containing the mean, minimum, maximum and standard deviation of the similarity metric across iterations for each pair of groups. If the similarity metric is not symmetric (e.g. the proportional area of A that overlaps B is not necessarily the same as the area of B that overlaps A, see \code{\link{space_similarity}}) separated columns are supplied for the two comparisons. 
#' @export
#' @name rarefact_space_similarity
#' @details The function applies a rarefaction sub-sampling procedure for evaluating pairwise space similarity (internally using \code{\link{space_similarity}}). The spread and shape of a phenotypic space might change as a function of number of samples. Hence, ideally, spaces should be compared between groups of similar sample sizes. Rarefaction allows to compare groups of unbalanced sample sizes by randomly re-sampling observations using the same number samples across groups iteratively.       
#' @examples {
#' # load data
#' data("example_space")
#' 
#' # get proportion of space that overlaps (try with more iterations on your own data)
#' prop_overlaps <- rarefact_space_similarity(
#'  formula = group ~ dimension_1 + dimension_2,
#'  data = example_space,
#' method = "proportional.overlap", 
#' iterations = 5)
#' 
#' # get minimum convex polygon overlap for each group (non-symmetric)
#' mcp_overlaps <- rarefact_space_similarity(
#'  formula = group ~ dimension_1 + dimension_2,
#'  data = example_space,
#'  iterations = 5)
#' 
#' # convert to non-symmetric triangular matrix
#' rectangular_to_triangular(mcp_overlaps, symmetric = FALSE)
#' }
#' @seealso \code{\link{space_similarity}}, \code{\link{rectangular_to_triangular}}
#' @author Marcelo Araya-Salas \email{marcelo.araya@@ucr.ac.cr})
#'
#' @references
#' Araya-Salas, M, & K. Odom. 2022, PhenotypeSpace: an R package to quantify and compare phenotypic trait spaces R package version 0.1.0.

rarefact_space_similarity <- function(formula,
                                      data,
                                      n = NULL,
                                      replace = FALSE,
                                      seed = NULL,
                                      cores = 1,
                                      pb = TRUE,
                                      iterations = 30,
                                      ...) {
  
  # get term names from formula 
  form_terms <- terms(formula)
  dimensions <- attr(form_terms, "term.labels")
  group <- as.character(form_terms[[2]])
  
  obs.n <- min(table(data[, group]))
  
  if (!is.null(n)) {
    if (obs.n < n) {
      
      if (!replace) 
        message("'n' higher than the minimum sample size for at least 1 group, running rarefaction with replacement (replace = TRUE)")
      replace <- TRUE
    } 
  } else {
    n <- obs.n
    message("Rarefaction sample size = ", n)
  }
  # run iterations
  space_similarities_list <- pblapply_phtpspc_int(1:iterations, cl = cores, pbar = pb, function(e){
    if (!is.null(seed))
      set.seed(seed + e)
    
    raref_indices <- unlist(lapply(sort(unique(data[, group])), function(x)
      sample(x = rownames(data)[data[, group] == x], size = n, replace = replace)
    ))
    
    overlaps <- space_similarity(data[raref_indices, ], formula = formula, pb = FALSE, cores = 1, seed = seed, ...)
    
    return(overlaps)
  })
  
  space_similarities_mat <- do.call(cbind, lapply(space_similarities_list, '[', -c(1, 2)))
  
  results <- space_similarities_list[[1]][, 1:2]
  
  # get type from arguments set in the call
  call.argms <- as.list(base::match.call())[-1]

  if (any(names(call.argms) == "method"))
    method <- call.argms$method else
      method <- "mcp.overlap" # default value in space_similarity()
  
  # create output data frame
  if (!any(names(space_similarities_list[[1]]) %in% c("overlap.1in2", "probability.1in2"))){
  results$mean.similarity <- rowMeans(space_similarities_mat)
  results$min.similarity <- apply(X = space_similarities_mat, MARGIN = 1, FUN = min)  
  results$max.similarity <- apply(X = space_similarities_mat, MARGIN = 1, FUN = max)  
  results$sd.similarity <- apply(X = space_similarities_mat, MARGIN = 1, FUN = stats::sd)
  } else {
    results$mean.similarity.1in2 <- rowMeans(space_similarities_mat[ , seq(1, (iterations * 2) - 1, 2)])
    results$mean.similarity.2in1 <- rowMeans(space_similarities_mat[ , seq(2, (iterations * 2), 2)])
    results$min.similarity.1in2 <- apply(X = space_similarities_mat[ , seq(1, (iterations * 2) - 1, 2)], MARGIN = 1, FUN = min)
    results$min.similarity.2.in1 <- apply(X = space_similarities_mat[ , seq(2, (iterations * 2), 2)], MARGIN = 1, FUN = min)  
    results$max.similarity.1in2 <- apply(X = space_similarities_mat[ , seq(1, (iterations * 2) - 1, 2)], MARGIN = 1, FUN = max)
    results$max.similarity.2.in1 <- apply(X = space_similarities_mat[ , seq(2, (iterations * 2), 2)], MARGIN = 1, FUN = max)  
    results$sd.similarity.1in2 <- apply(X = space_similarities_mat[ , seq(1, (iterations * 2) - 1, 2)], MARGIN = 1, FUN = sd)
    results$sd.similarity.2.in1 <- apply(X = space_similarities_mat[ , seq(2, (iterations * 2), 2)], MARGIN = 1, FUN = sd)  
  }
  
  # rename columns
  if (grepl("distance", method) & method != "probability") {
    names(results) <-gsub("similarity", "distance", names(results))
  } 
  
  if (!grepl("overlap", method) & method != "probability") {
    names(results) <- gsub("similarity", "overlap", names(results))
  }
  
  if (method == "probability") {
    names(results) <-gsub("similarity", "probability", names(results))
  } 
  
  
    return(results)
}
