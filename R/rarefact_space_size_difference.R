#' @title Calculates rarefacted space size differences
#'
#' @description \code{rarefact_space_size_difference}
#' @inheritParams template_params
#' @param n Integer vector of length 1 indicating the number of samples to be use for rarefaction (i.e. how many samples per group will be gather at each iteration). Default is the minimum sample size across groups.
#' @param replace Logical argument to control if sampling is done with replacement. Default is \code{FALSE}.
#' @param iterations Integer vector of length 1. Controls how the number of times the rarefaction routine is iterated. Default is 30.
#' @param ... Additional arguments to be passed to \code{\link{space_size_difference}}.
#' @return A data frame containing the mean, minimum, maximum and standard deviation of the space size difference across iterations for each pair of groups. 
#' @export
#' @name rarefact_space_size_difference
#' @details The function applies a rarefaction sub-sampling procedure for evaluating pairwise space size differences (internally using \code{\link{space_size_difference}}). The size of a phenotypic space might change as a function of number of samples. Hence, ideally, spaces should be compared between groups of similar sample sizes. Rarefaction allows to compare groups of unbalanced sample sizes by randomly re-sampling observations using the same number samples across groups iteratively. 
#' @examples {
#' # load data
#' data("example_space")
#' 
#' # get rarefied size difference using MCP (try with more iterations on your own data)
#' mcp_size_diff <- rarefact_space_size_difference(
#'  formula = group ~ dimension_1 + dimension_2,
#'  data = example_space,
#'  method = "mcp", 
#'  iterations = 5)
#' 
#' # convert to non-symmetric triangular matrix
#' rectangular_to_triangular(mcp_size_diff, symmetric = FALSE)
#' }
#' @seealso \code{\link{rarefact_space_similarity}}, \code{\link{space_size_difference}}
#' @author Marcelo Araya-Salas \email{marcelo.araya@@ucr.ac.cr})
#'
#' @references {
#' Araya-Salas, M, & K. Odom. 2022, PhenotypeSpace: an R package to quantify and compare phenotypic trait spaces R package version 0.1.0.
#' }

rarefact_space_size_difference <-
  function(formula,
           data,
           n = NULL,
           replace = FALSE,
           seed = NULL,
           cores = 1,
           pb = TRUE,
           iterations = 30,
           ...
  ){
    # get term names from formula 
    form_terms <- terms(formula)
    dimensions <- attr(form_terms, "term.labels")
    group <- as.character(form_terms[[2]])
    
  obs.n <- min(table(data[, group]))
  
  if (!is.null(n)) {
    if (obs.n < n) {
      
      if (replace) message("'n' higher than the minimum sample size for at least 1 group, running rarefaction with replacement (replace = TRUE)")
      replace <- TRUE
    } 
  } else 
    n <- obs.n
  
  # run iterations
  space_size_diffs_list <- pblapply_phtpspc_int(1:iterations, cl = cores, pbar = pb, function(e){
    if (!is.null(seed))
      set.seed(seed + e)
    
    raref_indices <- unlist(lapply(sort(unique(data[, group])), function(x)
      sample(x = rownames(data)[data[, group] == x], size = n, replace = replace)
      ))
    
    size_diffs <- space_size_difference(formula = formula, data = data[raref_indices, ], pb = FALSE, output = "rectangular", cores = 1, ...)
    
    return(size_diffs)
  })
  
  space_size_diffs_mat <- do.call(cbind, lapply(space_size_diffs_list, '[', -c(1, 2)))
  
  results <- space_size_diffs_list[[1]][, 1:2]
  results$mean.difference <- rowMeans(space_size_diffs_mat)
  results$min.difference <- apply(X = space_size_diffs_mat, MARGIN = 1, FUN = min)  
  results$max.difference <- apply(X = space_size_diffs_mat, MARGIN = 1, FUN = max)  
  results$sd.difference <- apply(X = space_size_diffs_mat, MARGIN = 1, FUN = stats::sd)  
  
  
  return(results)
}
