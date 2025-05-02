#' @title Estimates rarefacted size of phenotypic spaces
#'
#' @description \code{rarefact_space_size}
#' @param n Integer vector of length 1 indicating the number of samples to be use for rarefaction (i.e. how many samples per group will be gather at each iteration). Default is the minimum sample size across groups. Integer vector of length 1 indicating the number of samples to be use for rarefaction (i.e. how many samples per group will be gather at each iteration). Default is the minimum sample size across groups.
#' @param replace Logical argument to control if sampling is done with replacement. Default is \code{FALSE}. 
#' @inheritParams template_params
#' @param iterations Integer vector of length 1. Controls how the number of times the rarefaction routine is iterated. Default is 30.
#' @param ... Additional arguments to be passed to \code{\link{space_size}}.
#' @return A data frame containing the mean, minimum, maximum and standard deviation of the size difference across iterations for each pair of groups.
#' @export
#' @name rarefact_space_size
#' @details The function applies a rarefaction sub-sampling procedure for evaluating pairwise space size differences (internally using \code{\link{space_size}}). The size of a phenotypic space might change as a function of number of samples. Hence, ideally, spaces should be compared between groups of similar sample sizes. Rarefaction allows to compare groups of unbalanced sample sizes by randomly re-sampling observations using the same number samples across groups iteratively.  
#' @examples {
#' # load data
#' data("example_space")
#' 
#' # get rarefacted MCP space size 
#' # (try with more iterations on your own data)
#' rarefact_space_size(
#'  formula = group ~ dimension_1 + dimension_2,
#'  data = example_space,
#'  method = "mcp")
#' 
#' # mst rarefacted
#' rarefact_space_size(
#'  formula = group ~ dimension_1 + dimension_2,
#'  data = example_space,
#'  method = "mst")
#' }
#' @seealso \code{\link{rarefact_space_similarity}}, \code{\link{space_size_difference}}
#' @author Marcelo Araya-Salas \email{marcelo.araya@@ucr.ac.cr})
#'
#' @references
#' Araya-Salas, M, & K. Odom. 2022, PhenotypeSpace: an R package to quantify and compare phenotypic trait spaces R package version 0.1.0.

rarefact_space_size <-
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
    if (obs.n < n & !replace) {
      
      if (replace) message("'n' higher than the minimum sample size for at least 1 group, running rarefaction with replacement (replace = TRUE)")
      replace <- TRUE
    } 
  } else 
    n <- obs.n
  
  # run iterations
  space_sizes_list <- pblapply_phtpspc_int(1:iterations, cl = cores, pbar = pb, function(e){
    if (!is.null(seed))
      set.seed(seed + e)
    
    # randomly choose the index of the observations to use for each group levels
    raref_indices <- unlist(lapply(sort(unique(data[, group])), function(x)
      sample(x = rownames(data)[data[, group] == x], size = n, replace = replace)
    ))
    
    overlaps <- space_size(formula = formula, data = data[raref_indices, ], pb = FALSE, cores = 1, ...)
    
    return(overlaps)
  })
  
  space_sizes_mat <- do.call(cbind, lapply(space_sizes_list, '[', -c(1, 2)))
  
  results <- space_sizes_list[[1]][, 1:2]
  results$mean.size <- rowMeans(space_sizes_mat, na.rm = TRUE)
  results$min.size <- apply(X = space_sizes_mat, MARGIN = 1, FUN = min, na.rm = TRUE)  
  results$max.size <- apply(X = space_sizes_mat, MARGIN = 1, FUN = max, na.rm = TRUE)  
  results$sd.size <- apply(X = space_sizes_mat, MARGIN = 1, FUN = stats::sd, na.rm = TRUE)  
  
  return(results)
}
