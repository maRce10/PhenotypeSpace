#' @title Estimates rarefacted size of phenotypic spaces
#'
#' @description \code{rarefact_space_size}
#' @usage rarefact_space_size(X, dimensions, group,  n = NULL, replace = FALSE, 
#' seed = NULL, parallel = 1, pb = TRUE, iterations = 30,...)
#' @param X Data frame containing columns for the dimensions of the phenotypic space (numeric) and a categorical or factor column with group labels. 
#' @param dimensions Character vector with the names of columns containing the dimensions of the phenotypic space.
#' @param group Character vector with the name of the column (character or factor) containing group labels.
#' @param n Integer vector of length 1 indicating the number of samples to be use for rarefaction (i.e. how many samples per group will be gather at each iteration). Default is the minimum sample size across groups. Integer vector of length 1 indicating the number of samples to be use for rarefaction (i.e. how many samples per group will be gather at each iteration). Default is the minimum sample size across groups.
#' @param replace Logical argument to control if sampling is done with replacement. Default is \code{FALSE}. 
#' @param seed Integer vector of length 1 setting the seed (see \code{\link[base]{set.seed}}). If used results should be the same on different runs, so it makes them replicable.
#' @param parallel Integer vector of length 1. Controls whether parallel computing is applied. It specifies the number of cores to be used. Default is 1 (i.e. no parallel computing).
#' @param pb Logical argument to control if progress bar is shown. Default is \code{TRUE}.
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
#' X = example_space,
#' dimensions =  c("Dimension_1", "Dimension_2"),
#' group = "ID",
#' type = "mcp")
#' 
#' # mst rarefacted
#' rarefact_space_size(
#' X = example_space,
#' dimensions =  c("Dimension_1", "Dimension_2"),
#' group = "ID",
#' type = "mst")
#' }
#' @seealso \code{\link{rarefact_space_similarity}}, \code{\link{space_size_difference}}
#' @author Marcelo Araya-Salas \email{marcelo.araya@@ucr.ac.cr})
#'
#' @references {
#' Araya-Salas, M, & K. Odom. 2022, PhenotypeSpace: an R package to quantify and compare phenotypic trait spaces R package version 0.1.0.
#' }
# last modification on jan-2022 (MAS)

rarefact_space_size <- function(X, dimensions, group, n = NULL, replace = FALSE, seed = NULL, parallel = 1, pb = TRUE, iterations = 30, ...){
  
  obs.n <- min(table(X[, group]))
  
  if (!is.null(n)) {
    if (obs.n < n) {
      
      if (replace) message("'n' higher than the minimum sample size for at least 1 group, running rarefaction with replacement (replacement = TRUE)")
      replace <- TRUE
    } 
  } else 
    n <- obs.n
  
  X$...rownames <-  1:nrow(X) 
  
  # run iterations
  space_sizes_list <- warbleR:::pblapply_wrblr_int(1:iterations, cl = parallel, pb = pb, function(e){
    if (!is.null(seed))
      set.seed(seed + e)
    
    raref_indices <- unlist(lapply(sort(unique(X[, group])), function(x)
      sample(X$...rownames[X[, group] == x], n, replace = replace)
    ))
    
    overlaps <- space_size(X[raref_indices, ], group = group, dimensions = dimensions, pb = FALSE, parallel = 1, ...)
    
    return(overlaps)
  })
  
  space_sizes_mat <- do.call(cbind, lapply(space_sizes_list, '[', -c(1, 2)))
  
  results <- space_sizes_list[[1]][, 1:2]
  results$mean.size <- rowMeans(space_sizes_mat)
  results$min.size <- apply(X = space_sizes_mat, MARGIN = 1, FUN = min)  
  results$max.size <- apply(X = space_sizes_mat, MARGIN = 1, FUN = max)  
  results$sd.size <- apply(X = space_sizes_mat, MARGIN = 1, FUN = stats::sd)  
  
  return(results)
}
