#' @title Calculates rarefacted space overlaps
#'
#' @description \code{rarefact_space_overlap}
#' @usage rarefact_space_overlap(X, dimensions, group, n = NULL, replace = FALSE, 
#' seed = NULL, parallel = 1, pb = TRUE, iterations = 30, ...)
#' @param X Data frame containing columns for the dimensions of the phenotypic space (numeric) and a categorical or factor column with group labels. 
#' @param dimensions Character vector with the names of columns containing the dimensions of the phenotypic space.
#' @param group Character vector with the name of the column (character or factor) containing group labels.
#' @param n Integer vector of length 1 indicating the number of samples to be use for rarefaction (i.e. how many samples per group will be gather at each iteration). Default is the minimum sample size across groups.
#' @param replace Logical argument to control if sampling is done with replacement. Default is \code{FALSE}.
#' @param seed Integer vector of length 1 setting the seed (see \code{\link[base]{set.seed}}). If used results should be the same on different runs, so it makes them replicable.
#' @param parallel Integer vector of length 1. Controls whether parallel computing is applied. It specifies the number of cores to be used. Default is 1 (i.e. no parallel computing).
#' @param pb Logical argument to control if progress bar is shown. Default is \code{TRUE}.
#' @param iterations Integer vector of length 1. Controls how the number of times the rarefaction routine is iterated. Default is 30.
#' @param ...
#' @return
#' @export
#' @name rarefact_space_overlap
#' @details   
#' @examples {
#' # load data
#' data("example_space")
#' 
#' # get proportion of space that overlaps (try with more iterations on your own data)
#' prop_overlaps <- rarefact_space_overlap(
#' X = example_space,
#' dimensions =  c("Dimension_1", "Dimension_2"),
#' group = "ID",
#' type = "proportional", 
#' iterations = 5)
#' 
#' # get minimum convex polygon overlap for each group (non-symmetric)
#' mcp_overlaps <- rarefact_space_overlap(
#' X = example_space,
#' dimensions =  c("Dimension_1", "Dimension_2"),
#' group = "ID",
#' iterations = 5)
#' 
#' # convert to non-symmetric triangular matrix
#' rectangular_to_triangular(mcp_overlaps, symmetric = FALSE)
#' }
#' @seealso \code{\link{space_overlap}}, \code{\link{rectangular_to_triangular}}
#' @author Marcelo Araya-Salas \email{marcelo.araya@@ucr.ac.cr})
#'
#' @references {
#' Araya-Salas, M, & K. Odom. 2022, PhenotypeSpace: quantifying phenotypic trait spaces. R package version 0.1.0.
#' }
# last modification on jan-2022 (MAS)

rarefact_space_overlap <- function(X, dimensions, group, n = NULL, replace = FALSE, seed = NULL, parallel = 1, pb = TRUE, iterations = 30, ...){
  
  obs.n <- min(table(X[, group]))
  
  if (!is.null(n)) {
    if (obs.n < n) {
      
      if (replace) message("'n' higher than the minimum sample size for at least 1 group, running rarefaction with replacement (replace = TRUE)")
      replace <- TRUE
    } 
  } else 
    n <- obs.n
  
  X$...rownames <-  1:nrow(X) 
  
  # run iterations
  space_overlaps_list <- warbleR:::pblapply_wrblr_int(1:iterations, cl = parallel, pb = pb, function(e){
    if (!is.null(seed))
      set.seed(seed + e)
    
    raref_indices <- unlist(lapply(unique(X[, group]), function(x)
      sample(X$...rownames[X[, group] == x], n, replace = replace)
    ))
    
    overlaps <- space_overlap(X[raref_indices, ], group = group, dimensions = dimensions, pb = FALSE, parallel = 1, ...)
    
    return(overlaps)
  })
  
  space_overlaps_mat <- do.call(cbind, lapply(space_overlaps_list, '[', -c(1, 2)))
  
  results <- space_overlaps_list[[1]][, 1:2]
  
  # get type from arguments set in the call
  call.argms <- as.list(base::match.call())[-1]

  if (any(names(call.argms) == "type"))
    type <- call.argms$type else
      type <- "mcp" # default value in space_overlap()
  
  # create output data frame
  if (type %in% c("mean.density", "mean.mcp", "std.centroid.dist", "proportional")){
  results$mean.overlap <- rowMeans(space_overlaps_mat)
  results$min.overlap <- apply(X = space_overlaps_mat, MARGIN = 1, FUN = min)  
  results$max.overlap <- apply(X = space_overlaps_mat, MARGIN = 1, FUN = max)  
  results$sd.overlap <- apply(X = space_overlaps_mat, MARGIN = 1, FUN = stats::sd)
  } else {
    results$mean.overlap.1in2 <- rowMeans(space_overlaps_mat[ , seq(1, (iterations * 2) - 1, 2)])
    results$mean.overlap.2in1 <- rowMeans(space_overlaps_mat[ , seq(2, (iterations * 2), 2)])
    results$min.overlap.1in2 <- apply(X = space_overlaps_mat[ , seq(1, (iterations * 2) - 1, 2)], MARGIN = 1, FUN = min)
    results$min.overlap.2.in1 <- apply(X = space_overlaps_mat[ , seq(2, (iterations * 2), 2)], MARGIN = 1, FUN = min)  
    results$max.overlap.1in2 <- apply(X = space_overlaps_mat[ , seq(1, (iterations * 2) - 1, 2)], MARGIN = 1, FUN = max)
    results$max.overlap.2.in1 <- apply(X = space_overlaps_mat[ , seq(2, (iterations * 2), 2)], MARGIN = 1, FUN = max)  
    results$sd.overlap.1in2 <- apply(X = space_overlaps_mat[ , seq(1, (iterations * 2) - 1, 2)], MARGIN = 1, FUN = sd)
    results$sd.overlap.2.in1 <- apply(X = space_overlaps_mat[ , seq(2, (iterations * 2), 2)], MARGIN = 1, FUN = sd)  
  }
  
  return(results)
}
