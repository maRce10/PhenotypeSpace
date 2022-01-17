#' @title Estimates pairwise size differences of phenotypic spaces
#'
#' @description \code{space_size_difference}
#' @usage space_size_difference(X, dimensions, group, parallel = 1, type = "mcp", 
#' pb = TRUE, outliers = 0.95, proportional = FALSE, ...)
#' @param X Data frame containing columns for the dimensions of the phenotypic space (numeric) and a categorical or factor column with group labels. 
#' @param dimensions Character vector with the names of columns containing the dimensions of the phenotypic space.
#' @param group Character vector with the name of the column (character or factor) containing group labels.
#' @param parallel Integer vector of length 1. Controls whether parallel computing is applied. It specifies the number of cores to be used. Default is 1 (i.e. no parallel computing).
#' @param type
#' @param pb Logical argument to control if progress bar is shown. Default is \code{TRUE}.
#' @param outliers
#' @param proportional
#' @param ...
#' @return
#' @export
#' @name space_size_difference
#' @details   
#' @examples {
#' # load data
#' data("example_space")
#' 
#' # MCP size (try with more iterations on your own data)
#' mcp_size <- space_size_difference(
#' X = example_space,
#' dimensions =  c("Dimension_1", "Dimension_2"),
#' group = "ID",
#' type = "mcp", 
#' iterations = 5)
#' 
#' # MST size
#' mcp_size <- space_size_difference(
#' X = example_space,
#' dimensions =  c("Dimension_1", "Dimension_2"),
#' group = "ID",
#' type = "mst", 
#' iterations = 5)
#' }
#' @seealso \code{\link{space_size}}, \code{\link{space_overlap}}, \code{\link{rarefact_space_size_difference}}
#' @author Marcelo Araya-Salas \email{marcelo.araya@@ucr.ac.cr})
#'
#' @references {
#' Araya-Salas, M, & K. Odom. 2022, PhenotypeSpace: quantifying phenotypic trait spaces. R package version 0.1.0.
#' }
# last modification on jan-2022 (MAS)

space_size_difference <- function(X, dimensions, group, parallel = 1, type = "mcp", pb = TRUE, outliers = 0.95, proportional = FALSE, ...){
  
  group_combs <- t(utils::combn(unique(X[, group]), 2))
  
  propspace_size_diff_list <- lapply(1:nrow(group_combs), function(x){
    
    W <- X[X[, group] == group_combs[x, 1], ]
    Z <- X[X[, group] == group_combs[x, 2], ]     
    
    # combined both
    both <- rbind(W, Z)
    both[, group] <- as.character(both[, group])
    both[, group] <- "both"
    
    Y <- rbind(W, Z, both)
    
    sizes <- space_size(X = Y, dimensions = dimensions, group = group,  parallel = 1, type = type, pb = FALSE, outliers = outliers, proportional = proportional)
    
    size.diff <- sizes$size[sizes$group == group_combs[x, 1]] - sizes$size[sizes$group == group_combs[x, 2]]
    
    out <- data.frame(group.1 = group_combs[x, 1], group.2 = group_combs[x, 2], size.difference = size.diff)
    
    return(out)
    
  })
  
  result <- do.call(rbind, propspace_size_diff_list)
  
  
  return(result)  
}
