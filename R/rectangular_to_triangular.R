#' @title Convert rectangular pairwise matrices to triangular matrices
#'
#' @description \code{rectangular_to_triangular} converts rectangular pairwise matrices as those output by many PhenotypeSpace functions into triangular pairwise matrices. 
#' @param X Data frame containing three columns. The first two columns must contain group labels which will appear as rows (1 column) and column names (2 column) in the output triangular matrix. The third column (and fourth column if \code{symmetric = FALSE}) must have the numeric values to be included in the output triangular matrix. 
#' @param distance Logical argument to control if the input data contains pairwise distances (dissimilarities) or similarities. If \code{TRUE} then diagonal values are filled with 0, otherwise they are filled with 1. Note that diagonal values can be set with \code{\link[base]{diag}}.
#' @param symmetric Logical argument to define if values are duplicated on both off-diagonal triangles (a symmetric triangular matrix, \code{symmetric = TRUE}, default) or each triangle is filled with values from different columns (a non-symmetric triangular matrix, \code{symmetric = FALSE}). In the latter the upper triangle is filled with the first column and the lower triangle with the second column. In this case, a fourth column with numeric values should be supplied.    
#' @return A pairwise triangular matrix in which labels from the first group column in 'X' are shown in the columns and labels from the second group are shown in the rows. If \code{symmetric = FALSE} the same information is shown below and above the diagonal.
#' @export
#' @name rectangular_to_triangular
#' @details The function converts rectangular pairwise matrices as those output by many PhenotypeSpace functions into triangular pairwise matrices. It takes a data frame in which each observation (row) contains the pairwise value and related labels of the 'groups' being compared. The first two columns must contain group labels which will appear as rows (1 column) and column names (2 column) in the output triangular matrix. The third column (and fourth column if \code{symmetric = FALSE}) must have the numeric values to be included in the output triangular matrix.  
#' @examples {
#' # load data
#' data("example_space")
#' 
#' # get proportion of space that overlaps 
#' prop_overlaps <- space_similarity(
#'  formula = group ~ dimension_1 + dimension_2,
#'  data = example_space,
#'  method = "proportional.overlap")
#' 
#' # get symmetric triangular matrix
#' rectangular_to_triangular(prop_overlaps)
#' 
#' # get minimum convex polygon overlap for each group (non-symmetric)
#' mcp_overlaps <- space_similarity(
#'  formula = group ~ dimension_1 + dimension_2,
#'  data = example_space,
#'  method = "mcp.overlap")
#' 
#' # get a non-symmetric triangular matrix
#' rectangular_to_triangular(mcp_overlaps, symmetric = FALSE)
#' }
#' @seealso \code{\link{distance_to_rectangular}}, \code{\link{binary_triangular_matrix}}
#' @author Marcelo Araya-Salas \email{marcelo.araya@@ucr.ac.cr})
#'
#' @references
#' Araya-Salas, M, & K. Odom. 2022, PhenotypeSpace: an R package to quantify and compare phenotypic trait spaces R package version 0.1.0.

rectangular_to_triangular <- function(X, distance = TRUE, symmetric = TRUE){
  
  if (!symmetric & ncol(X) < 4) 
    stop2("'X' must have at least 2 numeric columns when `symmetric = FALSE`")

  # get groups
  groups <- sort(unique(c(X[, 1], X[, 2])))
  
  mat <- matrix(nrow = length(groups), ncol = length(groups))
  
  colnames(mat) <- rownames(mat) <- groups
  
  for (i in 1:nrow(X)) {
    
    mat[colnames(mat) == X[i, 1], colnames(mat) ==X[i, 2]] <- X[i, 3]
    
    
    mat[colnames(mat) == X[i, 2], colnames(mat) ==X[i, 1]] <- if (symmetric) X[i, 3] else 
        X[i, 4]
  }
  
  # fill out diagonal
  for (i in groups)
    mat[i, i] <- if (distance) 0 else 1
  
  return(mat)
  
}
