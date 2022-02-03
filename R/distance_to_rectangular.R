#' @title Convert pairwise distance matrices to rectangular matrices 
#'
#' @description \code{distance_to_rectangular} converts binary triangular matrices to rectangular matrices using Multidimensional Scaling. 
#' @usage distance_to_rectangular(distance.matrix, labels = names(distance.matrix), 
#' n.dimensions = 2, metric = TRUE, ...)
#' @param distance.matrix Distance matrix (i.e. object of class 'dist'). Can be created using the function \code{\link[stats]{dist}} or converted to using  \code{\link[stats]{as.dist}}. 
#' @param labels Character vector or factor containing labels to be used for rows/columns in the output data frame. Default is \code{names(distance.matrix)}. Must be the same length as the number of observations in 'distance.matrix'.
#' @param n.dimensions Integer vector of length 1 indicating the number of of dimensions to represent distances in a new space. Default is 2.
#' @param metric Logical argument to control if Metric (a.k.a. Classical, \code{TRUE}, default) or Non-Metric MUltidimensional Scaling (\code{FALSE}) is used to project in a new n-dimension space. Non-Metric MDS is conducted using the function \code{\link[MASS]{isoMDS}} while Classical MDS uses the function \code{\link[stats]{cmdscale}}. So yes, it is a silly wrapper over those 2 functions.
#' @param ... Additional arguments to be passed to \code{\link[MASS]{isoMDS}}.
#' @return  A data frame with the new dimensions representing the position of observations in a new n-dimension space. If \code{metric = FALSE} the output data frame is embedded in a list that also includes the stress value.
#' @export
#' @name distance_to_rectangular
#' @details It is a silly wrapper over 2 multidimensional scaling functions (\code{\link[MASS]{isoMDS}} and \code{\link[stats]{cmdscale}}) that simplifies the calculation of Multidimensional Scaling and formatting of its output to be used with other functions in the package.
#' @examples {
#' data("example_space")
#' 
#' dist_example <- dist(example_space[example_space$ID %in% c("G1", "G2"), 
#' c("Dimension_1", "Dimension_2")])
#' 
#' # convert into a 2-dimension space
#' rect_example <- distance_to_rectangular(distance.matrix = dist_example,
#' metric = TRUE)
#' 
#' head(rect_example)
#' 
#' # convert into a 2-dimension space with non-metric MDS
#' rect_example <- distance_to_rectangular(distance.matrix = dist_example, 
#' metric = FALSE, maxit = 3)
#' }
#' @seealso \code{\link{distance_to_rectangular}}, \code{\link{rectangular_to_triangular}}
#' @author Marcelo Araya-Salas \email{marcelo.araya@@ucr.ac.cr})
#'
#' @references {
#' Araya-Salas, M, & K. Odom. 2022, PhenotypeSpace: an R package to quantify and compare phenotypic trait spaces R package version 0.1.0.
#' }
# last modification on jan-2022 (MAS)

distance_to_rectangular <- function(distance.matrix, labels = names(distance.matrix), n.dimensions = 2, metric = TRUE, ...)
{
  
  # mds for 2 dimensions
  mds <- stats::cmdscale(distance.matrix, k = n.dimensions)
  
  if (length(labels) != attr(distance.matrix, "Size"))
    stop("'labels' should have the same length as the number of observations in 'distance.matrix' (see 'attributes(distance.matrix)$Size')")
  
  # non-metric
  if (!metric){
    cat("Calculating non-metric multidimensional scaling: ")
  nm_mds <- MASS::isoMDS(distance.matrix, y = mds, k = n.dimensions, ...)
  
  # extract results
  mds <- nm_mds$points
  stress <- nm_mds$stress
  }
  
  # add stress if distance matrix was supplied
  
  # put in a data frame
  output <- data.frame(labels, mds)
  names(output) <- c("labels", "MDS1", "MDS2")

  if (!metric)
    output <- list(stress = stress, points = output)
    
  return(output)
} 
