#' @title Example bi-dimensional space data
#'
#' @usage data(example_space)
#'
#' @description \code{example_space} is a data frame with 1550 observations of a simulated phenotypic bi-dimensional space including 5 groups.
#'
#' @example {
#' #plot space 
#' xs <- tapply(example_space$Dimension_1, example_space$ID, mean)
#' ys <- tapply(example_space$Dimension_2, example_space$ID, mean)
#' plot(example_space[, c("Dimension_1", "Dimension_2")], 
#'    col = example_space$color, pch = 20, cex = 1.8)
#' text(xs, ys, labels = names(xs), cex = 2.5)
#' }
#'
#' @source Marcelo Araya Salas, PhenotypeSpace
"example_space"
