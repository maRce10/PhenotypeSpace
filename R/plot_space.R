#' @title Plot bidimensional trait spaces
#'
#' @description \code{plot_space} plots bidimensional trait spaces   
#' @param X Data frame containing columns for the dimensions of the phenotypic space (numeric) of the total trait space. This is required so the extent of the plotting area represents the overall trait space in which the sub-space is found.  
#' @param dimensions Character vector of length 2 with the names of the columns containing the dimensions of the phenotypic space.
#' @param indices Numeric vector with the indices of the rows in 'X' to be used as sub-space for plotting. 
#' @param basecex Numeric vector of length 1 controlling the relative size of the axis labels and tick labels, legend and title. Legend and title are multiply by 1.5 (\code{basecex * 1.5}) to increase size compare to axis text.
#' @param title Character vector of length 1 to be used as the plot title. Default is \code{NULL}.
#' @param colors Character vector with the colors to use for density plotting. 2 values must be supplied if 'background.indices' is supplied. Default is \code{c("#3E4A89FF", "#35B779FF")}.
#' @param point.colors Character vector with the colors to use for point plotting. 2 values must be supplied if 'background.indices' is supplied. Default is the same as "colors".
#' @param point.alpha Numeric vector of length 1 >= 0 and <= 1 with the alpha value for color transparency. Default is 0.7. If 0 points are not plotted.
#' @param point.cex Numeric vector of length 1 controlling the relative size of the points. Default is 1. If 0 points are not plotted.
#' @param background.indices Numeric vector with the indices of the rows in 'X' to be used as background traits space for plotting. Points from 'indices' will be plotted on top of these points.
#' @param pch Either an integer specifying a symbol or a single character to be used as the default in plotting points. See \code{\link[graphics]{points}} for possible values and their interpretation. 
#' @param labels Character vector with the labels to be used in the legend. Not used if \code{legend.pos = NULL} or if 'background.indices' is not supplied. Default is \code{c("sub-space", "total space")}.
#' @param legend.pos Controls the position of the legend. Can take the following values: "bottomright", "bottom", "bottomleft", "left", "topleft", "top", "topright" (default), "right" and "center". If \code{NULL} the legend is not plotted.
#' @param density.alpha Numeric vector of length 1 >= 0 and <= 1 with the alpha value for color transparency to be used in the highest density regions. Lower density regions will gradually increase in transparency starting from the supplied value. Default is 0.6. If 0 densities are not plotted.
#' @param ... Additional arguments to be passed to \code{\link[graphics]{plot}} for plot customization.
#' @return A single panel plot in the active graphic device.
#' @export
#' @name plot_space
#' @details The function plots a sub-group of data (i.e. sub-space) within the overall trait space. The total trait space can also be plotted in the background. By default both points and kernel densities are shown. Graphs are returned in the active graphic device.   
#' @examples {
#' data("example_space")
#' 
#' # no background
#' plot_space(X = example_space, dimensions = c("dimension_1", "dimension_2"), 
#' indices = which(example_space$group == "G2"))
#' 
#' # add background
#' plot_space(X = example_space, dimensions = c("dimension_1", "dimension_2"), 
#' indices = which(example_space$group == "G2"), 
#' background.indices = which(example_space$group != "G2"))
#' 
#' # change legend labels
#' plot_space(X = example_space, dimensions = c("dimension_1", "dimension_2"), 
#' indices = which(example_space$group == "G2"), 
#' background.indices = which(example_space$group != "G2"), 
#' labels = c("G3", "trait space"))
#' 
#' # change legend position
#' plot_space(X = example_space, dimensions = c("dimension_1", "dimension_2"), 
#' indices = which(example_space$group == "G2"), 
#' background.indices = which(example_space$group != "G2"), 
#' labels = c("G3", "trait space"), legend.pos = "left")
#' 
#' # with title
#' plot_space(X = example_space, dimensions = c("dimension_1", "dimension_2"), 
#' indices = which(example_space$group == "G2"), 
#' background.indices = which(example_space$group != "G2"), 
#' labels = c("G3", "trait space"), legend.pos = "bottomleft", title = "G3")
#' }
#' @seealso \code{\link{distance_to_rectangular}}, \code{\link{rectangular_to_triangular}}
#' @author Marcelo Araya-Salas \email{marcelo.araya@@ucr.ac.cr})
#'
#' @references
#' Araya-Salas, M, K. Odom. & A. Rico-Guevara. 2022, PhenotypeSpace: an R package to quantify and compare phenotypic trait spaces R package version 0.1.0.
#' 
plot_space <-
  function(X,
           dimensions,
           indices,
           basecex = 1,
           title = NULL,
           colors = c("#3E4A89FF", "#35B779FF"),
           point.colors = colors,
           point.alpha = 0.7,
           point.cex = 1,
           background.indices = NULL,
           pch = 1,
           labels = c("sub-space", "total space"),
           legend.pos = "topright",
           density.alpha = 0.6, 
           ...
  ) {
  
  total_coors <- spatstat.geom::as.ppp(as.matrix(X[, dimensions]), c(range(X[, dimensions[1]]), range(X[, dimensions[2]])))
  total_space <- raster::raster(spatstat.explore::density.ppp(total_coors))
  
  xlim <- range(X[, dimensions[1]]) 
  ylim <-  range(X[, dimensions[2]])
  
  plot(x = X[, dimensions[1]], y = X[, dimensions[2]], col = "white", cex = basecex, xlab = dimensions[1], ylab= dimensions[2], xlim = xlim, ylim = ylim, cex.lab = basecex, xaxs="i", yaxs="i", ...) 
  
  # add background group
  if (!is.null(background.indices)){
    
    # get points  
    Y_bg <- X[background.indices, ]
    
    if (nrow(Y_bg) > 1 & point.alpha > 0)
      points(Y_bg[, dimensions[1]], Y_bg[, dimensions[2]], col = adjustcolor(point.colors[2], alpha.f = point.alpha), pch = pch, cex = point.cex)
    
    if (nrow(Y_bg) > 4){
      coors_bg <- spatstat.geom::as.ppp(as.matrix(Y_bg[, dimensions]), c(range(Y_bg[, dimensions[1]]), range(Y_bg[, dimensions[2]])))
      
      raster_dens_bg <- raster::raster(spatstat.explore::density.ppp(coors_bg))
      
      cols_bg <- sapply(1:10, function(x) adjustcolor(col =  colorRampPalette(c("white", colors[2]))(10)[x], alpha.f = seq(0.1, density.alpha, length.out = 10)[x]))
      
      raster::image(raster_dens_bg, add = TRUE, col = cols_bg)
    }
    
  }
  
  Y <- X[indices, ]
  
  if (nrow(Y) > 1 & point.alpha > 0)
    points(Y[, dimensions[1]], Y[, dimensions[2]], col = adjustcolor(point.colors[1], alpha.f = point.alpha), pch = pch, cex = point.cex)
  
  if (nrow(Y) > 4){
    coors_focal <- spatstat.geom::as.ppp(as.matrix(Y[, dimensions]), c(range(Y[, dimensions[1]]), range(Y[, dimensions[2]])))
    
    raster_dens_focal <- raster::raster(spatstat.explore::density.ppp(coors_focal))
        
    cols_focal <- sapply(1:10, function(x) adjustcolor(col =  colorRampPalette(c("white", colors[1]))(10)[x], alpha.f = seq(0.1, density.alpha, length.out = 10)[x]))
    
    if (density.alpha > 0)
    raster::image(x = raster_dens_focal, add = TRUE, col = cols_focal)
  }
  
  usr <- par()$usr
  
  # legend  
  if (!is.null(background.indices) & !is.null(legend.pos))
    legend(legend.pos, col = colors, legend = labels, pch = c(pch, pch), cex = basecex * 1.3, xjust = 0, bg = adjustcolor("white", alpha.f = 0.5), box.col = "gray40", box.lwd = 0.8)
  
  if (!is.null(title))
  title(title, cex.main = basecex * 1.2)
}
