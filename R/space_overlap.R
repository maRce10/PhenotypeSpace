# overlap.type can be used to define how overlap is calculated. Options:
# "mean" = average of the overlaps of group 1 vs group 2 and group 2 vs group 1
# assymetric =  group 1 vs group 2 and group 2 vs group 1
# "proportional" = proportion of the joint area of both groups that overlaps (overlapped area / total area of both groups)
# standardize centroid distance = "std.centroid.dist"
#' @title Estimates pairwise overlaps of phenotype spaces
#'
#' @description \code{space_overlap}
#' @usage space_overlap(X, dimensions, group, parallel = 1, type = "mcp", 
#' pb = TRUE, outliers = 0.95, iterations = 30)
#' @param X Data frame containing columns for the dimensions of the phenotypic space (numeric) and a categorical or factor column with group labels. 
#' @param dimensions Character vector with the names of the columns containing the dimensions of the phenotypic space.
#' @param group Character vector with the name of the column (character or factor) containing group labels.
#' @param parallel Integer vector of length 1. Controls whether parallel computing is applied. It specifies the number of cores to be used. Default is 1 (i.e. no parallel computing).
#' @param type 
#' @param pb Logical argument to control if progress bar is shown. Default is \code{TRUE}.
#' @param outliers
#' @param iterations Integer vector of length 1. Controls how the number of times the rarefaction routine is iterated. Default is 30.
#' @return
#' @export
#' @name space_overlap
#' @details   
#' @examples {
#' # load data
#' data("example_space")
#' 
#' # get proportion of space that overlaps 
#' prop_overlaps <- space_overlap(
#' X = example_space,
#' dimensions =  c("Dimension_1", "Dimension_2"),
#' group = "ID",
#' type = "proportional")
#' 
#' #' # get symmetric triangular matrix
#' rectangular_to_triangular(prop_overlaps)
#' 
#' # get minimum convex polygon overlap for each group (non-symmetric)
#' mcp_overlaps <- space_overlap(
#' X = example_space,
#' dimensions =  c("Dimension_1", "Dimension_2"),
#' group = "ID",
#' type = "mcp")
#' 
#' # convert to non-symmetric triangular matrix
#' rectangular_to_triangular(mcp_overlaps, symmetric = FALSE)
#' }
#' @seealso \code{\link{rarefact_space_overlap}}, \code{\link{rectangular_to_triangular}}
#' @author Marcelo Araya-Salas \email{marcelo.araya@@ucr.ac.cr})
#'
#' @references {
#' Araya-Salas, M, & K. Odom. 2022, PhenotypeSpace: quantifying phenotypic trait spaces. R package version 0.1.0.
#' }
# last modification on jan-2022 (MAS)

space_overlap <- function(X, dimensions, group, parallel = 1, type = "mcp", pb = TRUE, outliers = 0.95, iterations = 30) {
  
  if (type == "std.centroid.dist") {
    
    X[, dimensions[1]] <- scale(X[, dimensions[1]])
    X[, dimensions[2]] <- scale(X[, dimensions[2]])
  
    } else {
  # force dimensions to range 0-1
  X[, dimensions[1]] <- (X[, dimensions[1]] - min( X[, dimensions[1]])) / max((X[, dimensions[1]] - min( X[, dimensions[1]])))
  X[, dimensions[2]] <- (X[, dimensions[2]] - min( X[, dimensions[2]])) / max((X[, dimensions[2]] - min( X[, dimensions[2]])))
  
    }
  
  # split in a list of vectors
  X_l <- split(x = X, f = X[, group])
  
  # stop if too small sample sizes
  if (min(sapply(X_l, nrow)) < 2)
    stop("There is at least one group with less than 2 observations which is the minimum needed for overlap estimation")
  
  # get densities
  if (type %in% c("density", "mean.density")){
    
    total_coors <-  spatstat.geom::as.ppp(as.matrix(X[, dimensions]), c(range(X[, dimensions[1]]), range(X[, dimensions[2]])))
    total_space <-  raster::raster(spatstat.core::density.ppp(total_coors))
    
    input_data <- lapply(X_l, function(Y, dims = dimensions){
      
      coors <- spatstat.geom::as.ppp(as.matrix(Y[, dimensions]), c(range(Y[, dimensions[1]]), range(Y[, dimensions[2]])))
      
      raster_dens <- raster::raster(spatstat.core::density.ppp(coors))
      
      raster_dens <- raster::crop(raster::extend(raster_dens, total_space), total_space)
      
      # convert to range 0-1
      raster::values(raster_dens) <- (raster::values(raster_dens) - min(raster::values(raster_dens), na.rm = TRUE)) 
      raster::values(raster_dens) <- raster::values(raster_dens)/ max(raster::values(raster_dens), na.rm = TRUE)
      
      # keep 95% interval
      raster::values(raster_dens)[raster::values(raster_dens) < 1 - outliers] <- NA 
      
      # same resolution and extend as total space
      raster_dens <- raster::resample(raster_dens, total_space)
      raster::extent(raster_dens)  <- c(0,1, 0, 1)
      
      return(raster_dens)
    })
    
    names(input_data) <- names(X_l)
    
  } else input_data <- X_l
  
  
  # function to calculate areas
  ovlp_fun <- function(W, Z, tp, dims) {
    
    # get area
    # kernel area
    if (tp %in% c("density", "mean.density")){
      
      # convert NAs to 0
      raster::values(W)[is.na(raster::values(W))] <- 0
      raster::values(Z)[is.na(raster::values(Z))] <- 0
      
      # 1 vs 2
      wgt_ovlp.WvZ <- sum((raster::values(W) * raster::values(Z)) / raster::values(Z), na.rm = TRUE) / sum(raster::values(Z), na.rm = TRUE)
      
      # 2 vs 1
      wgt_ovlp.ZvW <- sum((raster::values(Z) * raster::values(W)) / raster::values(W), na.rm = TRUE) / sum(raster::values(W), na.rm = TRUE) 
      
      # convert to 1 if higher than 1
      if (wgt_ovlp.WvZ > 1) wgt_ovlp.WvZ <- 1
      if (wgt_ovlp.ZvW > 1) wgt_ovlp.ZvW <- 1
      
      
      # put results in am matrix     
      out <- matrix(c(wgt_ovlp.WvZ, wgt_ovlp.ZvW), nrow = 1)     
    } 
    
    # distance among points
    if (tp == "distance"){
      
      U <- rbind(W, Z)
      
      dists <- as.matrix(stats::dist(U[ , dims]))
      
      dist.1v2 <- dists[U[, group] == W[, group][1], U[, group] == Z[, group][1]]
      
      out <- matrix(mean(dist.1v2), nrow = 1)
    } 
    
    if (tp %in% c("mcp", "mean.mcp")){
      
      sp::coordinates(W) <- stats::as.formula(paste("~", paste(dims, collapse = "+")))
      sp::coordinates(Z) <- stats::as.formula(paste("~", paste(dims, collapse = "+")))
      
      
      # get intersect
      intrsct <- raster::intersect(adehabitatHR::mcp(W), adehabitatHR::mcp(Z))
      if (!is.null(intrsct))
      {
        intrsctspace_size <- raster::area(intrsct)
        
        ovlp1in2 <- intrsctspace_size / raster::area(adehabitatHR::mcp(W))
        ovlp2in1 <- intrsctspace_size / raster::area(adehabitatHR::mcp(Z))
      }      else ovlp1in2 <- ovlp2in1 <- 0
      
      out <- matrix(c(ovlp2in1, ovlp1in2), nrow = 1)  
    }
    
    if (tp == "proportional") {
      
      sp::coordinates(W) <- stats::as.formula(paste("~", paste(dims, collapse = "+")))
      sp::coordinates(Z) <- stats::as.formula(paste("~", paste(dims, collapse = "+"))) 
      
      # both combined
      Y <- rbind(W, Z)

      # get intersect
      intrsct <- raster::intersect(adehabitatHR::mcp(W, percent = outliers * 100), adehabitatHR::mcp(Z, percent = outliers * 100))
      
      if (!is.null(intrsct))
      {
        intrsctspace_size <- raster::area(intrsct)
        totalspace_size <- raster::area(adehabitatHR::mcp(Y, percent = outliers * 100))
        ovlp1in2 <- ovlp2in1 <- intrsctspace_size / totalspace_size
        
        }      else ovlp1in2 <- ovlp2in1 <- 0
      
      out <- matrix(c(ovlp2in1, ovlp1in2), nrow = 1)  
    }
   
    if (tp == "std.centroid.dist") {
      
      # both combined
      Y <- rbind(as.data.frame(W), as.data.frame(Z))
  
      frmla <- stats::as.formula(paste("cbind(", dims[1], ",", dims[2], ") ~ ",group))
      centroid.dist <- stats::dist(stats::aggregate(formula = frmla, data = Y, FUN = mean)[, -1])
      
      ovlp1in2 <- ovlp2in1 <- centroid.dist
      
      out <- matrix(c(ovlp2in1, ovlp1in2), nrow = 1)  
    }
    
    return(out)
  }

  
  # get all combinations to get pairwise overlaps
  group_combs <- t(utils::combn(sort(unique(X[, group])), 2))
  
  # calculate all areas
  ovlps_l <- warbleR:::pblapply_wrblr_int(1:nrow(group_combs), pbar = pb,  cl = parallel, function(i, data = input_data, gc = group_combs, dims = dimensions, typ = type) {
    
    if (type %in% c("density", "mean.density")){
      W <- data[[which(names(data) == gc[i, 1])]]
      Z <- data[[which(names(data) == gc[i, 2])]]
    } 
    
    if (type %in% c("distance", "mcp", "mean.mcp", "proportional", "std.centroid.dist")){
      W <- X_l[[which(names(X_l) == gc[i, 1])]]
      Z <- X_l[[which(names(X_l) == gc[i, 2])]]
    }
    
    suppressWarnings(overlaps <- ovlp_fun(W, Z, typ, dims))

    # put in a data frame
    
    out_df <- if (type %in% c("mean.density", "mean.mcp", "std.centroid.dist", "proportional", "distance")) 
      data.frame(group.1 = gc[i, 1], group.2 = gc[i, 2], overlap = mean(overlaps)) else
        data.frame(group.1 = gc[i, 1], group.2 = gc[i, 2], overlap.1in2 = overlaps[2], overlap.2in1 = overlaps[1]) 
    
    return(out_df)
  }) 
  
  # put in a data frame
  space_ovlp <- do.call(rbind, ovlps_l)

  return(space_ovlp)
}

