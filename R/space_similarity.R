#' @title Pairwise similarities of phenotype spaces
#'
#' @description \code{space_similarity} estimate pairwise similarities of phenotype spaces
#' @usage space_similarity(X, dimensions, group, parallel = 1, type = "mcp.overlap", 
#' pb = TRUE, outliers = 0.95, pairwise.scale = FALSE, distance.method = "Euclidean")
#' @param X Data frame containing columns for the dimensions of the phenotypic space (numeric) and a categorical or factor column with group labels. 
#' @param dimensions Character vector with the names of the columns containing the dimensions of the phenotypic space. For similarity metrics using overlap 2 dimensions are required. Distance metrics can used any number of dimensions.
#' @param group Character vector with the name of the column (character or factor) containing group labels.
#' @param parallel Integer vector of length 1. Controls whether parallel computing is applied. It specifies the number of cores to be used. Default is 1 (i.e. no parallel computing).
#' @param type Character vector of length 1. Controls the type of (di)similarity metric to be compare the phenotypic sub-spaces of two groups at the time. Seven metrics are available which quantify as pairwise sub-space overlap ('similarity') or pairwise distance between sub-spaces ('dissimilarity'):
#' \itemize{
#'  \item \code{density.overlap}: proportion of the phenotypic sub-spaces area that overlap, taking into account the irregular densities of the sub-spaces. Two groups that share their higher density areas will be more similar than similar sub-spaces that only share their lower density areas. Two values are supplied as the proportion of the space of A that overlaps B is not necessarily the same as the proportion of B that overlaps A. Similarity metric (higher values means more similar). The minimum sample size (per group) must be 6 observations. 
#'  \item \code{mean.density.overlap}: similar to 'density.overlap' but the two values are merged into a single pairwise mean overlap. Similarity metric (higher values means more similar). The minimum sample size (per group) must be 6 observations.
#'  \item \code{mcp.overlap}: proportion of the phenotypic sub-spaces area that overlap, in which areas are calculated as the minimum convex polygon of all observations for each su-space. Two values are supplied as the proportion of the space of A that overlaps B is not necessarily the same as the proportion of B that overlaps A. Similarity metric (higher values means more similar). The minimum sample size (per group) must be 5 observations.
#'  \item \code{mean.mcp.overlap}: similar to 'mcp.overlap' but the two values are merged into a single pairwise mean overlap. Similarity metric (higher values means more similar). The minimum sample size (per group) must be 5 observations.
#'  \item \code{proportional.overlap}: proportion of the joint area of both sub-spaces that overlaps (overlapped area / total area of both groups). Sub-space areas are calculated as the minimum convex polygon. Similarity metric (higher values means more similar). The minimum sample size (per group) must be 5 observations.
#'  \item \code{distance}: mean euclidean pairwise distance between all observations of the compared sub-spaces. Dissimilarity metric (higher values means less similar). The minimum sample size (per group) must be 1 observation.
#'  \item \code{centroid.distance}: euclidean distance between the centroid of the compared sub-spaces. Dissimilarity metric (higher values means less similar). The minimum sample size (per group) must be 1 observation.
#'  }
#' @param pb Logical argument to control if progress bar is shown. Default is \code{TRUE}.
#' @param outliers Numeric vector of length 1. A value between 0 and 1 controlling the proportion of outlier observations to be excluded. Outliers are determined as those farthest away from the sub-space centroid.
#' @param pairwise.scale Logical argument to control if pairwise phenotypic spaces are scaled (i.e. z-transformed) prior to similarity estimation. If so (\code{TRUE}) similarities are decoupled from the size of the global phenotypic space. Useful to compare similarities coming from different phenotypic spaces. Default is \code{FALSE}. Not available for 'density.overlap' and 'mean.density.overlap'.
#' @param distance.method Character vector of length 1 indicating the method to be used for measuring distances (hence only applicable when distances are calculated). Default is 'Euclidean'. All distance and similarity measures available in \code{\link[proxy]{dist}} can be used (but note that not all of them apply to continuous data). Check available metrics by running \code{summary(pr_DB)}. If a similarity measure is used similarities are converted to distances.
#' @return A data frame containing the similarity metric for each pair of groups. If the similarity metric is not symmetric (e.g. the proportional area of A that overlaps B is not necessarily the same as the area of B that overlaps A, see \code{\link{space_similarity}}) separated columns are supplied for the two comparisons.  
#' @export
#' @name space_similarity
#' @details The function quantifies pairwise similarity between phenotypic sub-spaces. Similarity is evaluated as the overlap (similarity) or distance (dissimilarity) between group.    
#' @examples {
#' # load data
#' data("example_space")
#' 
#' # get proportion of space that overlaps 
#' prop_overlaps <- space_similarity(
#' X = example_space,
#' dimensions =  c("Dimension_1", "Dimension_2"),
#' group = "ID",
#' type = "proportional.overlap")
#' 
#' #' # get symmetric triangular matrix
#' rectangular_to_triangular(prop_overlaps)
#' 
#' # get minimum convex polygon overlap for each group (non-symmetric)
#' mcp_overlaps <- space_similarity(
#' X = example_space,
#' dimensions =  c("Dimension_1", "Dimension_2"),
#' group = "ID",
#' type = "mcp.overlap")
#' 
#' # convert to non-symmetric triangular matrix
#' rectangular_to_triangular(mcp_overlaps, symmetric = FALSE)
#' 
#' # check available distance measures 
#' summary(pr_DB)
#' 
#' # get eculidean distances (default)
#' area_dist <- space_similarity(
#' X = example_space,
#' dimensions =  c("Dimension_1", "Dimension_2"),
#' group = "ID",
#' type = "distance",
#' distance.method = "Euclidean")
#'
#' # get Canberra distances 
#' area_dist <- space_similarity(
#' X = example_space,
#' dimensions =  c("Dimension_1", "Dimension_2"),
#' group = "ID",
#' type = "distance",
#' distance.method = "Canberra")
#' }
#' @seealso \code{\link{rarefact_space_similarity}}, \code{\link{space_size_difference}}
#' @author Marcelo Araya-Salas \email{marcelo.araya@@ucr.ac.cr})
#'
#' @references {
#' Araya-Salas, M, & K. Odom. 2022, PhenotypeSpace: quantifying phenotypic trait spaces. R package version 0.1.0.
#' }
# last modification on jan-2022 (MAS)

space_similarity <- function(X, dimensions, group, parallel = 1, type = "mcp.overlap", pb = TRUE, outliers = 0.95, pairwise.scale = FALSE, distance.method = "Euclidean") {
  
  if (!type %in% c("density.overlap", "mean.density.overlap", "mcp.overlap", "mean.mcp.overlap", "proportional.overlap", "distance", "centroid.distance"))
    stop("Unsupported 'type' declared")
  # 
  # if (type == "centroid.distance") {
  #   
  #   X[, dimensions[1]] <- scale(X[, dimensions[1]])
  #   X[, dimensions[2]] <- scale(X[, dimensions[2]])
  # 
  #   } else {
  # # force dimensions to range 0-1
  # X[, dimensions[1]] <- (X[, dimensions[1]] - min( X[, dimensions[1]])) / max((X[, dimensions[1]] - min( X[, dimensions[1]])))
  # X[, dimensions[2]] <- (X[, dimensions[2]] - min( X[, dimensions[2]])) / max((X[, dimensions[2]] - min( X[, dimensions[2]])))
  # 
  #   }
  
  # split in a list of vectors
  X_l <- split(x = X, f = X[, group])
  
  # stop if too small sample sizes
  if (min(sapply(X_l, nrow)) < 2)
    stop("There is at least one group with less than 2 observations which is the minimum needed for overlap estimation")
  
  # get densities
  if (type %in% c("density.overlap", "mean.density.overlap")){
    
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
  ovlp_fun <- function(W, Z, tp, dims, dist.meth) {
    
    # get area
    # raster density
    if (tp %in% c("density.overlap", "mean.density.overlap")){
      
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
    
    if (tp %in% c("mcp.overlap", "mean.mcp.overlap")){
      
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
    
    if (tp == "proportional.overlap") {
      
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
   
    # distance among points
    if (tp == "distance"){
      
      U <- rbind(W, Z)
      
      dists <- as.matrix(proxy::dist(U[ , dims], method = dist.meth,  convert_similarities = TRUE))
      
      dists <- dists[U[, group] == W[, group][1], U[, group] == Z[, group][1]]
      
      out <- matrix(mean(dists), nrow = 1)
    } 
    
    if (tp == "centroid.distance") {
      
      # both combined
      Y <- rbind(as.data.frame(W), as.data.frame(Z))
      
      frmla <- stats::as.formula(paste("cbind(", paste(dims, collapse = ","), ") ~ ", group))
      
      centroid.dist <- proxy::dist(stats::aggregate(formula = frmla, data = Y, FUN = mean)[, -1], method = dist.meth,  convert_similarities = TRUE)
      
      out <- matrix(centroid.dist, nrow = 1)  
    }
    
    return(out)
  }

  
  # get all combinations to get pairwise overlaps
  group_combs <- t(utils::combn(sort(unique(X[, group])), 2))
  
  # calculate all similarities
  similarities_l <- warbleR:::pblapply_wrblr_int(1:nrow(group_combs), pbar = pb,  cl = parallel, function(i, data = input_data, gc = group_combs, dims = dimensions, typ = type, pair.scale = pairwise.scale, dist.meth = distance.method) {
    
    if (type %in% c("density.overlap", "mean.density.overlap")){
      W <- data[[which(names(data) == gc[i, 1])]]
      Z <- data[[which(names(data) == gc[i, 2])]]
    } 
    
    if (type %in% c("distance", "mcp.overlap", "mean.mcp.overlap", "proportional.overlap", "centroid.distance")){
      W <- X_l[[which(names(X_l) == gc[i, 1])]]
      Z <- X_l[[which(names(X_l) == gc[i, 2])]]
      
      # pairwise scale
      if (pair.scale) {
        nrow_W <- nrow(W)
        
        # bind and scale
        U <- rbind(W, Z)
        U[, dimensions[1]] <- scale(U[, dimensions[1]])
        U[, dimensions[2]] <- scale(U[, dimensions[2]])
        
        # split back
        W <- U[1:nrow_W, ]
        Z <- U[(nrow_W + 1):nrow(U), ]
      }
    }
    
    suppressWarnings(similarities <- ovlp_fun(W, Z, typ, dims, dist.meth))

    # put in a data frame
    
    out_df <- if (type %in% c("mean.density.overlap", "mean.mcp.overlap", "centroid.distance", "proportional.overlap", "distance")) 
      data.frame(group.1 = gc[i, 1], group.2 = gc[i, 2], similarity = mean(similarities)) else
        data.frame(group.1 = gc[i, 1], group.2 = gc[i, 2], similarity.1in2 = similarities[2], similarity.2in1 = similarities[1]) 
    
    return(out_df)
  }) 
  
  # put in a data frame
  space_similarities <- do.call(rbind, similarities_l)

  # rename columns
  names(space_similarities) <- if (grepl("distance", type)) gsub("similarity", "distance", names(space_similarities)) else
    gsub("similarity", "overlap", names(space_similarities))

  return(space_similarities)
}

