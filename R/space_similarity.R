#' @title Pairwise similarities of phenotype spaces
#'
#' @description \code{space_similarity} estimate pairwise similarities of phenotype spaces
#' @inheritParams template_params
#' @param method Character vector of length 1. Controls the method of (di)similarity metric to be compare the phenotypic sub-spaces of two groups at the time. Seven built-in metrics are available which quantify as pairwise sub-space overlap ('similarity') or pairwise distance between bi-dimensional sub-spaces ('dissimilarity'):
#' \itemize{
#'  \item \code{density.overlap}: proportion of the phenotypic sub-spaces area that overlap, taking into account the irregular densities of the sub-spaces. Two groups that share their higher density areas will be more similar than similar sub-spaces that only share their lower density areas. Two values are supplied as the proportion of the space of A that overlaps B is not necessarily the same as the proportion of B that overlaps A. Similarity metric (higher values means more similar). The minimum sample size (per group) must be 6 observations.
#'  \item \code{mean.density.overlap}: similar to 'density.overlap' but the two values are merged into a single pairwise mean overlap. Similarity metric (higher values means more similar). The minimum sample size (per group) must be 6 observations.
#'  \item \code{mcp.overlap}: proportion of the phenotypic sub-spaces area that overlap, in which areas are calculated as the minimum convex polygon of all observations for each su-space. Two values are supplied as the proportion of the space of A that overlaps B is not necessarily the same as the proportion of B that overlaps A. Similarity metric (higher values means more similar). The minimum sample size (per group) must be 5 observations.
#'  \item \code{mean.mcp.overlap}: similar to 'mcp.overlap' but the two values are merged into a single pairwise mean overlap. Similarity metric (higher values means more similar). The minimum sample size (per group) must be 5 observations.
#'  \item \code{proportional.overlap}: proportion of the joint area of both sub-spaces that overlaps (overlapped area / total area of both groups). Sub-space areas are calculated as the minimum convex polygon. Similarity metric (higher values means more similar). The minimum sample size (per group) must be 5 observations.
#'  \item \code{distance}: mean euclidean pairwise distance between all observations of the compared sub-spaces. Dissimilarity metric (higher values means less similar). The minimum sample size (per group) must be 1 observation.
#'  \item \code{centroid.distance}: euclidean distance between the centroid of the compared sub-spaces. Dissimilarity metric (higher values means less similar). The minimum sample size (per group) must be 1 observation.
#'  }
#'  In addition, machine learning classification models can also be used for quantify dissimilarity as a measured of how discriminable two groups are. These models can use more than two dimensions to represent phenotyypic spaces. The following classification models can be used: "AdaBag",           "avNNet", "bam", "C5.0", "C5.0Cost", "C5.0Rules", "C5.0Tree", "gam", "gamLoess", "glmnet",   "glmStepAIC", "kernelpls", "kknn", "lda", "lda2", "LogitBoost", "msaenet", "multinom", "nnet",     "null", "ownn", "parRF", "pcaNNet", "pls", "plsRglm", "pre", "qda", "randomGLM", "rf", "rFerns", "rocc", "rotationForest", "rotationForestCp", "RRF", "RRFglobal", "sda", "simpls", "slda", "smda", "snn", "sparseLDA", "svmLinear2", "svmLinearWeights", "treebag", "widekernelpls" and "wsrf". See \url{https://topepo.github.io/caret/train-models-by-tag.html} for details on each of these models. Additional arguments can be pased using \code{...}. Note that some machine learning methods can significantly affect computational efficiency (i.e. take a long time to compute). 
#' @param outliers Numeric vector of length 1. A value between 0 and 1 controlling the proportion of outlier observations to be excluded. Outliers are determined as those farthest away from the sub-space centroid. Ignored when using machine learning methods.
#' @param pairwise.scale Logical argument to control if pairwise phenotypic spaces are scaled (i.e. z-transformed) prior to similarity estimation. If so (\code{TRUE}) similarities are decoupled from the size of the global phenotypic space. Useful to compare similarities coming from different phenotypic spaces. Default is \code{FALSE}. Not available for 'density.overlap', 'mean.density.overlap' or any machine learning model.
#' @param distance.method Character vector of length 1 indicating the method to be used for measuring distances (hence only applicable when distances are calculated). Available distance measures are: "Euclidean" (default), "Manhattan", "supremum", "Canberra", "Wave", "divergence", "Bray", "Soergel", "Podani", "Chord", "Geodesic" and "Whittaker". If a similarity measure is used similarities are converted to distances.
#' @param seed Integer number containing the random number generator (RNG) state for random number generation in order to make results from the machine learning stochastic methods replicable.
#' @param ... Additional arguments to be passed to \code{\link[caret]{train}}.
#' @return A data frame containing the similarity metric for each pair of groups. If the similarity metric is not symmetric (e.g. the proportional area of A that overlaps B is not necessarily the same as the area of B that overlaps A, see \code{\link{space_similarity}}) separated columns are supplied for the two comparisons.
#' @export
#' @name space_similarity
#' @details The function quantifies pairwise similarity between phenotypic sub-spaces. The built-in methods quantify similarity as the overlap (similarity, or machine learning based discriminability) or distance (dissimilarity) between group. Machine learning methods implemented in the caret package function \code{\link[caret]{train}} are available to assess the similarity of spaces as the proportion of observations that are incorrectly classified. In this case group overlaps are the class-wise errors (if available) while the mean overlap is calculated as \code{1- model accuracy}.      
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
#' #' # get symmetric triangular matrix
#' rectangular_to_triangular(prop_overlaps)
#'
#' # get minimum convex polygon overlap for each group (non-symmetric)
#' mcp_overlaps <- space_similarity(
#'  formula = group ~ dimension_1 + dimension_2,
#'  data = example_space,
#'  method = "mcp.overlap")
#'
#' # convert to non-symmetric triangular matrix
#' rectangular_to_triangular(mcp_overlaps, symmetric = FALSE)
#'
#' # check available distance measures
#' summary(proxy::pr_DB)
#'
#' # get eculidean distances (default)
#' area_dist <- space_similarity(
#'  formula = group ~ dimension_1 + dimension_2,
#'  data = example_space,
#'  method = "distance",
#'  distance.method = "Euclidean")
#'
#' # get Canberra distances
#' area_dist <- space_similarity(
#'  formula = group ~ dimension_1 + dimension_2,
#'  data = example_space,
#'  method = "distance",
#'  distance.method = "Canberra")
#' 
#' ## using machine learning classification methods
#' 
#' # check if caret package and needed dependencies are available
#'  rlang::check_installed("caret")
#'  rlang::check_installed("randomForest")
#'  
#' # random forest 3 dimension data, using 5 repeats and repeated CV resampling
#' # extract data subset
#' sub_data <- example_space[example_space$group %in% c("G1", "G2", "G3"), ]
#' 
#' # set method parameters
#' ctrl <- caret::trainControl(method = "repeatedcv", repeats = 5)
#'
#' # get similarities ("overlap")
#' space_similarity(
#'  formula = group ~ dimension_1 + dimension_2 + dimension_3,
#'  data = sub_data,
#'  method = "rf",
#'  trControl = ctrl, 
#'  tuneLength = 4,
#'  seed = 123
#' )
#'
#' # Single C5.0 Tree using boot resampling
#' ctrl <- caret::trainControl(method = "boot")
#'
#' space_similarity(
#'  formula = group ~ dimension_1 + dimension_2,
#'  data = sub_data,
#'  method = "C5.0Tree",
#'  trControl = ctrl,
#'  tuneLength =  3
#')
#' }
#' @seealso \code{\link{rarefact_space_similarity}}, \code{\link{space_size_difference}}
#' @author Marcelo Araya-Salas \email{marcelo.araya@@ucr.ac.cr})
#'
#' @references {
#' Araya-Salas, M, & K. Odom. 2022, PhenotypeSpace: an R package to quantify and compare phenotypic trait spaces R package version 0.1.0.
#' }

space_similarity <-
  function(formula,
           data,
           cores = 1,
           method = "mcp.overlap",
           pb = TRUE,
           outliers = 0.95,
           pairwise.scale = FALSE,
           distance.method = "Euclidean",
           seed = NULL,
           ...) {
    
    # supported machine learning methods
    ml_methods <- c("AdaBag",
                    "avNNet",
                    "bam",
                    "C5.0",
                    "C5.0Cost",
                    "C5.0Rules",
                    "C5.0Tree",
                    "gam",
                    "gamLoess",
                    "glmnet",
                    "glmStepAIC",
                    "kernelpls",
                    "kknn",
                    "lda",
                    "lda2",
                    "LogitBoost",
                    "msaenet",
                    "multinom",
                    "nnet",
                    "null",
                    "ownn",
                    "parRF",
                    "pcaNNet",
                    "pls",
                    "plsRglm",
                    "pre",
                    "qda",
                    "randomGLM",
                    "rf",
                    "rFerns",
                    "rocc",
                    "rotationForest",
                    "rotationForestCp",
                    "RRF",
                    "RRFglobal",
                    "sda",
                    "simpls",
                    "slda",
                    "smda",
                    "snn",
                    "sparseLDA",
                    "svmLinear2",
                    "svmLinearWeights",
                    "treebag",
                    "widekernelpls", 
                    "wsrf")
    
    if (!method %in% c(
      # phenotype space methods
      "density.overlap",
      "mean.density.overlap",
      "mcp.overlap",
      "mean.mcp.overlap",
      "proportional.overlap",
      "distance",
      "centroid.distance",
      ml_methods
    ))
    stop2("Unsupported 'method' declared")
    
    if (!distance.method %in% c("Euclidean", "Manhattan", "supremum", "Canberra", "Wave", "divergence", "Bray", "Soergel", "Podani", "Chord", "Geodesic", "Whittaker"))
      stop2("Unsupported 'distance.method' declared")
    
    # check if model package is installed and offer user to install it if not 
    if (method %in% ml_methods) {
      rlang::check_installed("caret")
      
      # install dependenciesfor machine learning method
      caret::checkInstall(caret::getModelInfo(method, regex = FALSE)[[1]]$library)
    }
    
    # get term names from formula 
    form_terms <- terms(formula, data = data)
    dimensions <- attr(form_terms, "term.labels")
    group <- as.character(form_terms[[2]])

    # center data
    # if (method == "centroid.distance") {
    #
    #   data[, dimensions[1]] <- scale(data[, dimensions[1]])
    #   data[, dimensions[2]] <- scale(data[, dimensions[2]])
    #
    #   } else {
    # # force dimensions to range 0-1
    # data[, dimensions[1]] <- (data[, dimensions[1]] - min( data[, dimensions[1]])) / max((data[, dimensions[1]] - min( data[, dimensions[1]])))
    # data[, dimensions[2]] <- (data[, dimensions[2]] - min( data[, dimensions[2]])) / max((data[, dimensions[2]] - min( data[, dimensions[2]])))
    #
    #   }
    
    # split in a list of vectors
    split_data_list <- split(x = data, f = data[, group])
    
    # stop if too small sample sizes
    if (min(sapply(split_data_list, nrow)) < 2)
      stop2(
        "There is at least one group with less than 2 observations which is the minimum needed for similarity estimation"
      )
    
    # get densities
    if (method %in% c("density.overlap", "mean.density.overlap")) {
      total_coors <-
        spatstat.geom::as.ppp(as.matrix(data[, dimensions]), c(range(data[, dimensions[1]]), range(data[, dimensions[2]])))
      total_space <-
        raster::raster(spatstat.explore::density.ppp(total_coors))
      
      input_data <-
        lapply(split_data_list, function(Y, dims = dimensions) {
          coors <-
            spatstat.geom::as.ppp(as.matrix(Y[, dimensions]), c(range(Y[, dimensions[1]]), range(Y[, dimensions[2]])))
          
          raster_dens <-
            raster::raster(spatstat.explore::density.ppp(coors))
          
          raster_dens <-
            raster::crop(raster::extend(raster_dens, total_space), total_space)
          
          # convert to range 0-1
          raster::values(raster_dens) <-
            (raster::values(raster_dens) - min(raster::values(raster_dens), na.rm = TRUE))
          raster::values(raster_dens) <-
            raster::values(raster_dens) / max(raster::values(raster_dens), na.rm = TRUE)
          
          # keep 95% interval
          raster::values(raster_dens)[raster::values(raster_dens) < 1 - outliers] <-
            NA
          
          # same resolution and extend as total space
          raster_dens <- raster::resample(raster_dens, total_space)
          raster::extent(raster_dens)  <- c(0, 1, 0, 1)
          
          return(raster_dens)
        })
      
      names(input_data) <- names(split_data_list)
      
    } else {
      input_data <- split_data_list
    }
    
    # function to calculate areas
    ovlp_fun <- function(W, Z, WZ, tp, dims, dist.meth, seed) {
      # get area
      # raster density
      if (tp %in% c("density.overlap", "mean.density.overlap")) {
        # convert NAs to 0
        raster::values(W)[is.na(raster::values(W))] <- 0
        raster::values(Z)[is.na(raster::values(Z))] <- 0
        
        # 1 vs 2
        wgt_ovlp.WvZ <-
          sum((raster::values(W) * raster::values(Z)) / raster::values(Z),
              na.rm = TRUE) / sum(raster::values(Z), na.rm = TRUE)
        
        # 2 vs 1
        wgt_ovlp.ZvW <-
          sum((raster::values(Z) * raster::values(W)) / raster::values(W),
              na.rm = TRUE) / sum(raster::values(W), na.rm = TRUE)
        
        # convert to 1 if higher than 1
        if (wgt_ovlp.WvZ > 1)
          wgt_ovlp.WvZ <- 1
        if (wgt_ovlp.ZvW > 1)
          wgt_ovlp.ZvW <- 1
        
        
        # put results in am matrix
        out <- matrix(c(wgt_ovlp.WvZ, wgt_ovlp.ZvW, mean(c(wgt_ovlp.WvZ, wgt_ovlp.ZvW))), nrow = 1)
      }
      
      if (tp %in% c("mcp.overlap", "mean.mcp.overlap")) {
        sp::coordinates(W) <-
          stats::as.formula(paste("~", paste(dims, collapse = "+")))
        sp::coordinates(Z) <-
          stats::as.formula(paste("~", paste(dims, collapse = "+")))
        
        # Create the MCPs (Minimum Convex Polygons) for W and Z
        mcp_W <- sf::st_as_sf(adehabitatHR::mcp(W))
        mcp_Z <- sf::st_as_sf(adehabitatHR::mcp(Z))
        
        # Perform the intersection
        intrsct <- sf::st_intersection(mcp_W, mcp_Z)

        if (nrow(intrsct) == 0) {
          ovlp1in2 <- ovlp2in1 <- 0
        } else {
          # Calculate the area of the intersection (in square meters by default)
          intrsctspace_size <- sf::st_area(intrsct)
          
          ovlp1in2 <-
            intrsctspace_size / sf::st_area(mcp_W)
          ovlp2in1 <-
            intrsctspace_size / sf::st_area(mcp_Z)
          
        } 
        
        out <- matrix(c(ovlp2in1, ovlp1in2, mean(c(ovlp2in1, ovlp1in2))), nrow = 1)
      }
      
      if (tp == "proportional.overlap") {
        sp::coordinates(W) <-
          stats::as.formula(paste("~", paste(dims, collapse = "+")))
        sp::coordinates(Z) <-
          stats::as.formula(paste("~", paste(dims, collapse = "+")))
        
        # both combined
        Y <- rbind(W, Z)
        
        # get intersect
        mcp_W <- sf::st_as_sf(adehabitatHR::mcp(W, percent = outliers * 100))
        mcp_Z <- sf::st_as_sf(adehabitatHR::mcp(Z, percent = outliers * 100))

        # Perform the intersection
        intrsct <- sf::st_intersection(mcp_W, mcp_Z)
        
        if (nrow(intrsct) == 0) {
          ovlp1in2 <- ovlp2in1 <- 0
        } else {
          # intrsctspace_size <- raster::area(intrsct)
          intrsctspace_size <- sf::st_area(intrsct)
          
          # totalspace_size <-
            # raster::area(adehabitatHR::mcp(Y, percent = outliers * 100))
            
            mcp_Y <- sf::st_as_sf(adehabitatHR::mcp(Y), percent = outliers * 100)  
            totalspace_size <- sf::st_area(mcp_Y)
          ovlp1in2 <-
            ovlp2in1 <- intrsctspace_size / totalspace_size
          
        } 
        
        out <- matrix(c(ovlp2in1, ovlp1in2, mean(c(ovlp2in1, ovlp1in2))), nrow = 1)
      }
      
      # distance among points
      if (tp == "distance") {
        U <- rbind(W, Z)
        
        dists <-
          as.matrix(proxy::dist(U[, dims], method = dist.meth,  convert_similarities = TRUE))
        
        dists <-
          dists[U[, group] == W[, group][1], U[, group] == Z[, group][1]]
        
        out <- matrix(mean(dists), nrow = 1)
      }
      
      if (tp == "centroid.distance") {
        # both combined
        Y <- rbind(as.data.frame(W), as.data.frame(Z))
        
        frmla <-
          stats::as.formula(paste("cbind(", paste(dims, collapse = ","), ") ~ ", group))
        
        centroid.dist <-
          proxy::dist(
            stats::aggregate(
              x = frmla,
              data = Y,
              FUN = mean
            )[, -1],
            method = dist.meth,
            convert_similarities = TRUE
          )
        
        out <- matrix(centroid.dist, nrow = 1)
      }
      
      # if any machine learning model
      if (tp %in% ml_methods) {
        if (!is.null(seed))
          set.seed(seed)
        
        capture.output(train_caret <-
          caret::train(form = formula,
                data = WZ,
                method = tp,
                verbose = FALSE,
                ...
                ))
        
        out <- matrix(1 - train_caret$results$Accuracy[nrow(train_caret$results)], nrow = 1)
        
        if (!is.null(train_caret$finalModel$confusion)){
        conf_matrix <- train_caret$finalModel$confusion
        out <- matrix(c(conf_matrix[rownames(conf_matrix) == W[1, group], 3], conf_matrix[rownames(conf_matrix) == Z[1, group], 3], 1 - train_caret$results$Accuracy[nrow(train_caret$results)]), nrow = 1)
        }
      }
      return(out)
    }
    
    # get all combinations to get pairwise overlaps
    group_combs <- t(utils::combn(sort(unique(data[, group])), 2))
    
    # calculate all similarities
    similarities_l <-
      pblapply_phtpspc_int(1:nrow(group_combs), pbar = pb,  cl = cores, function(i,
    data = input_data,
    gc = group_combs,
    dims = dimensions,
    typ = method,
    pair.scale = pairwise.scale,
    dist.meth = distance.method,
    sed = seed
    ) {
        if (typ %in% c("density.overlap", "mean.density.overlap")) {
          W <- data[[which(names(data) == gc[i, 1])]]
          Z <- data[[which(names(data) == gc[i, 2])]]
        }
        
        if (typ %in% c(
          "distance",
          "mcp.overlap",
          "mean.mcp.overlap",
          "proportional.overlap",
          "centroid.distance",
          ml_methods)) {
          W <- split_data_list[[which(names(split_data_list) == gc[i, 1])]]
          Z <-
            split_data_list[[which(names(split_data_list) == gc[i, 2])]]
          
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
        
          if (typ %in% ml_methods) {
            WZ <- rbind(W, Z)
            WZ <- droplevels(WZ)
          } else {
              WZ <- NULL
            }
        }
        
        suppressWarnings(similarities <-
                           ovlp_fun(W, Z, WZ, typ, dims, dist.meth, seed = sed))
        
        # put in a data frame
        
        out_df <-
          if (typ %in% c(
            "mean.density.overlap",
            "mean.mcp.overlap",
            "centroid.distance",
            "proportional.overlap",
            "distance"
          ) | ncol(similarities) == 1) {
            data.frame(
              group.1 = gc[i, 1],
              group.2 = gc[i, 2],
              similarity = mean(similarities)
            )
          } else
            data.frame(
              group.1 = gc[i, 1],
              group.2 = gc[i, 2],
              similarity.1in2 = similarities[2],
              similarity.2in1 = similarities[1],
              mean.similarity = similarities[3]
            )
        
        return(out_df)
      })
    
    # put in a data frame
    space_similarities <- do.call(rbind, similarities_l)
    
    # rename columns
    names(space_similarities) <-
      if (grepl("distance", method)) {
        gsub("similarity", "distance", names(space_similarities))
      } else {
        gsub("similarity", "overlap", names(space_similarities))
      }
    
    return(space_similarities)
  }
