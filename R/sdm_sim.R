#' sdm_sim: Simulate single species dispersal dynamics using the BAM framework.

#' @param set_A A setA object returned by the function \code{\link[bam]{model2sparse}}
#' @param set_M A setM object cointaining the adjacency matrix of the study area.
#' See \code{\link[bam]{adj_mat}}
#' @param initial_points A sparse vector returned by the function
#' \code{\link[bam]{occs2sparse}}
#' @param nsteps Number of steps to run the simulation
#' @importFrom magrittr %>%
#'
#' @examples
#' \dontrun{
#' model_path <- system.file("extdata/Lepus_californicus_cont.tif",
#'                           package = "bam")
#' model <- raster::raster(model_path)
#'
#' sparse_mod <- bam::model2sparse(model,threshold=0.05)
#' adj_mod <- bam::adj_mat(sparse_mod,ngbs=1)
#' occs_lep_cal <- data.frame(longitude = c(-122.08880,
#'                                          -98.89638),
#'                            latitude = c(37.43455,
#'                                         30.19919))
#'
#' occs_sparse <- bam::occs2sparse(modelsparse = sparse_mod,
#'                                 occs = occs_lep_cal)
#' smd_lep_cal <- bam::sdm_sim(set_A = sparse_mod,
#'                             set_M = adj_mod,
#'                             initial_points = occs_sparse,
#'                             nsteps = 10)
#'}
#' @export

sdm_sim <- function(set_A,set_M,initial_points,nsteps){
  if(!inherits(set_A,"setA")){
    stop("set_A should be of class setA")
  }
  if(!inherits(set_M,"setM")){
    stop("set_M should be of class setM")
  }
  AMA <- set_A@sparse_model %*% set_M@adj_matrix
  nsteps <- nsteps+1
  iter_vec <- 2:nsteps

  g0 <- initial_points
  sdm <- list(g0)
  for (i in iter_vec) {
    g0 <- AMA%*%g0
    g0[g0>1] <- 1
    #g0 <- (g0 - min(g0))/(max(g0)-min(g0))
    sdm[[i]] <- g0
  }
  bam_sim <- bam(sdm_sim =sdm,
                 bin_model=set_A@bin_model,
                 cellIDs=set_A@cellIDs,
                 sparse_model = set_A@sparse_model,
                 coordinates =set_A@coordinates,
                 adj_matrix = set_M@adj_matrix,
                 initial_points = initial_points,
                 sim_steps = nsteps-1)
  return(bam_sim)
}
