#' sim2Raster: Convert a BAM simulation object to RasterStack
#'
#' @description Convert a BAM simulation object to RasterStack.
#' @param sdm_simul A bam object. See \code{\link[bam]{sdm_sim}}
#' @param which_steps A numeric vector indicating the simulation steps that
#'  are going to be converted into raster layers.
#' @return A RasterStack of species' distribution at each simulation step
#' @export
#' @examples
#' \dontrun{
#' model_path <- system.file("extdata/Lepus_californicus_cont.tif",
#'                           package = "bam")
#' model <- raster::raster(model_path)

#' sparse_mod <- bam::model2sparse(model)
#' adj_mod <- bam::adj_mat(sparse_mod,nbgs=1)

#' occs_lep_cal <- data.frame(longitude = c(-115.10417,
#'                                          -104.90417),
#'                            latitude = c(29.61846,
#'                                         29.81846))

#' occs_sparse <- bam::occs2sparse(modelsparse = sparse_mod,
#'                                 occs = occs_lep_cal)
#' sdm_lep_cal <- bam::sdm_sim(set_A = sparse_mod,
#'                             set_M = adj_mod,
#'                             initial_points = occs_sparse,
#'                             nsteps = 100)
#' sdm_lep_cal_st <- bam::sim2Raster(sdm_simul = sdm_lep_cal,
#'                                   which_steps = seq(1,100,by=5))
#' }

sim2Raster <- function(sdm_simul,which_steps = NULL){

  if(!inherits(sdm_simul,"bam")){
    stop("sdm_simul should be of class bam")
  }
  if(is.null(which_steps))
    stepsvec <- 1:sdm_simul@sim_steps
  if(max(which_steps)>sdm_simul@sim_steps)
    stop("One or more steps were not simulated by the sdm_sim function")
  if(!all(which_steps%% 1 %in% 0))
    stop("which_steps should be integers")
  else{
    stepsvec <- which_steps
  }

  sim_stack <- stepsvec %>% purrr::map(function(x){
    grid_base <- sdm_simul@niche_model * 0
    grid_base[sdm_simul@cellIDs] <- sdm_simul@sdm_sim[[x]]
    return(grid_base)
  })
  sim_stack <- raster::stack(sim_stack)
  names(sim_stack) <- paste0("sim_",stepsvec)
  return(sim_stack)
}
