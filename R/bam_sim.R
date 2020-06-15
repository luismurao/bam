#' bam_sim: Simulate dispersal dynamics using the set B of the BAM framework.

#' @param sp1 Niche model of the focal species (the one that disperses).
#' @param sp2 Niche model of the species with whom sp1 interacts (currently no dispersal dynamics for this species).
#' @param set_M A setM object cointaining the adjacency matrix for sp1.
#' See \code{\link[bam]{adj_mat}}
#' @param periods_negative  Time periods that sps2 takes to develop defense mechanisms (i.e. toxic).
#' @param periods_positive This the time that sp2 takes to become non-toxic
#' @param initial_points A sparse vector returned by the function
#' \code{\link[bam]{occs2sparse}}
#' @param nsteps Number of steps to run the simulation
#' @param progress_bar Show progress bar
#' @importFrom magrittr %>%
#'
#' @examples
#' \dontrun{
#' ura <- raster::raster(system.file("extdata/urania_omph/urania_guanahacabibes.tif",
#'                                   package = "bam"))
#' omp <- raster::raster(system.file("extdata/urania_omph/omphalea_guanahacabibes.tif",
#'                                   package = "bam"))
#' msparse <- bam::model2sparse(ura)
#' init_coordsdf <- data.frame(x=-84.38751, y= 22.02932)
#' initial_points <- bam::occs2sparse(modelsparse = msparse,init_coordsdf)
#' set_M <- bam::adj_mat(modelsparse = msparse,ngbs = 1)
#' ura_sim <- bam_sim(sp1=ura, sp2=omp, set_M=set_M,
#'                    initial_points=initial_points,
#'                    periods_negative=5,periods_positive=1,
#'                    nsteps=40)
#' # Animation example
#' anp <-"C:/Users/l916o895/Dropbox/TeoriadeBAM/articulo_bam/ura_omp_sim.gif"
#' new_sim <- bam::sim2Animation(sdm_simul = ura_sim,
#'                               which_steps = 1:ura_sim@sim_steps,
#'                               fmt = "GIF",
#'                               filename = anp)
#'}
#' @export
#'
#'

bam_sim <- function(sp1,sp2,set_M,initial_points,
                    periods_negative,periods_positive,
                    nsteps,progress_bar=TRUE){
  if(!(class(sp1)  == "RasterLayer" && class(sp2)  == "RasterLayer")){
    stop("sp1 and sp2 should be of raster class")
  }
  if(!inherits(set_M,"setM")){
    stop("set_M should be of class setM")
  }
  sp1_sp2 <- sp1*sp2
  A <- bam::model2sparse(sp1_sp2)
  bin_model <- A@sparse_model

  AMA <- A@sparse_model %*% set_M@adj_matrix %*% A@sparse_model
  g0 <- initial_points
  nonzerosL <- list()
  to_off_vec <- seq(0,nsteps,by=periods_negative)[-1]
  to_on_vec <- to_off_vec + periods_positive
  index_off <- 1
  index_on <- 1
  sdm <- list(g0)
  if(progress_bar){
    pb <- utils::txtProgressBar(min = 0,
                                max = nsteps,
                                style = 3)
  }

  for(i in 1:nsteps){
    g0 <- AMA%*%g0
    pos <- .nonzero(g0)[,1]
    nonzerosL[[i]] <- pos

    if(i %in% to_off_vec){
      to_off <-  index_off:i
      nonzeros <- unlist(nonzerosL[to_off])
      tb_nonzeros <- table(nonzeros)
      to_converIDs <- which(tb_nonzeros>=periods_negative)
      cids <- as.numeric(names(to_converIDs))
      diag(bin_model[cids,cids]) <- 0
      AMA <-bin_model %*% set_M@adj_matrix
      index_off <- i+1
    }

    if(i %in%  to_on_vec && exists("cids")){
      to_on <- index_on:(index_off-1)
      nonzeros <- unlist(nonzerosL[to_on])
      tb_nonzeros <- table(nonzeros)
      to_converIDs <- which(tb_nonzeros>=periods_positive)
      cids <- as.numeric(names(to_converIDs))
      Matrix::diag(bin_model[cids ,cids ]) <- 1
      AMA <-bin_model %*% set_M@adj_matrix
    }
    g0[g0>1] <- 1
    sdm[[i+1]] <- g0
    #cat("\nDoing step ",i, "of ",nsteps,"\n")
    if(progress_bar){
      utils::setTxtProgressBar(pb, i)
    }

  }

  bamsim <- bam(sdm_sim =sdm,
                bin_model=sp1_sp2,
                cellIDs=A@cellIDs,
                sparse_model = A@sparse_model,
                coordinates =A@coordinates,
                adj_matrix = set_M@adj_matrix,
                initial_points = initial_points,
                sim_steps = nsteps)
  return(bamsim)

}


#' Helper function to compute the elements in g0
#' that have no zero values.The function is taken from the
#' Ringo package
#' @param x A matrix of class "dgCMatrix"

.nonzero <- function(x){
  stopifnot(inherits(x, "dgCMatrix"))
  if (all(x@p == 0))
    return(matrix(0, nrow=0, ncol=2,
                  dimnames=list(character(0), c("row","col"))))
  res <- cbind(x@i+1, rep(seq(dim(x)[2]), diff(x@p)))
  colnames(res) <- c("row", "col")
  res <- res[x@x != 0, , drop = FALSE]
  return(res)
}


