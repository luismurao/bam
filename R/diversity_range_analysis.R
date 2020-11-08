#' range_diversity_analysis: diversity analysis
#' @description diversity_range_analysis biodiversity indices related to
#' diversity-range plots
#' @param pam A Presence-Absence-Matrix of matrix class or sparse matrix.
#' @param xy_mat A two dimensional matrix with longitude and latitude data.
#' @param lower_interval Lower interval.
#' @param upper_interval Upper interval.
#' @param raster_templete Araster templete.
#' @param niter Number of iterations to obtain the distribution.
#' @param return_null_dfield If TRUE the null distribution of dispersal field
#' will be returned.
#' @param parallel If TRUE the computations will be performed in parallel.
#' @param n_cores Number of cores for the parallel computation.
#' @importFrom future plan tweak multisession sequential
#' @importFrom Rcpp sourceCpp
#' @useDynLib bam
#' @export
#' @examples
#' \dontrun{
#' set.seed(111)
#' pam <- matrix(rbinom(10000,1,0.5),nrow = 100,ncol = 1000)
#' rdivan <- diversity_range_analysis(pam=pam,parallel = FALSE,
#'                                    return_null_dfield=TRUE)
#' bam::plot(rdivan,plot_type="diversity_range")
#' }

diversity_range_analysis <- function(pam,xy_mat=NULL,lower_interval=0.05,
                                     upper_interval=0.95, raster_templete=NULL,
                                     niter=100,return_null_dfield=FALSE,parallel=TRUE,
                                     n_cores=4){


  results <- methods::new(Class = "diversity_range")

  distfield_rand <- null_dispersion_field_distribution(pam = pam,
                                                       n_iter=niter,
                                                       parallel=parallel,
                                                       n_cores =n_cores)

  if(return_null_dfield){
    results@null_dispersion_field_dist <- distfield_rand
  }


  bioind <- bam::pam2bioindex(pam=pam,
                             biodiv_index = c("alpha",
                                              "dispersion_field"),
                             as_sparse = F)

  results@alpha <- bioind@alpha
  results@omega <- bioind@omega
  results@nsps <- ncol(pam)
  results@nsites <- nrow(pam)
  results@n_iterations <- niter
  results@dispersion_field <- bioind@dispersion_field
  disfield_cat <-  null_dispersion_field_cat(dfield = bioind@dispersion_field,
                                             dfield_rand = distfield_rand,
                                             lower_interval=lower_interval,
                                             upper_interval=upper_interval)


  #cols <- c("#000000","#F6BDC0",
  #          "#F07470","#BBDFFA",
   #         "#DC1C13","#6987D5",
  #          "#1727AE")
  cols <- c("#000000","#F6BDC0",
            "#F1A13A","#BBDFFA",
            "#DC1C13","#6987D5",
            "#1727AE")
  names(cols) <- c("Random","HE/LR",
                   "HE/IR","LE/LR",
                   "HE/HR","LE/IR",
                   "LE/HR")

  richness_cat <-  bioind@alpha
  q_rich <-  stats::quantile(richness_cat)
  q_rich_low <- q_rich[2]
  q_rich_int <- q_rich[4]
  q_rich_hig <- q_rich[5]

  q_rich_low_ids <- which(richness_cat[,1] <  q_rich_low)
  q_rich_int_ids <- which(richness_cat[,1] >= q_rich_low &
                            richness_cat[,1] < q_rich_int)
  q_rich_hig_ids <- which(richness_cat[,1] >= q_rich_int &
                            richness_cat[,1] <= q_rich_hig)

  richness_cat[q_rich_low_ids,1] <- 1
  richness_cat[q_rich_int_ids,1] <- 2
  richness_cat[q_rich_hig_ids,1] <- 4
  dfalpha <- disfield_cat*richness_cat

  dfalpha <- disfield_cat*richness_cat
  range_div_cols <-  dfalpha
  #alpha_rasterC <-   range_div_cols
  range_div_cols <- as.factor(dfalpha)
  levels(range_div_cols) <- cols
  #range_div_cols <- as.character(range_div_cols)
  results@diversity_range_colors <- as.character(range_div_cols)
  if(is.matrix(xy_mat) || is.data.frame(xy_mat)){


    results@xy_coordinates <- data.matrix(xy_mat)

    if(class(raster_templete) == "RasterLayer"){

      r1 <- raster_templete*0
      alpha_raster <- r1
      names(alpha_raster) <- "alpha"
      dispersion_field_raster <- r1
      names(dispersion_field_raster) <- "dispersion_field"
      diversity_range_raster <- r1
      names(diversity_range_raster) <- "diversity_range"
      cellIDs <- raster::cellFromXY(r1,xy_mat[,1:2])
      alpha_raster[cellIDs]<- bioind@alpha
      results@alpha_raster <- alpha_raster
      dispersion_field_raster[cellIDs] <- bioind@dispersion_field
      diversity_range_raster[cellIDs] <- range_div_cols
      results@dispersion_field_raster <- dispersion_field_raster
      results@diversity_range_raster<- diversity_range_raster
    }

  }

  return(results)

}
