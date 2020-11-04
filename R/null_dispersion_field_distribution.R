#' null_dispersion_field_distribution: Null distribution of the dispersion field
#' @description null_dispersion_field_distribution estimates a random distribuion
#' of the dispersion field values.
#' @param pam A Presence-Absence-Matrix of matrix class or sparse matrix.
#' @param n_iter Number of iterations to obtain the distribution.
#' @param parallel If TRUE the computations will be performed in parallel.
#' @param n_cores Number of cores for the parallel computation.
#' @importFrom Rdpack reprompt
#' @importFrom Rcpp sourceCpp
#' @useDynLib bam
#' @return A data matrix of size nrow(pam) X n_iter with dispersion field values.
#' @details
#'  Estimates a random distribution of the dispersion field values. To obtein
#'          random values it uses the function code{\link[bam]{permute_pam}}
#'          at each step of the iterations.
#'
#' @references
#' \insertRef{Soberon2015}{bam}
#' @examples
#' set.seed(111)
#' pam <- matrix(rbinom(100,1,0.3),nrow = 10,ncol = 10)
#' dfield_rand <- bam::null_dispersion_field_distribution(pam,n_iter=10,
#'                                                        parallel=FALSE,
#'                                                        n_cores = 2)
#' @export

null_dispersion_field_distribution <- function(pam,n_iter=10,parallel=TRUE,n_cores=2){

  if(is.data.frame(pam))
    pam <- data.matrix(pam)
  if(!is.matrix(pam))
    stop("m should be a matrix or a data.frame")
  sniter <- 1:n_iter
  if(parallel){
    plan(tweak(multisession, workers = n_cores))
  } else{
    plan(sequential)
  }

  nms <-paste0("dfrand_",sniter)
  distfield_rand <- sniter %>% furrr::future_map_dfc(function(x){
    ppam <- bam::permute_pam(m=pam,as_sparse=T)
    distfield <-bam::pam2bioindex(pam=ppam,
                                  biodiv_index = "dispersion_field",
                                  as_sparse = F)

    y=data.frame(dfield =distfield@dispersion_field)
    return(y)
  },.progress = T)
  distfield_rand <- data.matrix(distfield_rand)
  colnames(distfield_rand) <- nms
  return(distfield_rand)
}
