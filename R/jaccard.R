#' jaccard: Estimates the jaccard index for comparing two binary maps
#'
#' @description Estimates the jaccard index for comparing two binary maps
#' @param m1 A binary raster A or an object of class setA returned by the function \code{\link[bam]{model2sparse}}.
#' @param m2 A binary raster A or an object of class setA returned by the function \code{\link[bam]{model2sparse}}.
#' @return Returns a numeric value representing the Jaccard index.
#' @details The Jaccard index is computed as follows
#' \deqn{J(A,B) = {{|A \cap B|}\over{|A \cup B|}} = {{|A \cap B|}\over{|A| + |B| - |A \cap B|}}.}
#' @export
#' @examples
#' m1_path <- system.file("extdata/conejos/Lepus_othus_cont.tif",
#'                        package = "bam")
#' m2_path <- system.file("extdata/conejos/Brachylagus_idahoensis_cont.tif",
#'                        package = "bam")
#' m1 <- raster::raster(m1_path) > 0.01
#' m2 <- raster::raster(m2_path) >0.01
#' jcc <- bam::jaccard(m1,m2)
#' print(jcc)
jaccard <- function(m1,m2){
  c1 <- class(m1)
  c2 <- class(m2)
  cl <- c("setA","RasterLayer")
  if(!c1 %in% cl){
    stop("m1 should be of class raster or setA")
  }
  if(!c2 %in% cl){
    stop("m2 should be of class raster or setA")
  }
  if(c1 == "setA"){
    m1 <- m1@sparse_model
  }
  if(c2 == "setA"){
    m2 <- m2@sparse_model
  }
  if(c1 == "RasterLayer"){
    m1 <- bam::model2sparse(m1)@sparse_model
  }
  if(c2 == "RasterLayer"){
    m2 <- bam::model2sparse(m2)@sparse_model
  }
  m3 <- m1 * m2
  A <- sum(m1)
  B <- sum(m2)
  int_AB <- sum(m3)
  jcc <-int_AB/(A+B- int_AB)

  return(jcc)
}
