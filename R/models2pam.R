#' models2pam: Converts binary rasters to a PAM
#' @description Function to convert binary raster models to a Presence Absences Matrix.
#' @param mods_stack A raster stack containing binay models of each species in the community.
#' @param sparse Logical. If TRUE the PAM will be returned as a sparse matrix.
#' @return A PAM.
#' @export
#' @examples
#' \dontrun{
#' lagos_path <- system.file("extdata/conejos",
#'                           package = "bam")
#' enm_path <- list.files(lagos_path,
#'                        pattern = ".tif",
#'                        full.names = TRUE)

#' en_models <- raster::stack(enm_path) >0.01
#' pam <- bam::models2pam(en_models,sparse=FALSE)
#' }
models2pam <- function(mods_stack,sparse=TRUE){
  cmod <-class(mods_stack)
  if(!cmod %in% c("RasterStack","RasterBrick")){
    stop("'mods_stack' should be of class 'RasterStack' or 'RasterBrick'")
  }
  else{
    mods_stack <- raster::stack(mods_stack)
    nsps <- raster::nlayers(mods_stack)
    rvals <- mods_stack[[1]][]
    cellIDs <- which(!is.na(rvals))
    if(!sparse){
      pam0 <- 1:nsps %>% purrr::map_dfc(function(x){
        m1 <- mods_stack[[x]]
        m2 <- m1[]
        m3 <- m2[cellIDs]
        df1 <- data.frame(m3)
        df1 <-setNames(df1,names(m1))
        return(df1)
      })
      pam0 <- data.matrix(pam0)
      return(pam0)
    }
    else{

      pamL <- 1:nsps %>% purrr::map(function(x){
        m1 <- mods_stack[[x]]
        m2 <- m1[]
        m3 <- m2[cellIDs]
        unosIDs <- which(m3==1)
        msparse0 <- Matrix::sparseVector(x = rep(1,length(unosIDs)),
                                         i = unosIDs,
                                         length=length(msparse@cellIDs))
        return(msparse0)
      })
      pamL <- lapply(pamL, as, "sparseMatrix")

      pam0 <- do.call(Matrix::cBind, pamL)
      colnames(pam0) <- names(mods_stack)
      return(pam0)
    }

  }
}
