#' model2sparse: Converts a niche model into a diagonal sparse matrix
#' @param model A raster object representing the geographic projection of a niche model.
#' @param threshold A threshold to convert a continuous model into a binary model.
#' @import Matrix
#' @return A diagonal sparse matrix representing the geographic projection of a niche model.
#' @export
#' @examples
#' \dontrun{
#' model_path <- system.file("extdata/Lepus_californicus_cont.tif",
#'                           package = "bam")
#' model <- raster::raster(model_path)
#'
#' sparse_mod <- bam::model2sparse(model,threshold=0.05)
#' }
model2sparse <- function(model, threshold=NULL){
  if(is.numeric(threshold)){
    model <- model > threshold
  }
  model_vals <- raster::getValues(model)
  in_niche <- which(!is.na(model_vals))
  cell_ids <- stats::complete.cases(model_vals)
  out_niche <- which(cell_ids[-in_niche])
  all_area <- which(cell_ids)

  ncols <- nrows <- length(all_area)
  mod_sparse <- Matrix::sparseMatrix(i=match(in_niche,all_area),
                                     j=match(in_niche,all_area),
                                     x=model_vals[in_niche],
                                     dims = c(nrows,ncols))*1
  mod_coords <- raster::coordinates(model)[all_area,]
  mod_atts <- setA(bin_model = model,
                   cellIDs = all_area,
                   sparse_model =  mod_sparse,
                   coordinates= mod_coords)
  return(mod_atts)

}

