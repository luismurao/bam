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
#' sparse_mod <- bam::model2sparse(model, threshold=0.05)
#' }
model2sparse <- function(model, threshold=NULL){

  is_continous <- function(model){
    ff <- raster::freq(model,digits=2, useNA="no")
    if(nrow(ff)>3)
      return(TRUE)
    else
      return(FALSE)
  }

  source_model <- model

  if(is.numeric(threshold)){
    model <- model > threshold
    source_model <- source_model*model
  }

  is_cont <- is_continous(model = model)
  if(is_cont){
    base::stop("Please provide a suitability value to binarize model")
  }
  model_vals <- raster::getValues(model)*1
  in_calArea <- which(!is.na(model_vals))
  ncols <- nrows <- length(in_calArea)
  mod_sparse <- Matrix::sparseMatrix(i=1:length(in_calArea),
                                     j=1:length(in_calArea),
                                     x=model_vals[in_calArea],
                                     dims = c(nrows,ncols))*1
  mod_coords <- raster::coordinates(model)[in_calArea,]
  mod_atts <- setA(niche_model = model,
                   suit_threshold = ifelse(is.numeric(threshold),
                                           threshold,as.numeric(NA)),
                   cellIDs = in_calArea,
                   suit_values = if(is_continous(source_model)){
                     source_model[in_calArea]} else as.numeric(NA),
                   sparse_model =  mod_sparse,
                   coordinates= mod_coords)
  return(mod_atts)

}
