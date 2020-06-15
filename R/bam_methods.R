
#' Show information in setA class \pkg{bam}.
#' @importFrom methods new
#' @param object An object of class setA
#' @rdname show
#' @export

methods::setMethod(f = "show",
                   signature = "setA",
                   function(object) {
                     slotsin <- methods::slotNames(object)

                     cat("Set A of the BAM digram it contains",
                         length(slotsin),"slots \n\n")
                     npixs <- length(object@cellIDs)
                     if("bin_model" %in% slotsin){
                       cat("@bin_model: a niche model:\n\n")
                       print(object@bin_model)
                     }
                     if("cellIDs" %in% slotsin){
                       cat("@cellIDs: ids of the cells that have values",
                           paste0("(",npixs," pixels)"),"\n\n")
                     }
                     if("sparse_model" %in% slotsin){
                       cat("@sparse_model:",
                           "A sparse square matrix of ",
                           npixs,"x",npixs,
                           "entries (showing only 4 entries) \n\n")
                       print(head(object@sparse_model[1:2,1:2]))
                     }
                    # if("occs_sparse" %in% slotsin){
                    #   cat("Number of occurences:",
                    #       (object@n_occs))
                    # }

                   })

#' Show information in csd class \pkg{bam}.
#' @importFrom methods new
#' @param object An object of class setA
#' @rdname show
#' @export

methods::setMethod(f = "show",
                   signature = "csd",
                   function(object) {
                     slotsin <- methods::slotNames(object)
                     cat("Object of class csd it contains",
                         length(slotsin),"slots: \n\n")
                     cat("@connections: Geographic clusters data.frame \n\n")
                     print(head(object@connections,4))
                     cat("@interactive_map: A leaflet map showing the geographic clusters\n\n")
                     cat("@raster_map: A raster map of the clusters")
                   })


#' Show information in pam class \pkg{bam}.
#' @importFrom methods new
#' @param object An object of class setA
#' @rdname show
#' @export

methods::setMethod(f = "show",
                   signature = "pam",
                   function(object) {
                     slotsin <- methods::slotNames(object)
                     cat("Object of class pam it contains",
                         length(slotsin),"slots: \n\n")
                     cat("@pams: A list of sparse PAM matrices\n\n")
                     print(head(object@pams[[1]]))
                     cat("\n@sp_names: the species in the pam with ",
                         length(object@sp_names), "species:\n",
                         paste0(object@sp_names,collapse = ", "))
                     cat("\n\n@which_steps: time steps",object@which_steps,"\n\n")
                     cat("@grid:", "raster grid of the studied area\n\n" )
                     print(object@grid)
                     cat("\n@cellIDs:", "site ids regarding the raster grid\n\n" )

                   })

#' Show information in setA class \pkg{bam}.
#' @importFrom methods new
#' @param object An object of class setA
#' @rdname show
#' @export

methods::setMethod(f = "show",
                   signature = "setM",
                   function(object) {
                     slotsin <- methods::slotNames(object)

                     cat("Set M of the BAM digram it contains",
                         length(slotsin),"slots \n\n")


                     cat("@coordinates: A matrix with longitude and latitude",
                         "values of each cell of the raster area\n\n")
                     print(head(object@coordinates))

                     cat("@initial_points: A list of inital coordinates where the",
                          "invasion process starts\n\n")
                     if(length(object@initial_points)>0L)
                       print(head(object@initial_points))

                     cat("@eigen_val: Eigen values of the connectivity matrix M\n\n")
                     if(length(object@eigen_val)>0L)
                       print(head(object@eigen_val))

                     cat("@eigen_vec: Eigen vector of the connectivity matrix M\n\n")
                     if(length(object@eigen_vec)>0L)
                       print(head(object@eigen_vec))
                   })


#' Predict method of the package \pkg{bam}.
#' @aliases predict bam-method predict
#' @description predicts species' distribution under suitability changes
#' @param object a of class bam.
#' @param niche_layers A raster or RasterStack with the niche models for
#' each time period
#' @param nbgs_vec A vector with the number of neighbors for the adjacency matrices
#' @param nsteps_vec Number of simulation steps for each time period.
#' @param animate Logical. If TRUE a dispersal animation on climate change scenarios will be created
#' @param period_names Character vector with the names of periods that will be animated. Default NULL.
#' @param ani.width Animation width unit in px
#' @param ani.height Animation height unit in px
#' @param ani.res Animation resolution unit in px
#' @param fmt Animation format. Posible values are GIF and HTML
#' @param filename File name.
#' @param png_keyword A keyword name for the png images generated by the function
#' @export
#' @rdname predict
#' @examples
#' \dontrun{
#' # Load R packages
#' library(bam)
#' library(raster)
#' # rm(list = ls())
#' # Read raster model for Lepus californicus
#' model_path <- system.file("extdata/Lepus_californicus_cont.tif",
#'                           package = "bam")
#' model <- raster::raster(model_path) >0.1
#' # Convert model to sparse
#' sparse_mod <- bam::model2sparse(model = model)
#' # Compute adjacency matrix
#' adj_mod <- bam::adj_mat(sparse_mod,ngbs=1)
#'
#' # Initial points to start dispersal process
#'
#' occs_lep_cal <- data.frame(longitude = c(-115.10417,
#'                                          -104.90417),
#'                            latitude = c(29.61846,
#'                                         29.81846))
#' # Convert to sparse the initial points
#' occs_sparse <- bam::occs2sparse(modelsparse = sparse_mod,
#'                                 occs = occs_lep_cal)
#'
#' # Run the bam (sdm) simultation for 100 time steps
#' smd_lep_cal <- bam::sdm_sim(set_A = sparse_mod,
#'                             set_M = adj_mod,
#'                             initial_points = occs_sparse,
#'                             nsteps = 10)
#' #----------------------------------------------------------------------------
#' # Predict species' distribution under suitability change
#' # scenarios (could be climate chage scenarios).
#' #----------------------------------------------------------------------------
#'
#' # Read suitability layers (two suitability change scenarios)
#' layers_path <- system.file("extdata/suit_change",
#'                            package = "bam")
#' niche_mods_stack <- raster::stack(list.files(layers_path,
#'                                              pattern = ".tif$",
#'                                              full.names = TRUE)) > 0.1
#' raster::plot(niche_mods_stack)
#' # Predict
#' new_preds <- predict(object = smd_lep_cal,
#'                      niche_layers = niche_mods_stack,
#'                      nsteps_vec = c(50,100))
#'
#' # Generate the dispersal animation for time period 1 and 2
#' new_preds <- predict(object = smd_lep_cal,
#'                      niche_layers = niche_mods_stack,
#'                      nsteps_vec = c(10,10),
#'                      animate=TRUE,
#'                      filename="/home/l916o895/Desktop/animacion_01.html",
#'                      fmt="HTML")
#' }

methods::setMethod(f = "predict",
                   signature = methods::signature(object = "bam"),
                   function(object,niche_layers,nbgs_vec=NULL,nsteps_vec,
                            animate=FALSE,
                            period_names= NULL,
                            fmt="GIF",filename,
                            png_keyword="sdm_sim",
                            ani.width = 1200,
                            ani.height = 1200,
                            ani.res=300){
                     if(!class(niche_layers) %in% c("RasterLayer",
                                                    "RasterStack",
                                                    "RasterBrick")){
                       stop("niche_layers must be a RasterLayer or RasterStack")
                     }
                     if(class(niche_layers) == "RasterLayer")
                       niche_layers <- raster::stack(niche_layers)
                     n_enm <- raster::nlayers(niche_layers)
                     n_steps_each <- length(nsteps_vec)
                     if(n_steps_each == 1)
                       nsteps_vec <- rep(nsteps_vec,n_enm)
                     if(n_steps_each > 1 && n_enm != n_steps_each){
                       stop(paste("nsteps_vec should have the same",
                                  "length as the number of niche_layers"))
                     }

                     if(!is.null(nbgs_vec)){
                       n_nbgs_each <- length(nbgs_vec)
                       if(n_nbgs_each == 1)
                         nbgs_vec <- rep(nbgs_vec,n_enm)
                       if(n_nbgs_each > 1 && n_enm != n_nbgs_each){
                         stop(paste("nbgs_vec should have the same",
                                    "length as the number of niche_layers"))
                       }
                       ad_mat <- lapply(1:n_nbgs_each, function(x){
                         ad <- bam::adj_mat(modelsparse = object,
                                            nbgs = nbgs_vec[x])
                       })
                     }
                     else{
                       ad_mat <- lapply(1:n_enm, function(x){
                         methods::new("setM",
                                      adj_matrix=object@adj_matrix)
                       })

                     }

                     sim_results <- list(object)

                     initial_points <- object@sdm_sim[[object@sim_steps]]
                     nsteps <- nsteps_vec[1]
                     niche_mod <- niche_layers[[1]]
                     sparse_mod <- bam::model2sparse(niche_mod)
                     sdm <- bam::sdm_sim(set_A = sparse_mod,
                                         set_M = ad_mat[[1]],
                                         initial_points = initial_points,
                                         nsteps = nsteps)


                     periods <- raster::stack(object@bin_model,niche_layers)

                     sim_results[[2]] <-sdm
                     if(n_enm > 1){
                       for(x in 2:length(nsteps_vec)){
                         nsteps <- nsteps_vec[x]
                         niche_mod <- niche_layers[[x]]
                         sparse_mod <- bam::model2sparse(niche_mod)
                         bam_object <- sim_results[[x]]
                         initial_points <- bam_object@sdm_sim[[bam_object@sim_steps]]

                         sdm <- bam::sdm_sim(set_A = sparse_mod,
                                             set_M = ad_mat[[x]],
                                             initial_points = initial_points,
                                             nsteps = nsteps)
                         sim_results[[x+1]] <- sdm
                       }
                     }

                     names(sim_results) <- paste0("time_period_",
                                                  1:length(sim_results))

                     if(animate){

                       nsteps_vec <- c(object@sim_steps,nsteps_vec)
                       nsteps <- sum(nsteps_vec)

                       if(nsteps> 80 && fmt=="GIF"){

                         #contri <- nsteps/length(sim_results)
                         which_steps <- round(seq(1,nsteps,
                                                  along.with = 1:80))
                         step <- max(diff(which_steps))

                        which_stepsL <-  lapply(1:length(nsteps_vec),
                                                function(x){
                          unique(c(seq(1,nsteps_vec[x],step),
                                   nsteps_vec[x]))
                        })

                       }
                       else{
                         which_stepsL <-  lapply(1:length(nsteps_vec),
                                                 function(x){
                                                   seq(1,nsteps_vec[x])
                                                 })
                       }

                       titles <- lapply(1:length(which_stepsL),
                                        function(x) {
                                          if(!is.null(period_names) &&
                                             length(period_names) == length(which_stepsL)){
                                            paste0(period_names[x],paste0(" (t_"),
                                                   which_stepsL[[x]],")")
                                          }
                                          else{
                                            paste0(paste0("Period ",x," (t_"),
                                                   which_stepsL[[x]],")")
                                          }

                                        })
                       titles <- unlist(titles)


                       sdm_st <- 1:length(sim_results) %>%
                         purrr::map(function(x){
                           #nsims <-length(sim_results[[x]]@sdm_sim) -1
                           sdm_ani <- bam::sim2Raster(sim_results[[x]],
                                                      which_steps = which_stepsL[[x]])
                           sdm_ani <- periods[[x]]*1 * (sdm_ani+1)

                           return(sdm_ani)
                         })
                       sdm_st <- raster::stack(sdm_st)
                       names(sdm_st) <- titles

                       fmt <- toupper(fmt)
                       if(!fmt %in% c("GIF",'HTML'))
                         stop("fmt should be GIF or HTML")

                       dir1 <- unlist(strsplit(filename,split = "[/]|[\\]"))
                       filename <- paste0(dir1,collapse = "/")
                       dir2 <- paste0(dir1[1:(length(dir1)-1)],collapse = '/')
                       dir2 <- normalizePath(dir2)
                       if(fmt == "GIF"){
                         animation::ani.options(ani.width = ani.width,
                                                ani.height = ani.height,
                                                ani.res = ani.res)

                         animation::saveGIF({
                           for (i in 1:raster::nlayers(sdm_st)) {
                             maxv <- raster::maxValue(sdm_st[[i]])
                             if(maxv<1.5) colores <- c("#F6F2E5","#0076BE")
                             else colores <- c("#F6F2E5","#0076BE","#03C33F")

                             graphics::par(xpd = FALSE)


                             raster::plot(sdm_st[[i]],main=titles[i],
                                          col=colores,legend=FALSE,
                                          xaxt = 'n',
                                          yaxt = 'n')

                             graphics::par(xpd = TRUE)
                             graphics::legend(
                               "bottom",
                               legend = c("Unsuitable", "Suitable", "Occupied"),
                               fill = c("#F6F2E5","#0076BE","#03C33F"),
                               horiz = TRUE,
                               inset = -0.2,
                               cex = 0.75,
                               bty="n"
                             )

                           }

                         },interval=0.8,ani.width = ani.width,
                         movie.name = filename)
                        }
                       if(fmt == "HTML"){

                         dir3 <- file.path(dir2,paste0("pngs_",png_keyword),
                                           fsep = '/')
                         dir3 <- gsub("[\\]","/",dir3)

                         dir3 <- gsub("[.]","_",dir3)

                         animation::saveHTML({
                           for (i in 1:raster::nlayers(sdm_st)) {
                             maxv <- raster::maxValue(sdm_st[[i]])
                             if(maxv<1.5) colores <- c("#F6F2E5","#0076BE")
                             else colores <- c("#F6F2E5","#0076BE","#03C33F")

                             graphics::par(xpd = FALSE)


                             raster::plot(sdm_st[[i]],main=titles[i],
                                          col=colores,legend=FALSE,
                                          xaxt = 'n',
                                          yaxt = 'n')

                             graphics::par(xpd = TRUE)
                             graphics::legend(
                               "bottom",
                               legend = c("Unsuitable", "Suitable", "Occupied"),
                               fill = c("#F6F2E5","#0076BE","#03C33F"),
                               horiz = TRUE,
                               inset = -0.2,
                               cex = 0.75,
                               bty="n"
                             )
                           }
                         },img.name = png_keyword,
                         imgdir = dir3 ,
                         htmlfile = filename,
                         ani.width=ani.width,
                         ani.height=ani.width,interval=0.1,
                         ani.dev = function(...){grDevices::png(res=ani.res,...)})
                       }


                     }

                     return(sim_results)
                   })



