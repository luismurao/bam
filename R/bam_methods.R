
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


#' Show information in pam class \pkg{bam}.
#' @importFrom methods new
#' @param object An object of class bioindiex_sparse
#' @rdname show
#' @export

methods::setMethod(f = "show",
                   signature = "bioindex_sparse",
                   function(object) {
                     slotsin <- methods::slotNames(object)
                     cat("Object of class bioindex it contains",
                         length(slotsin),"slots: \n\n")
                     cat("@alpha: A sparse matrix with the richness of species per site \n\n")
                     print(object@alpha,6)
                     cat("\n")
                     cat("@omega: A sparse matrix with the range size of every species \n\n")
                     print(object@omega,6)
                     cat("\n")
                     cat("@dispersion_field: A sparse with the set of ranges of all species that occur in at each locality \n\n")
                     print(object@omega,6)

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


#' Show information in diversity_range class \pkg{bam}.
#' @importFrom methods new
#' @param object An object of class diversity_range
#' @rdname show
#' @export

methods::setMethod(f = "show",
                   signature = "diversity_range",
                   function(object) {
                     slotsin <- methods::slotNames(object)
                     cat("diversity_range object it contains",
                         length(slotsin),"slots \n\n")


                     cat("@alpha: A column vector of size ",object@nsites,
                         " with values of alpha diversity at each site",
                         "\n\n")
                     print(head(object@alpha))

                     cat("@alpha_raster: Alpha diversity raster",
                         "\n\n")
                     print(object@alpha_raster)

                     cat("@omega: A comun vector of size ",object@nsps,
                         " with the range size of each species",
                         "\n\n")
                     print(head(object@omega))
                     cat("@dispersion_field: A column vector of size ",object@nsites ,
                         "with values of dispersion field at each site",
                         "\n\n")
                     print(head(object@dispersion_field))
                     cat("@dispersion_field_raster: Dispersion field raster",
                         "\n\n")
                     print(object@dispersion_field_raster)
                     cat("@null_dispersion_field_dist: null dispersion field distribution. ")
                     cat("It is matrix of\n n =",object@nsites,"sites x","m =",
                         object@n_iterations, "simulations",
                         "used to generate random values of dispersion field" ,
                         "\n\n")
                     if(ncol(object@null_dispersion_field_dist)>=2){
                       print(head(object@null_dispersion_field_dist[1:2,1:2]))
                     }

                     cat("@diversity_range_raster: Raster with diversity range categories",
                         "\n\n")
                     print(object@diversity_range_raster)
                     cat("@xy_coordinates: Geographical coordinates of the sites",
                         "\n\n")
                     print(head(object@xy_coordinates))

                   })


#if (!isGeneric("plot")) {setGeneric("plot", function(x,y,...)standardGeneric("plot"))}


#' Plot method for objects of class diversity_range \pkg{bam}.
#' @importFrom methods new
#' @param x An object of class diversity_range

#' @param plot_type Plot type: possible options: "diversity_range" (range-diversity plot),
#'                  "diversity_range_map" (a raster map with diversity_range categories),
#'                  "alpha" (a raster mapa with alpha diversity values),
#'                  "dispersion_field" (a raster with dispersion field)
#' @param legend Logical. If TRUE the legend of the categorical diversity range values will appear.
#' @param legend_position Lengend position.
#' @param xlab x label
#' @param ylab y label
#' @param col Plot colors.
#' @param pch Patch type.
#' @param pch_legend Patch type for legends.
#' @param ... Graphical parameters. Any argument that can be passed to base plot,
#'            such as axes=FALSE, main='title', ylab='latitude'.
#' @rdname plot
#' @export

methods::setMethod(f = "plot",
                   signature = c(x="diversity_range"),
                   function(x,xlab=NULL,plot_type="diversity_range",legend=TRUE,
                            legend_position = "bottomright",
                            ylab=NULL,col=NULL,pch=NULL,pch_legend=19,...) {
                     if(inherits(x, 'diversity_range')){
                       #slotsin <- methods::slotNames(x)
                       #zeros_alpha <- which(x@alpha ==0)
                       #zeros_disperfield <- which(x@dispersion_field==0)
                       nsites <- x@nsites
                       nsps <- x@nsps

                       if(x@nsites>0){
                         #cols <- c("#000000","#F6BDC0",
                         #         "#F07470","#BBDFFA",
                          #         "#DC1C13","#6987D5",
                           #        "#1727AE")

                         cols <- c("#000000","#F6BDC0",
                                   "#F1A13A","#BBDFFA",
                                   "#DC1C13","#6987D5",
                                   "#1727AE")

                         names(cols) <- c("Random","HE/LR",
                                          "HE/IR","LE/LR",
                                          "HE/HR","LE/IR",
                                          "LE/HR")
                         COLORES<- cols
                         COLORES <- COLORES[c(1:3,5,4,6,7)]
                         if(is.null(pch)){
                           pch <- 19
                         }

                         if("diversity_range" %in% plot_type){
                           alpha_norm <- x@alpha/nsps
                           dispersion_field <- x@dispersion_field/nsps
                           alpha_st <- x@alpha/x@nsps
                           betty <- round(1/mean(alpha_st),3);
                           fistprom <-x@dispersion_field/x@alpha;
                           rho  <-alpha_st*(fistprom-1/betty);
                           am <- min(alpha_st)
                           aM <- max(alpha_st)
                           fm <- min(x@dispersion_field,na.rm = T)/nsites
                           fM <- max(x@dispersion_field,na.rm =T)/nsites

                           rhom  <- min(rho,na.rm = T)
                           rhoM <- max(rho,na.rm = T)

                           vx <- c(am,aM,aM,am)

                           vy <- c(rhom+am/betty,rhom+aM/betty,
                                   rhoM+aM/betty,rhoM+am/betty)

                           if(is.null(col)){
                             col <- x@diversity_range_colors
                           }
                           if(is.null(xlab)){
                             xlab <- expression(alpha)
                           }
                           if(is.null(ylab)){
                             ylab <- "Dispersion field"
                           }
                           xmin1<- 0.97*min(vx)
                           xmax1 <- 1.03*max(vx)
                           ymin1 <- 0.97*min(vy)
                           ymax1 <- 1.03*max(vy)
                           plot(alpha_norm,dispersion_field,xlim=c(xmin1,xmax1),
                                ylim=c(ymin1,ymax1),
                                xlab=xlab,ylab=ylab,pch=pch,col=col,...)
                           graphics::lines(graphics::polygon(vx,vy));
                           if(legend){
                             graphics::legend(legend_position,
                                              legend = names(COLORES),
                                              pch =pch_legend,
                                              col = COLORES,bty = "n",...)
                           }
                         }

                         if("alpha" %in% plot_type && raster::hasValues(x@alpha_raster)){
                           raster::plot(x@alpha_raster,...)
                         }

                         if("dispersion_field" %in% plot_type &&
                            raster::hasValues(x@dispersion_field_raster)){
                           raster::plot(x@dispersion_field_raster,...)
                         }

                         if("diversity_range_map" %in% plot_type &&
                            raster::hasValues(x@diversity_range_raster)){
                           if(is.null(col)){
                             col1 <- cols
                           }
                           else{
                             col1 <- rev(grDevices::terrain.colors(7))
                           }
                           randiv <- x@diversity_range_raster
                           raster::plot(randiv,col=col1,legend=FALSE,...)
                           if(legend){
                             graphics::legend(legend_position,legend = names(COLORES),
                                              pch=15,col = COLORES,bty = "n",...)
                           }
                         }
                         if("diversity_range_map" %in% plot_type &&
                            !raster::hasValues(x@diversity_range_raster) &&
                            nrow(x@xy_coordinates) == nsites){

                           plot(x@xy_coordinates,col=x@diversity_range_colors,pch=15)
                           graphics::legend("bottomleft",legend = names(COLORES),
                                            pch=15,col = COLORES,bty = "n")
                         }
                       }

                     }

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
#' @param bg_color Color for unsuitable pixels. Default "#F6F2E5".
#' @param suit_color Color for suitable pixels. Default "#0076BE".
#' @param occupied_color Color for occupied pixels. Default "#03C33F".
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
                            bg_color="#F6F2E5",
                            suit_color="#0076BE",
                            occupied_color="#03C33F",
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
                             if(maxv<1.5) colores <- c(bg_color,suit_color)
                             else colores <- c(bg_color,suit_color,occupied_color)

                             graphics::par(xpd = FALSE)


                             raster::plot(sdm_st[[i]],main=titles[i],
                                          col=colores,legend=FALSE,
                                          xaxt = 'n',
                                          yaxt = 'n')

                             graphics::par(xpd = TRUE)
                             graphics::legend(
                               "bottom",
                               legend = c("Unsuitable", "Suitable", "Occupied"),
                               fill = colores,
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
                             if(maxv<1.5) colores <- c(bg_color,suit_color)
                             else colores <- c(bg_color,suit_color,occupied_color)

                             graphics::par(xpd = FALSE)


                             raster::plot(sdm_st[[i]],main=titles[i],
                                          col=colores,legend=FALSE,
                                          xaxt = 'n',
                                          yaxt = 'n')

                             graphics::par(xpd = TRUE)
                             graphics::legend(
                               "bottom",
                               legend = c("Unsuitable", "Suitable", "Occupied"),
                               fill = colores,
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



