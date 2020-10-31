# R classes for bam package
# December 2019
# Version 0.3.2
# Licence GPL v3

#' S4 classes to organize data and results of \code{bam} objects
#' @aliases coordinates-class
#' @importFrom methods new
#' @slot coordinates A two column matrix with coordinates
#' @slot eigen_vec Eigen vector of adjacency matrix
#' @slot eigen_val Eigen value of adjacency matrix
#' slot g_model A raster representing the geographic area
#' slot g_sparse A sparse matrix of the geographic area
#' @author Luis Osorio-Olvera & Jorge Soberón
#' @exportClass g_area
#' @export
#'
g_area <- methods::setClass(Class = "g_area",
                            slots = c(coordinates="matrix",
                                      eigen_vec = "matrix",
                                      eigen_val = "numeric"
                                      #g_raster = "RasterLayer",
                                      #g_sparse = "dgCMatrix"
                                      ))

#' Class for the A set of the BAM diagram
#' @aliases A-class
#' @importClassesFrom raster RasterLayer
#' @importClassesFrom raster RasterStack
#' @importClassesFrom Matrix dgCMatrix
#' @importFrom methods new
#' @author Luis Osorio-Olvera & Jorge Soberón
#' @slot bin_model A niche model in raster format
#' @slot cellIDs A numeric vector with the IDs of the cells with prediction values
#' @slot sparse_model A niche model in sparse matrix format
#' @exportClass setA
#' @export
#' @examples
#' showClass("setA")
#' @rdname setA
setA <- methods::setClass(Class = "setA",
                          slots = c(bin_model = "RasterLayer",
                                    cellIDs = "numeric",
                                    sparse_model = "dgCMatrix"),
                          contains = "g_area",
                          validity = function(object){
                            bin_class <- class(object@bin_model)
                            if(!bin_class %in% c("RasterLayer",
                                                 "RasterStack")){
                              stop("bin_model should be a raster layer
                                   or raster stack")
                            }})

#' Class for the M set of the \code{bam} digram
#' @aliases M-class
#' @slot zone A raster representing the study area (could be the niche model)
#' @slot initial_points A presence-absence vector with species' occurrences
#' @slot n_initial_points Number of initial points used to start the dispersal process
#' @importClassesFrom raster RasterLayer
#' @importClassesFrom raster RasterStack
#' @importClassesFrom Matrix dgCMatrix
#' @author Luis Osorio-Olvera & Jorge Soberón
#' @exportClass setM
#' @export
#' @examples
#' showClass("setM")
#' @rdname setM

setM <- methods::setClass(Class = "setM",
                          contains = "g_area",
                          slots =  c(adj_matrix = "dgCMatrix",
                                     initial_points = "dgCMatrix"))

#' Class  \code{bam} digram
#' @aliases bam-class
#' @author Luis Osorio-Olvera & Jorge Soberón
#' @slot sdm_sim A list of sparse vectors representing the area occupied
#' @slot sim_steps Number of simulation steps
#' by the species
#' @exportClass bam
#' @export
#'
bam <- methods::setClass(Class = "bam",
                         contains = c("setA","setM"),
                         slots = c(sdm_sim = "list",
                                   sim_steps="numeric"))


#' Class  \code{community_sim} digram
#' @aliases community-class
#' @author Luis Osorio-Olvera & Jorge Soberón
#' @slot community_sim A list of sparse vectors representing the area occupied
#' by the species

community_bam <- methods::setClass(Class = "community_bam",
                                   #contains = c("bam"),
                                   slots = c(community_sim = "list"))


#' Class  \code{pam} Presence-Absence Matrix
#' @aliases PAM
#' @author Luis Osorio-Olvera & Jorge Soberón
#' @slot pams A list of sparse matices representing Presence-Absecen Matrix for
#' each simulation time
#' @slot which_steps Simulation steps
#' @slot sp_names Names of species in the PAM
#' @slot grid Raster grid of the studied area
#' @slot cellIDs Cells with ids of the PAM sites
#' @exportClass pam
#' @export
#'
pam <- methods::setClass(Class = "pam",
                         #contains = c("bam"),
                         slots = c(pams = "list",
                                   which_steps= "numeric",
                                   sp_names = "character",
                                   grid = "RasterLayer",
                                   cellIDs = "numeric"))




leaflet <- methods::setClass(Class = "leaflet")

#' Class  \code{csd}
#' @aliases csd-class
#' @author Luis Osorio-Olvera & Jorge Soberón
#' @slot connections A data.frame with four columns: x, y, clusterID and cluster_size
#' @slot interactive_map A leaflet map with markers showing the geographical clusters
#' @slot raster_map A raster map with cluster IDs as values.
#' @exportClass csd
#' @export
#' @import leaflet
csd <- methods::setClass(Class = "csd",
                         representation(connections="data.frame",
                                        interactive_map = "leaflet",
                                        raster_map = "RasterLayer"))




#setClassUnion("sparse_or_matrix", c("dgeMatrix", "matrix"))


#' Class  \code{bioindex_sparse}
#' @aliases biodiversity_index_sparse
#' @author Luis Osorio-Olvera & Jorge Soberón
#' @slot alpha A sparse matrix with the richness of species per site
#' @slot omega A sparse matrix with the range size of every species
#' @slot dispersion_field A sparse matrix with the set of ranges of all species that occur in at each locality
#' @exportClass bioindex
#' @export

bioindex_sparse <- methods::setClass(Class = "bioindex_sparse",
                              representation (
                                alpha = "dgeMatrix",
                                omega = "dgeMatrix",
                                dispersion_field = "dgeMatrix"
                              ))


#' Class  \code{bioindex}
#' @aliases biodiversity_index
#' @author Luis Osorio-Olvera & Jorge Soberón
#' @slot alpha A matrix with the richness of species per site
#' @slot omega A matrix with the range size of every species
#' @slot dispersion_field A matrix with the set of ranges of all species that occur in at each locality
#' @exportClass bioindex
#' @export

bioindex <- methods::setClass(Class = "bioindex",
                                     representation (
                                       alpha = "matrix",
                                       omega = "matrix",
                                       dispersion_field = "matrix"
                                     ))
