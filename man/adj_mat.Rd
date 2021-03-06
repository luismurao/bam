% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/adj_mat.R
\name{adj_mat}
\alias{adj_mat}
\title{adj_mat: Function to compute the adjacency matrix of an area.}
\usage{
adj_mat(modelsparse, ngbs = 1, eigen_sys = FALSE, which_eigs = 1)
}
\arguments{
\item{modelsparse}{A setA object returned by the function \code{\link[bam]{model2sparse}}.}

\item{ngbs}{Numeric. Number of neighbors (see details).}

\item{eigen_sys}{Logical. If TRUE the eigen analsys of the adjancency matrix will be returned.}

\item{which_eigs}{Numeric. Which eigen value and eigen vector will be returned.}
}
\value{
Returns an adjacency matrix of class sparseMatrix of n+2 x n columns (n number of the non-NA cells of grid_base) with the coordinates of the non-NA cells of grid_base.
}
\description{
Creates an adjacency matrix of an area of interest. This could be the accessible area (M) of a species or any geographic region of interest.
}
\details{
The grid_base raster object is the area where the dispersal process will occur.
The number of neighbors depends on the dispersal abilities of the species and the spatial resolution
of the grid_base; for example, a species's with big dispersal abilities will move throughout more
than 1 km^2 per day, so the idea is to give an approximate number of moving neighbors (pixels) per
unit of time.
}
\examples{
\dontrun{
data("wrld_simpl", package = "maptools")
mx <- wrld_simpl[wrld_simpl$NAME=="Mexico",]
mx_grid <- shape2Grid(mx,0.5)
mx_sparse <- bam::model2sparse(mx_grid)
adj_mx <- adj_mat(mx_sparse,ngbs=1)
# Adjacency matrix from a niche model
model_path <- system.file("extdata/Lepus_californicus_cont.tif",
                          package = "bam")
model <- raster::raster(model_path)

sparse_mod <- bam::model2sparse(model,threshold=0.05)
adj_mod <- adj_mat(sparse_mod,ngbs=1)
}
}
