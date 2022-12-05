#' @importFrom SpatialExperiment spatialCoords<-

rescale <- function(binned_coordinates,
                    full_coordinates,
                    nbins_x = 200,
                    nbins_y = nbins_x) {
  binned_coordinates[, 1] <-
    quantile(range(full_coordinates[, 1]), binned_coordinates[, 1]/nbins_x)
  binned_coordinates[, 2] <-
    quantile(range(full_coordinates[, 2]), binned_coordinates[, 2]/nbins_y)
  binned_coordinates
}

find_affine_matrix <- function(source_coordinates,
                               reference_coordinates,
                               source_spe = NULL,
                               reference_spe = NULL,
                               nbins_x = 200,
                               nbins_y = nbins_x) {
  if (!is.null(source_spe))
    source_coordinates <- rescale(source_coordinates,
                                  spatialCoords(source_spe),
                                  nbins_x = nbins_x,
                                  nbins_y = nbins_y)
  if (!is.null(reference_spe))
    reference_coordinates <- rescale(reference_coordinates,
                                  spatialCoords(reference_spe),
                                  nbins_x = nbins_x,
                                  nbins_y = nbins_y)

  if (is.data.frame(reference_coordinates))
    reference_coordinates <- t(reference_coordinates) |> as.matrix()
  if (is.data.frame(source_coordinates))
    source_coordinates <- t(source_coordinates) |> as.matrix()

  reference_coordinates %*% MASS::ginv(source_coordinates |> rbind(1))
}

apply_affine_to_spe <- function(spe,
                                affine_matrix) {
  spatialCoords(spe) <- cbind(spatialCoords(spe), u = 1) %*% t(affine_matrix)

  spe
}
