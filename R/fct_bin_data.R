NBINS <- 300
NBINS_PRECISE <- NBINS / 10

#' @import data.table
#' @importFrom SpatialExperiment spatialCoords
#' @importFrom SummarizedExperiment colData colData<- assay assay<-

get_bin_boundaries <- function(coordinates, nbins) {
  seq(min(coordinates),
      max(coordinates),
      length.out = nbins + 1)
}
map_coordinates_to_bins <- function(coordinates, bin_boundaries) {
  cut(coordinates,
      breaks = bin_boundaries,
      labels = FALSE,
      include.lowest = TRUE)
}

add_bins <- function(spe,
                     nbins_x = NBINS,
                     nbins_y = NBINS) {
  binx <- map_coordinates_to_bins(spatialCoords(spe)[, 1],
                                  get_bin_boundaries(spatialCoords(spe)[, 1], nbins_x))
  biny <- map_coordinates_to_bins(spatialCoords(spe)[, 2],
                                  get_bin_boundaries(spatialCoords(spe)[, 2], nbins_y))

  colData(spe) <- cbind(colData(spe), binx = binx, biny = biny)

  spe
}

#' Convert a Spatial Experiment object to a rasterized image
#'
#' @description
#' This function bins the total counts and returns a rasterized image.
#'
#' @param spe A [SpatialExperiment::SpatialExperiment] object containing the
#'     data to rasterize.
#' @param nbins_x Number of bins to use on the X axis (default is 300).
#' @param nbins_y Number of bins to use on the Y axis (default is 300).
#' @param nas_to_zero Whether or not to convert missing values to 0 (default is
#'     `TRUE`).
#'
#' @returns An array of dimensions `nbins_x` by `nbins_y` containing the binned
#'     counts for each "pixel".
#'
#' @export
convert_spe_to_image <- function(spe,
                                 nbins_x = NBINS,
                                 nbins_y = NBINS,
                                 nas_to_zero = TRUE) {
  add_bins(spe, nbins_x = nbins_x, nbins_y = nbins_y) |>
    convert_binned_spe_to_image()
}

convert_binned_spe_to_image <- function(spe,
                                 replace_nas = 0) {
  binned_counts <- data.table::as.data.table(t(assay(spe))) |>
    cbind(data.table::as.data.table(colData(spe)[, c("binx", "biny")]))

  counts_image <- binned_counts[
    ,
    lapply(.SD, sum),
    .SDcols = -c("binx", "biny"),
    by = .(binx, biny)
  ][
    CJ(biny = seq_len(max(binned_counts$biny)),
       binx = seq_len(max(binned_counts$binx))),
    .(rowSums(.SD)),
    on = .(binx, biny),
    .SDcols = -c("binx", "biny")
  ] |>
    unlist() |>
    array(dim = c(max(binned_counts$binx),
                  max(binned_counts$biny)),
          dimnames = list(x = seq_len(max(binned_counts$binx)),
                          y = seq_len(max(binned_counts$biny))))

  if (isFALSE(replace_nas))
    return(counts_image)
  replace(counts_image, is.na(counts_image), replace_nas)
}
