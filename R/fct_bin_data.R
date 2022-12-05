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

# sum_counts_per_feature_per_bin <- function(spe) {
#   counts <- as.data.table(t(assay(spe)))[
#     ,
#     `:=`(binx = colData(spe)$binx,
#          biny = colData(spe)$biny)
#   ][
#     ,
#     lapply(.SD, sum),
#     .SDcols = -c("binx", "biny"),
#     by = .(binx, biny)
#   ]
#
#   SpatialExperiment::SpatialExperiment(
#     assay = t(counts[, -c("binx", "biny")]),
#     spatialCoords = as.matrix(counts[, .(binx, biny)])
#   )
# }

# sum_features_per_spot <- function(spe) {
#   counts <- assay(spe) |>
#     t() |>
#     as.data.table() |>
#     cbind(spatialCoords(spe))
#
#   counts <- counts[
#     ,
#     .(count = rowSums(.SD)),
#     by = .(binx, biny),
#     .SDcols = -c("binx", "biny")
#   ] |>
#     t()
#
#   SpatialExperiment::SpatialExperiment(
#     assay = counts,
#     spatialCoords = spatialCoords(spe),
#     colData = colData(spe)
#   )
# }


# bin_counts <- function(spe,
#                        nbins_x = NBINS,
#                        nbins_y = NBINS) {
#   add_bins(spe, nbins_x = nbins_x, nbins_y = nbins_y) |>
#     sum_counts_per_feature_per_bin() |>
#     sum_features_per_spot() |>
#     add_missing_bins_to_spe(nbins_x = nbins_x, nbins_y = nbins_y)
# }

# add_missing_bins_to_spe <- function(spe,
#                                     nbins_x = NBINS,
#                                     nbins_y = NBINS ) {
#   counts <- assay(spe) |>
#     t() |>
#     as.data.table() |>
#     cbind(spatialCoords(spe))
#
#   counts <- counts[
#     CJ(binx = seq_len(nbins_x),
#        biny = seq_len(nbins_y)),
#     on = .(binx, biny)
#   ]
#
#   SpatialExperiment::SpatialExperiment(
#     assay = counts[, .SD, .SDcols = -c("binx", "biny")] |> t(),
#     spatialCoords = counts[, .(binx, biny)] |>
#       as.matrix()
#   )
# }

convert_spe_to_image <- function(spe,
                                 nbins_x = NBINS,
                                 nbins_y = NBINS,
                                 nas_to_zero = TRUE) {
  add_bins(spe, nbins_x = nbins_x, nbins_y = nbins_y) |>
    convert_binned_spe_to_image()
}

# bin_transcript_count_per_sample <- function(
#     spe,
#     nbins_x = NBINS,
#     nbins_y = NBINS
# ) {
#   spe <- add_bins(spe, nbins_x, nbins_y)
#
#   unique(colData(spe)$sample_id) |>
#     lapply(
#       function(sample_id)
#         sum_counts_per_bin(spe[, colData(spe)$sample_id == sample_id])
#     )
# }

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
