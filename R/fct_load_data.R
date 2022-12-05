#' @import data.table
load_as_spatial_experiment <- function(
    counts_file,
    metadata_file,
    metadata = NULL,
    id_columns = c("fov", "cell_ID"),
    coord_names = c("CenterX_global_px",
                    "CenterY_global_px"),
    ...
) {
  counts <- data.table::fread(counts_file)
  if (is.null(metadata))
    metadata <- data.table::fread(metadata_file)

  if (!is.null(id_columns)) {
    counts <- merge(metadata[, ..id_columns], counts)
    metadata <- merge(counts[, ..id_columns], metadata)
    counts <- counts[, -..id_columns]
  }

  SpatialExperiment::SpatialExperiment(
    assay = t(counts),
    colData = metadata,
    spatialCoordsNames = coord_names,
    ...
  )
}

load_cosmx_as_spatial_experiment <- function(
    counts_file,
    metadata_file,
    ...
) {
  load_as_spatial_experiment(
    counts_file,
    metadata_file,
    id_columns = c("fov", "cell_ID"),
    coord_names = c("CenterX_global_px",
                    "CenterY_global_px"),
    ...
  )
}

load_merscope_as_spatial_experiment <- function(
    counts_file,
    metadata_file,
    ...
) {
  load_as_spatial_experiment(
    counts_file,
    metadata = data.table::fread(metadata_file) |>
      data.table::setnames("V1", "cell"),
    id_columns = c("cell"),
    coord_names = c("center_x", "center_y"),
    ...
  )
}

load_xenium_as_spatial_experiment <- function(
    counts_file,
    metadata_file,
    ...
) {
  load_as_spatial_experiment(
    counts_file,
    metadata_file,
    id_columns = NULL,
    coord_names = c("x_centroid", "y_centroid"),
    ...
  )
}
