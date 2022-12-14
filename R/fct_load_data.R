#' Load a counts file and metadata file as a SpatialExperiment object
#'
#' @description
#' `load_as_spatial_experiment` loads the corresponding files from disk,
#' combines them based on the provided identifying columns and returns a
#' SpatialExperiment object.
#'
#' @param counts_file Path to the file containing the feature counts per
#'     spot/cell.
#' @param metadata_file Path to the file containing the spot/cell metadata. Must
#'     include the spatial coordinates.
#' @param metadata Object containing the spot/cell metadata. If one is supplied,
#'     `metadata_file` is ignored.
#' @param id_columns Names of the columns in common between `counts` and
#'     `metadata`, used for aligning the spots/cells. Default is
#'     `c("fov", "cell_ID")`.
#' @param coord_names A vector of length 2, containing the names of the columns
#'     containing the x and y coordinates in the metadata. Default is
#'     `c("CenterX_global_px", "CenterY_global_px")` (which are the default
#'     Cosmx column names).
#' @param ... Additional parameters passed to the
#'     [SpatialExperiment::SpatialExperiment] constructor.
#'
#' @returns A `SpatialExperiment` object.
#'
#' @import data.table
#' @export
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

#' @export
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

#' @export
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

#' @export
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
