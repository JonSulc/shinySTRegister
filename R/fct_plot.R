#' Plots the image produced by rasterizing the counts data.
#'
#' @param counts_image The array of binned counts, as returned by
#'     [convert_spe_to_image].
#' @param flip Whether or not to reverse the x and/or y axes. Should be one of
#'     "x", "X", "y", "Y", `TRUE`, or `FALSE` (default).
#' @param points Optional data.frame containing the coordinates of points to
#'     plot on the image.
#' @param log1p If `TRUE`, the counts are normalized prior to plotting. Default
#'     is `FALSE`.
#' @param minimal Whether to remove the axis text and legend. Default is
#'     `FALSE`.
#'
#' @returns A ggplot object.
#'
#' @import ggplot2
#' @export
plot_counts <- function(counts_image,
                        flip = FALSE,
                        points = data.frame(),
                        log1p = FALSE,
                        minimal = FALSE) {
  counts_plot <- reshape2::melt(counts_image) |>
    ggplot(aes(x = x, y = y, fill = value)) +
    geom_raster() +
    viridis::scale_fill_viridis(trans = ifelse(log1p, "log1p", "identity")) +
    theme_minimal()
  if (nrow(points) != 0)
    points <- points[points[, 1] < max(as.numeric(dimnames(counts_image)$x)) &
                       points[, 1] > min(as.numeric(dimnames(counts_image)$x)) &
                       points[, 2] < max(as.numeric(dimnames(counts_image)$y)) &
                       points[, 2] > min(as.numeric(dimnames(counts_image)$y)), ]

  if (nrow(points) != 0)
    counts_plot <- counts_plot +
      geom_point(data = points |> cbind(point = letters[1:nrow(points)]),
                 aes(fill = NULL),
                 shape = 3, color = "white", size = 3, stroke = 2)
  if (minimal)
    counts_plot <- counts_plot +
      theme(axis.text = element_blank(), legend.position = "none")

  if (flip %in% c("x", "X") || isTRUE(flip))
    counts_plot <- counts_plot + scale_x_reverse()
  if (flip %in% c("y", "Y") || isTRUE(flip))
    counts_plot <- counts_plot + scale_y_reverse()

  counts_plot
}

ensure_interval_within_image <- function(center_position, size, image_dim) {
  stopifnot(image_dim > size)
  max(center_position, size / 2) |>
    min(image_dim - size / 2)
}

subset_image <- function(counts_image,
                         center,
                         size_x = NBINS_PRECISE,
                         size_y = size_x) {
  x <- ensure_interval_within_image(center$x, size_x+1, ncol(counts_image))
  y <- ensure_interval_within_image(center$y, size_y+1, nrow(counts_image))

  counts_image[x - size_x/2 + seq_len(size_x),
               y - size_y/2 + seq_len(size_y)]
}
