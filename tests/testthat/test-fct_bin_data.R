test_spe <- SpatialExperiment::SpatialExperiment(
  assay = data.frame(a = (1:100)/10, b = (1:100)/6) |> t(),
  spatialCoords = data.frame(x = 1:100, y = 1:100) |> as.matrix()
)

test_that("Bin boundaries work", {
  expect_equal(get_bin_boundaries(2:500, 2), c(2, 251, 500))
  expect_equal(get_bin_boundaries(1:10, 4), c(1, 3.25, 5.5, 7.75, 10))
})

test_that("Mapping to bins works", {
  expect_equal(map_coordinates_to_bins(c(1,3,4,5,6,5,10), c(1, 5, 7, 10)),
               c(1, 1, 1, 1, 2, 1, 3))
})


