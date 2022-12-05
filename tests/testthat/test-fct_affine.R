library(testthat)
test_that("Finding the affine transform matrix works", {
  # Scaling
  expect_equal(find_affine_matrix(matrix(c(0, 1, 1, 0, 1, 1), nrow = 2),
                                  matrix(c(0, 2, 3, 0, 3, 2), nrow = 2)),
               matrix(c(3, 0, 0, 2, 0, 0), nrow = 2))

  # Translation
  expect_equal(find_affine_matrix(matrix(c(0, 1, 1, 0, 1, 1), nrow = 2),
                                  matrix(c(1, 3, 2, 2, 2, 3), nrow = 2)),
               matrix(c(1, 0, 0, 1, 1, 2), nrow = 2))

  # Rotation
  # theta <- pi/3
  for (theta in runif(3))
    expect_equal(find_affine_matrix(matrix(c(0, 1, 1, 0, 1, 1), nrow = 2),
                                    matrix(c(-sin(theta), cos(theta),
                                             cos(theta), sin(theta),
                                             cos(theta)-sin(theta), sin(theta) + cos(theta)), nrow = 2)),
                 matrix(c(cos(theta), sin(theta), -sin(theta), cos(theta), 0, 0), nrow = 2))

  # Shear
  k <- .5
  expect_equal(find_affine_matrix(matrix(c(0, 1, 1, 0,   1, 1), nrow = 2),
                                  matrix(c(k, 1, 1, 0, 1+k, 1), nrow = 2)),
               matrix(c(1, 0, k, 1, 0, 0), nrow = 2))
  expect_equal(find_affine_matrix(matrix(c(0, 1, 1, 0, 1, 1), nrow = 2),
                                  matrix(c(0, 1, 1, k, 1, 1+k), nrow = 2)),
               matrix(c(1, k, 0, 1, 0, 0), nrow = 2))
})
