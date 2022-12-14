#' @import shiny
app_ui <- function() {
  fluidPage(
    fluidRow(
      radioButtons(
        "point_selected",
        "Select point to set",
        choices = letters[1:3],
        inline = TRUE
      )
    ),
    fluidRow(
      column(5,
             plotOutput("source_image", click = "source_location"),
             fluidRow(
               column(6, checkboxInput("flip_source", "Rotate image")),
               column(6, checkboxInput("log1p_source", "Normalize"))
             )
      ),
      column(2,
        plotOutput("source_precise",
                   click = "source_location_precise",
                   height = "200px"),
        plotOutput("reference_precise",
                   click = "reference_location_precise",
                   height = "200px")
      ),
      column(5,
             plotOutput("reference_image", click = "reference_location"),
             fluidRow(
               column(6, checkboxInput("flip_reference", "Rotate image")),
               column(6, checkboxInput("log1p_reference", "Normalize"))
             )
      )
    ),
    fluidRow(
      column(
        4,
        DT::DTOutput("affine_matrix"),
        actionButton("align", "Align"),
        fluidRow(
          column(4,
                 textInput(
                   "rds_filename",
                   "Filename (RDS)",
                   value = "output.rds"
                 )
          ),
          column(8,
                 actionButton(
                   "save_rds_button",
                   "Save"
                 )
          )
        )
      ),
      column(
        8,
        plotOutput("overlay_image")
      )
    )
  )
}

app_server <- function(source_spe,
                       reference_spe,
                       nbins = NBINS,
                       nbins_precise = NBINS_PRECISE) {
  function(input, output, session) {
    source_image <- reactive({
      source_spe |>
        convert_spe_to_image(nbins_x = nbins, nbins_y = nbins)
    }) |>
      bindEvent(source_spe, nbins)
    reference_image <- reactive({
      reference_spe |>
        convert_spe_to_image(nbins_x = nbins, nbins_y = nbins)
    }) |>
      bindEvent(reference_spe, nbins)

    output$source_image <- renderPlot({
      plot_counts(source_image(), input$flip_source,
                  points = source_coordinates$xy,
                  log1p = input$log1p_source)
    })
    output$reference_image <- renderPlot({
      plot_counts(reference_image(), input$flip_reference,
                  points = reference_coordinates$xy,
                  log1p = input$log1p_reference)
    })

    source_coordinates <- reactiveValues(
      xy = data.frame(x = c(.75*nbins, .75*nbins, .25*nbins),
                      y = c(.25*nbins, .75*nbins, .5*nbins),
                      row.names = letters[1:3]),
      approx_xy = data.frame(x = c(.75*nbins, .75*nbins, .25*nbins),
                             y = c(.25*nbins, .75*nbins, .5*nbins),
                             row.names = letters[1:3])
    )
    observe({
      source_coordinates$approx_xy[input$point_selected, ] <-
        source_coordinates$xy[input$point_selected, ] <-
        input$source_location[c("x", "y")] |>
        unlist()
    }) |>
      bindEvent(input$source_location)
    observe({
      source_coordinates$xy[input$point_selected, ] <-
        input$source_location_precise[c("x", "y")] |>
        unlist()
    }) |>
      bindEvent(input$source_location_precise)
    output$source_precise <- renderPlot({
      source_image() |>
        subset_image(source_coordinates$approx_xy[input$point_selected, ],
                     nbins_precise) |>
        plot_counts(input$flip_source,
                    points = source_coordinates$xy,
                    log1p = input$log1p_source,
                    minimal = TRUE)
    })

    reference_coordinates <- reactiveValues(
      xy = data.frame(x = c(.75*nbins, .75*nbins, .25*nbins),
                      y = c(.25*nbins, .75*nbins, .5*nbins),
                      row.names = letters[1:3]),
      approx_xy = data.frame(x = c(.75*nbins, .75*nbins, .25*nbins),
                             y = c(.25*nbins, .75*nbins, .5*nbins),
                             row.names = letters[1:3])
    )
    observe({
      reference_coordinates$approx_xy[input$point_selected, ] <-
      reference_coordinates$xy[input$point_selected, ] <-
        input$reference_location[c("x", "y")] |>
        unlist()
    }) |>
      bindEvent(input$reference_location)
    observe({
      reference_coordinates$xy[input$point_selected, ] <-
        input$reference_location_precise[c("x", "y")] |>
        unlist()
    }) |>
      bindEvent(input$reference_location_precise)
    output$reference_precise <- renderPlot({
      reference_image() |>
        subset_image(reference_coordinates$approx_xy[input$point_selected, ],
                     nbins_precise) |>
        plot_counts(input$flip_reference,
                    points = reference_coordinates$xy,
                    log1p = input$log1p_reference,
                    minimal = TRUE)
    })

    affine_matrix <- reactive({
      find_affine_matrix(source_coordinates$xy,
                         reference_coordinates$xy,
                         source_spe = source_spe,
                         reference_spe = reference_spe,
                         nbins_x = nbins)
    }) |>
      bindEvent(source_coordinates$xy, reference_coordinates$xy)

    output$affine_matrix <- DT::renderDT({
      affine_matrix()
    })
    source_rotated <- reactive({
      source_spe |>
        apply_affine_to_spe(affine_matrix())
    }) |>
      bindEvent(input$align)

    output$overlay_image <- renderPlot({
      req(source_rotated())
      source_rotated() |>
        SingleCellExperiment::cbind(reference_spe) |>
        convert_spe_to_image(nbins_x = nbins, nbins_y = nbins) |>
        plot_counts()
    })

    observe({
      source_rotated() |>
        saveRDS(input$rds_filename)
    }) |>
      bindEvent(input$save_rds_button)
  }
}

#' Launch the shiny interface to align two spatial transcriptomics samples
#'
#' @description
#' This launches a shiny interface displaying a rasterized version of the
#' samples (total counts per bin). Three points must be placed on the image of
#' the sample to rotate, with their corresponding projections on the reference
#' image. The coordinates of the first sample are then transformed such that the
#' points are lined up in both. The affine matrix and the resulting
#' superposition of the samples are displayed. When the result is satisfactory,
#' the `SpatialExperiment` object with the transformed coordinates can be saved
#' as an `RDS` file.
#'
#' @param source A [SpatialExperiment::SpatialExperiment] object containing the
#'     sample to rotate.
#' @param reference A [SpatialExperiment::SpatialExperiment] object containing
#'     the reference sample.
#' @param nbins The number of bins on each axis. This is mainly for the image
#'     display for the definition of the points and has no impact on the final
#'     object.
#'
#' @returns A [SpatialExperiment::SpatialExperiment] object containing the
#'     rotated sample with the new coordinates.
#'
#' @export
shiny_register <- function(
    source,
    reference,
    nbins = 200
) {
  runApp(list(ui = app_ui(), server = app_server(source, reference, nbins)))
}
