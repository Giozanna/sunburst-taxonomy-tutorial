#' app.R
#'
#' Minimal Shiny wrapper around the sunburst pipeline from
#' sunburst_tutorial.Rmd: adds a dropdown to explore the taxonomy of a single
#' sample, or all samples pooled together.
#'
#' Run:
#'   shiny::runApp()

library(shiny)
library(plotly)
library(dplyr)
library(tidyr)
library(readr)

source("R/utils.R")

feature_table <- read_csv("data/feature_table.csv", show_col_types = FALSE)
taxonomy <- parse_taxonomy(read_csv("data/taxonomy.csv", show_col_types = FALSE))

long_data <- feature_table %>%
  pivot_longer(-Feature_ID, names_to = "Sample", values_to = "Abundance") %>%
  left_join(taxonomy, by = "Feature_ID")

sample_choices <- c("All samples (pooled)", sort(unique(long_data$Sample)))

ui <- fluidPage(
  titlePanel("Interactive Taxonomy Explorer"),
  sidebarLayout(
    sidebarPanel(
      selectInput("sample", "Sample:", choices = sample_choices, selected = sample_choices[1]),
      helpText("Demo data — see README.md for details.")
    ),
    mainPanel(
      plotlyOutput("sunburst", height = "650px")
    )
  )
)

server <- function(input, output, session) {
  output$sunburst <- renderPlotly({
    df <- if (input$sample == "All samples (pooled)") {
      long_data
    } else {
      long_data %>% filter(Sample == input$sample)
    }

    family_totals <- df %>%
      group_by(Phylum, Class, Order, Family) %>%
      summarise(Abundance = sum(Abundance), .groups = "drop")

    sun_df <- build_sunburst_df(family_totals, c("Phylum", "Class", "Order", "Family"), "Abundance") %>%
      color_by_top_level()

    plot_ly(
      sun_df, ids = ~id, labels = ~label, parents = ~parent, values = ~value,
      type = "sunburst", branchvalues = "total",
      marker = list(colors = ~color),
      hovertemplate = "<b>%{label}</b><br>Reads: %{value}<br>%{percentParent} of parent<extra></extra>"
    )
  })
}

shinyApp(ui, server)
