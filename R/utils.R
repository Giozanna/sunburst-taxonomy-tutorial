#' utils.R — shared helpers for the taxonomy sunburst tutorial

library(dplyr)
library(tidyr)
library(stringr)

#' Parse a QIIME2/SILVA-style semicolon taxonomy string into separate rank
#' columns, e.g. "d__Fungi;p__Ascomycota;c__Sordariomycetes;o__Hypocreales;f__Nectriaceae"
#' becomes Kingdom = "Fungi", Phylum = "Ascomycota", ...
#'
#' @param taxonomy_df a data frame with a `Taxon` column
#' @return the same data frame with Kingdom/Phylum/Class/Order/Family columns added
parse_taxonomy <- function(taxonomy_df) {
  taxonomy_df %>%
    tidyr::separate(
      Taxon,
      into = c("Kingdom", "Phylum", "Class", "Order", "Family"),
      sep = ";",
      fill = "right",
      extra = "drop"
    ) %>%
    dplyr::mutate(dplyr::across(
      Kingdom:Family,
      ~ stringr::str_trim(stringr::str_remove(.x, "^[a-z]__"))
    ))
}

#' Build the ids/labels/parents/values table a plotly sunburst needs, from a
#' data frame that has one row per LEAF taxon (e.g. Family) with a numeric
#' abundance column. Aggregates a value for every ancestor node too, since a
#' sunburst needs a value at *every* level of the hierarchy, not just the
#' leaves.
#'
#' @param df data frame with the rank columns + a value column
#' @param rank_cols character vector of column names, in hierarchical order
#'   e.g. c("Phylum", "Class", "Order", "Family")
#' @param value_col name of the numeric column to aggregate (e.g. "Abundance")
#' @return data frame with columns id, label, parent, value — ready for
#'   plot_ly(..., ids = ~id, labels = ~label, parents = ~parent, values = ~value)
build_sunburst_df <- function(df, rank_cols, value_col) {
  levels_list <- vector("list", length(rank_cols))

  for (i in seq_along(rank_cols)) {
    cols_i <- rank_cols[1:i]

    agg <- df %>%
      dplyr::group_by(dplyr::across(dplyr::all_of(cols_i))) %>%
      dplyr::summarise(value = sum(.data[[value_col]]), .groups = "drop") %>%
      tidyr::unite("id", dplyr::all_of(cols_i), sep = "/", remove = FALSE)

    if (i == 1) {
      agg$parent <- ""
    } else {
      parent_cols <- rank_cols[1:(i - 1)]
      agg <- agg %>%
        tidyr::unite("parent", dplyr::all_of(parent_cols), sep = "/", remove = FALSE)
    }

    agg$label <- agg[[cols_i[i]]]
    levels_list[[i]] <- agg %>% dplyr::select(id, label, parent, value)
  }

  dplyr::bind_rows(levels_list) %>% dplyr::distinct(id, .keep_all = TRUE)
}

#' Fixed color per phylum, so the same phylum always gets the same color
#' across different sunbursts (e.g. when filtering by sample in the Shiny app).
phylum_palette <- function() {
  c(
    "Ascomycota"    = "#4C72B0",
    "Basidiomycota" = "#C44E52",
    "Mucoromycota"  = "#55A868"
  )
}

#' Assign each node in a sunburst data frame the color of its top-level
#' (Phylum) ancestor, so classes/orders/families inherit their phylum's hue.
color_by_top_level <- function(sun_df, palette = phylum_palette()) {
  sun_df %>%
    dplyr::mutate(
      top_level = stringr::word(id, 1, sep = stringr::fixed("/")),
      color = unname(palette[top_level])
    )
}
