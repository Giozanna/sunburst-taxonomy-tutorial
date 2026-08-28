#' 00_generate_demo_data.R
#'
#' Generates a SYNTHETIC demo dataset shaped like a typical QIIME2 /
#' amplicon-sequencing output:
#'   - data/feature_table.csv : Feature_ID (ASV) x Sample read counts
#'   - data/taxonomy.csv      : Feature_ID -> semicolon-delimited taxonomy
#'                               string (SILVA/QIIME2 style: d__;p__;c__;o__;f__)
#'                               + a confidence score, as exported by
#'                               `qiime taxa` / DADA2 assignTaxonomy()
#'
#' Taxa are real fungal Phylum/Class/Order/Family names (standard reference
#' taxonomy); the READ COUNTS are simulated. This is demo data for the
#' tutorial, not a real sequencing run.
#'
#' Run:
#'   Rscript R/00_generate_demo_data.R

library(dplyr)
library(tibble)
library(readr)

set.seed(42)

# Phylum > Class > Order > Family
tree <- list(
  Ascomycota = list(
    Sordariomycetes = list(
      Hypocreales = c("Nectriaceae", "Hypocreaceae"),
      Xylariales  = c("Xylariaceae", "Diatrypaceae")
    ),
    Dothideomycetes = list(
      Pleosporales       = c("Didymellaceae", "Pleosporaceae"),
      Botryosphaeriales  = c("Botryosphaeriaceae", "Phyllostictaceae")
    )
  ),
  Basidiomycota = list(
    Agaricomycetes = list(
      Agaricales  = c("Marasmiaceae", "Psathyrellaceae"),
      Polyporales = c("Polyporaceae", "Meruliaceae")
    ),
    Tremellomycetes = list(
      Tremellales = c("Tremellaceae", "Bulleraceae")
    )
  ),
  Mucoromycota = list(
    Mucoromycetes = list(
      Mucorales = c("Mucoraceae", "Cunninghamellaceae")
    )
  )
)

# flatten the nested list into a Phylum/Class/Order/Family lookup table
# (plain nested loops, deliberately avoiding purrr's map*_dfr — deprecated
# in recent purrr versions)
rows <- list()
for (phylum in names(tree)) {
  for (klass in names(tree[[phylum]])) {
    for (order in names(tree[[phylum]][[klass]])) {
      families <- tree[[phylum]][[klass]][[order]]
      rows[[length(rows) + 1]] <- tibble(
        Phylum = phylum, Class = klass, Order = order, Family = families
      )
    }
  }
}
tax_levels <- bind_rows(rows)

cat(nrow(tax_levels), "leaf families across", n_distinct(tax_levels$Phylum), "phyla\n")

# assign 60 ASVs to families (weighted, so some families are more diverse than others)
n_features <- 60
family_weights <- rgamma(nrow(tax_levels), shape = 2, rate = 1)
family_idx <- sample(seq_len(nrow(tax_levels)), n_features, replace = TRUE, prob = family_weights)

feature_family <- tax_levels[family_idx, ] %>%
  mutate(Feature_ID = sprintf("ASV_%03d", row_number())) %>%
  relocate(Feature_ID)

taxonomy_df <- feature_family %>%
  mutate(
    Taxon = paste0("d__Fungi;p__", Phylum, ";c__", Class, ";o__", Order, ";f__", Family),
    Confidence = round(runif(n(), 0.85, 0.999), 4)
  ) %>%
  select(Feature_ID, Taxon, Confidence)

# per-feature baseline abundance (log-normal), Poisson counts across 8 samples
samples <- paste0("Sample_", 1:8)
baseline <- rlnorm(n_features, meanlog = 3.0, sdlog = 1.1)

counts <- sapply(seq_along(samples), function(i) {
  rpois(n_features, lambda = baseline * runif(n_features, 0.6, 1.4))
})
colnames(counts) <- samples

feature_table <- as_tibble(counts) %>%
  mutate(Feature_ID = feature_family$Feature_ID) %>%
  relocate(Feature_ID)

dir.create("data", showWarnings = FALSE)
write_csv(taxonomy_df, "data/taxonomy.csv")
write_csv(feature_table, "data/feature_table.csv")

cat("Total simulated reads:", sum(counts), "\n")
cat("Wrote data/taxonomy.csv and data/feature_table.csv\n")
