# Portions of this script were adapted from Noble-Lab/FDRBench:
# https://github.com/Noble-Lab/FDRBench
# Licensed under the Apache License 2.0.
# Modifications were made for the DIA-NN analyses in this repository.

library(tidyverse)
library(fs)
library(ggpubr)

args <- commandArgs(trailingOnly = TRUE)
base_dir <- if (length(args) >= 1) {
  fs::path_expand(args[[1]])
} else {
  fs::path(getwd(), "output", "diann")
}
results_dir <- if (length(args) >= 2) {
  fs::path_expand(args[[2]])
} else {
  fs::path(getwd(), "results")
}
fs::dir_create(results_dir)

peptide_dir <- file.path(
  base_dir,
  "peptide_entrapment/accurate-design/fdp-calculation"
)

protein_dir <- file.path(
  base_dir,
  "protein_entrapment/accurate-design/fdp-calculation"
)

# Keep fdp at a constant rate of 1%
cutoff_grid <- tribble(
  ~fdp_cutoff, ~fdr_cutoff,
  0.01,        0.01,
  0.01,        0.02,
  0.01,        0.05,
  0.01,        0.10,
  # 0.05,        0.10
) %>%
  mutate(
    cutoff_label = paste0(
      "FDR ",
      scales::percent(fdr_cutoff, accuracy = 1),
      " / FDP ",
      scales::percent(fdp_cutoff, accuracy = 1)
    ),
    file_stub = paste0(
      "fdr_",
      str_replace(as.character(fdr_cutoff), "\\.", "p"),
      "_fdp_",
      str_replace(as.character(fdp_cutoff), "\\.", "p")
    )
  )

get_ratio <- function(fdp_file,
                      dataset_name = "",
                      level = "",
                      score_lower_is_better = TRUE,
                      fdr_cutoff = 0.1,
                      fdp_cutoff = 0.1,
                      fdr_decimal_place = 2) {

  a <- readr::read_csv(fdp_file, show_col_types = FALSE)

  required_cols <- c("q_value", "score", "paired_fdp")
  missing_cols <- setdiff(required_cols, names(a))

  if (length(missing_cols) > 0) {
    warning(
      "Skipping file because required columns are missing: ",
      fdp_file,
      "\nMissing: ",
      paste(missing_cols, collapse = ", ")
    )
    return(NULL)
  }

  a <- a %>%
    filter(!is.na(q_value), !is.na(score), !is.na(paired_fdp))

  if (nrow(a) == 0) {
    warning("Skipping empty FDP table after NA filtering: ", fdp_file)
    return(NULL)
  }

  # Original DIA-NN discoveries at 1% FDR
  n_original_discoveries <- a %>%
    filter(q_value <= fdr_cutoff) %>%
    nrow()

  if (n_original_discoveries == 0) {
    warning("No discoveries at ", fdr_cutoff, " FDR in: ", fdp_file)
    return(NULL)
  }

  # FDP at the DIA-NN 1% FDR threshold
  fdp_at_fdr <- a %>%
    filter(q_value <= fdr_cutoff) %>%
    arrange(desc(q_value)) %>%
    slice(1)

  fdp_label <- sprintf(
    paste0("%.", fdr_decimal_place, "f%%"),
    fdp_at_fdr$paired_fdp * 100
  )

  # Find the most permissive score threshold that still has paired FDP <= 1%
  passing_fdp <- a %>%
    filter(paired_fdp <= fdp_cutoff)

  if (nrow(passing_fdp) == 0) {
    warning("No discoveries pass paired FDP <= ", fdp_cutoff, " in: ", fdp_file)
    return(NULL)
  }

  if (score_lower_is_better) {
    max_score_fdp_cutoff <- max(passing_fdp$score, na.rm = TRUE)

    n_fdp_discoveries <- a %>%
      filter(score <= max_score_fdp_cutoff) %>%
      nrow()

  } else {
    max_score_fdp_cutoff <- min(passing_fdp$score, na.rm = TRUE)

    n_fdp_discoveries <- a %>%
      filter(score >= max_score_fdp_cutoff) %>%
      nrow()
  }

  if (n_fdp_discoveries == 0) {
    ratio <- NA_real_
  } else {
    # Inflation is expressed relative to the original FDR-based discoveries:
    # how many discoveries we'd "lose" or "gain" by using the FDP threshold
    # instead of the FDR threshold, with original discoveries as baseline.
    ratio <- (n_original_discoveries / n_fdp_discoveries) - 1
  }

  tibble(
  original_discoveries = n_original_discoveries,
  fdp_discoveries = n_fdp_discoveries,
  ratio = ratio,
  fdp_at_fdr = fdp_at_fdr$paired_fdp,
  fdp_at_fdr_label = fdp_label,
  observed_fdp_cutoff = max(passing_fdp$paired_fdp, na.rm = TRUE),
  score_cutoff = max_score_fdp_cutoff
  )
}

parse_run_name <- function(run_name) {
  instrument <- case_when(
    str_detect(run_name, "^astral") ~ "Astral",
    str_detect(run_name, "^qef") ~ "QExHF",
    TRUE ~ "Unknown"
  )

  library_type <- case_when(
    str_detect(run_name, "free") ~ "Library-free",
    str_detect(run_name, "prosit") ~ "Prosit",
    str_detect(run_name, "ms2pip") ~ "MS2PIP",
    str_detect(run_name, "carafe") ~ "CARAFE",
    TRUE ~ "Unknown"
  )

  # NOTE: matches both "random_shuffling" and "random-shuffling" so this
  # stays consistent with whatever separator is used in the filter step below.
  entrapment_type <- case_when(
    str_detect(run_name, "(?i)foreign[_-]c[_-]?elegans") ~ "Foreign C. elegans",
    str_detect(run_name, "(?i)random[_-]shuffling") ~ "Random shuffling",
    TRUE ~ "Unknown"
  )

  corrected <- str_detect(run_name, "corrected")

  tibble(
    instrument = instrument,
    library_type = library_type,
    entrapment_type = entrapment_type,
    corrected = corrected,
    dataset_name = paste(instrument, library_type, entrapment_type, sep = " - ")
  )
}

find_fdp_files <- function(root_dir, level) {
  dirs <- fs::dir_ls(root_dir, type = "directory")

  tibble(run_dir = dirs) %>%
    mutate(
      run_name = basename(run_dir),
      fdp_file = map_chr(run_dir, function(x) {
        candidates <- fs::dir_ls(
          x,
          regexp = paste0("diann_fdp_", level, "\\.csv$"),
          type = "file"
        )

        if (length(candidates) == 0) {
          candidates <- fs::dir_ls(
            x,
            regexp = paste0(".*fdp.*", level, ".*\\.csv$"),
            type = "file"
          )
        }

        if (length(candidates) == 0) {
          NA_character_
        } else {
          candidates[1]
        }
      })
    ) %>%
    filter(!is.na(fdp_file))
}

peptide_files <- find_fdp_files(peptide_dir, level = "precursor") %>%
  mutate(level = "peptide")

protein_files <- find_fdp_files(protein_dir, level = "protein") %>%
  mutate(level = "protein")

# Use the same separator-agnostic pattern here as in parse_run_name() above,
# so this filter and the entrapment_type detection always agree with each
# other regardless of whether folder names use "_" or "-".
all_files <- bind_rows(peptide_files, protein_files) %>%
  filter(str_detect(str_to_lower(run_name), "(?i)random[_-]shuffling")) %>%
  mutate(parsed = map(run_name, parse_run_name)) %>%
  unnest(parsed)

if (nrow(all_files) == 0) {
  stop(
    "all_files has 0 rows after filtering for 'random_shuffling' / ",
    "'random-shuffling'. Check that this pattern actually matches your run ",
    "directory names, e.g.:\n",
    "  fs::dir_ls(peptide_dir, type = 'directory') %>% basename() %>% str_subset('(?i)random')"
  )
}

summary_dat <- all_files %>%
  tidyr::crossing(cutoff_grid) %>%
  mutate(
    result = pmap(
      list(fdp_file, fdr_cutoff, fdp_cutoff),
      function(fdp_file, fdr_cutoff, fdp_cutoff) {
        get_ratio(
          fdp_file = fdp_file,
          score_lower_is_better = TRUE,
          fdr_cutoff = fdr_cutoff,
          fdp_cutoff = fdp_cutoff,
          fdr_decimal_place = 2
        )
      }
    )
  ) %>%
  filter(!map_lgl(result, is.null)) %>%
  unnest(result)

summary_dat %>%
  write_tsv(fs::path(results_dir, "diann_fdp_vs_fdr_all_cutoff_combinations.tsv"))

# If you want one representative run per (level, instrument, library_type,
# entrapment_type) group instead of plotting every individual run, use
# plot_dat_collapsed below. Otherwise plot_dat (all runs with a valid ratio)
# is used for the figure.
plot_dat_collapsed <- summary_dat %>%
  group_by(level, instrument, library_type, entrapment_type) %>%
  arrange(desc(corrected), desc(str_detect(run_name, "_v2")), run_name) %>%
  slice(1) %>%
  ungroup()

plot_dat <- summary_dat %>%
  filter(!is.na(ratio)) %>%
  mutate(
    level = factor(level, levels = c("peptide", "protein")),
    cutoff_label = factor(
      cutoff_label,
      levels = cutoff_grid$cutoff_label
    )
  )

plot_dat %>%
  write_tsv(fs::path(results_dir, "diann_fdp_vs_fdr_all_cutoff_combinations_plotting.tsv"))

# make_inflation_plot <- function(dat, level_label = NULL) {
#   if (!is.null(level_label)) {
#     dat <- dat %>% filter(level == level_label)
#   }

#   dat <- dat %>%
#     mutate(
#       instrument = factor(
#         instrument,
#         levels = c("Astral", "QExHF")
#       ),
#       library_type = factor(
#         library_type,
#         levels = c("Library-free", "Prosit", "MS2PIP", "CARAFE")
#       )
#     ) %>%
#     arrange(library_type, instrument) %>%
#     mutate(
#       dataset_name = factor(dataset_name, levels = unique(dataset_name))
#     )

#   ggplot(dat, aes(x = dataset_name, y = ratio, fill = library_type)) +
#     geom_col(width = 0.65) +
#     geom_text(
#       aes(
#         label = paste0(
#           fdp_at_1pct_fdr,
#           "\n",
#           scales::percent(ratio, accuracy = 0.1),
#           "\n",
#           original_discoveries,
#           " \u2192 ",
#           fdp_discoveries
#         )
#       ),
#       hjust = 0,
#       size = 6,
#       lineheight = 0.9
#     ) +
#     facet_grid(level ~ ., scales = "free_y", space = "free_y") +
#     coord_flip() +
#     scale_y_continuous(
#       labels = scales::percent
#     ) +
#     theme_pubr(base_size = 14, border = TRUE) +
#     theme(
#       axis.text.x = element_text(angle = 0, hjust = 0.5),
#       axis.text.y = element_text(angle = 0, hjust = 1),
#       legend.position = "bottom",
#       plot.margin = unit(c(0.1, 0.5, 0.2, 0.1), "cm"),
#       strip.text = element_blank(),
#       strip.background = element_blank()
#     ) +
#     labs(
#       x = "Dataset",
#       y = "Discovery inflation at 1% paired FDP",
#       fill = "Library type",
#       title = level_label
#     )
# }

make_inflation_plot <- function(dat, level_label = NULL) {
  if (!is.null(level_label)) {
    dat <- dat %>% filter(level == level_label)
  }

  dat <- dat %>%
    mutate(
      instrument = factor(instrument, levels = c("Astral", "QExHF")),
      library_type = factor(
        library_type,
        levels = c("Library-free", "Prosit", "MS2PIP", "CARAFE")
      )
    ) %>%
    arrange(cutoff_label, library_type, instrument) %>%
    mutate(
      dataset_name = factor(dataset_name, levels = rev(unique(dataset_name)))
    )

  ggplot(dat, aes(x = dataset_name, y = ratio, fill = library_type)) +
    geom_col(width = 0.65) +
    geom_text(
      aes(
        label = paste0(
          fdp_at_fdr_label,
          "\n",
          scales::percent(ratio, accuracy = 0.1),
          "\n",
          original_discoveries,
          " -> ",
          fdp_discoveries
        )
      ),
      hjust = 0,
      size = 3,
      lineheight = 0.9
    ) +
    facet_grid(cutoff_label ~ ., scales = "free_y", space = "free_y") +
    coord_flip() +
    scale_y_continuous(
      labels = scales::percent,
      expand = expansion(mult = c(0, 0.25))
    ) +
    theme_pubr(base_size = 12, border = TRUE) +
    theme(
      axis.text.x = element_text(angle = 0, hjust = 0.5),
      axis.text.y = element_text(angle = 0, hjust = 1),
      legend.position = "bottom",
      plot.margin = unit(c(0.1, 0.5, 0.2, 0.1), "cm"),
      strip.text = element_text(size = 11),
      strip.background = element_rect(fill = "grey95")
    ) +
    labs(
      x = "Dataset",
      y = "Discovery inflation",
      fill = "Library type",
      title = level_label
    )
}

levels <- unique(plot_dat$level)

plots_by_level <- set_names(levels) %>%
  map(~ make_inflation_plot(plot_dat, .x))

# Print + save one PDF/SVG per level
walk2(plots_by_level, names(plots_by_level), function(gg, level_label) {
  print(gg)

  file_stub <- str_replace_all(str_to_lower(level_label), "[^a-z0-9]+", "-")

  ggsave(
    fs::path(results_dir, sprintf("Figure-diann_fdp_fdr_all_cutoff_combinations-%s.pdf", file_stub)),
    gg,
    width = 9,
    height = 14
  )

  ggsave(
    fs::path(results_dir, sprintf("Figure-diann_fdp_fdr_all_cutoff_combinations-%s.svg", file_stub)),
    gg,
    width = 9,
    height = 14
  )

})

# Plot all inflation rates in the same plot wihout facetting
plot_dat_stacked <- plot_dat %>%
  filter(!is.na(ratio)) %>%
  mutate(
    level = factor(level, levels = c("peptide", "protein")),
    instrument = factor(instrument, levels = c("Astral", "QExHF")),
    library_type = factor(
      library_type,
      levels = c("Library-free", "Prosit", "MS2PIP", "CARAFE")
    ),
    cutoff_label = factor(
      cutoff_label,
      levels = c(
        "FDR 1% / FDP 1%",
        "FDR 2% / FDP 1%",
        "FDR 5% / FDP 1%",
        "FDR 10% / FDP 1%"
      )
    ),
    dataset_label = paste(instrument, library_type, sep = " - ")
  ) %>%
  arrange(library_type, instrument, cutoff_label) %>%
  mutate(
    dataset_label = factor(dataset_label, levels = rev(unique(dataset_label)))
  )

cutoff_colours <- c(
  "FDR 1% / FDP 1%" = "#2166AC",
  "FDR 2% / FDP 1%" = "#67A9CF",
  "FDR 5% / FDP 1%" = "#F4A582",
  "FDR 10% / FDP 1%" = "#B2182B"
)

make_stacked_diverging_plot <- function(dat, level_label = NULL) {
  if (!is.null(level_label)) {
    dat <- dat %>% filter(level == level_label)
  }

  ggplot(
    dat,
    aes(
      x = dataset_label,
      y = ratio,
      fill = cutoff_label
    )
  ) +
    geom_hline(
      yintercept = 0,
      linewidth = 0.4,
      colour = "grey30"
    ) +
    geom_col(
      width = 0.65,
      position = position_dodge(width = 0.75)
    ) +
    coord_flip() +
    facet_grid(level ~ ., scales = "free_y", space = "free_y") +
    scale_fill_manual(values = cutoff_colours) +
    scale_y_continuous(
      labels = scales::percent,
      expand = expansion(mult = c(0.12, 0.12))
    ) +
    theme_pubr(base_size = 12, border = TRUE) +
    theme(
      legend.position = "bottom",
      axis.text.x = element_text(angle = 0, hjust = 0.5),
      axis.text.y = element_text(angle = 0, hjust = 1),
      plot.margin = unit(c(0.1, 0.5, 0.2, 0.1), "cm"),
      strip.text = element_blank(),
      strip.background = element_blank()
    ) +
    labs(
      x = "Dataset",
      y = "Cumulative discovery inflation",
      fill = "Cutoff combination",
      title = level_label
    )
}

levels_to_plot <- unique(plot_dat$level)

plots_by_level <- set_names(levels_to_plot) %>%
  map(~ make_stacked_diverging_plot(plot_dat_stacked, .x))

walk2(plots_by_level, names(plots_by_level), function(gg, level_label) {
  print(gg)

  file_stub <- str_replace_all(str_to_lower(level_label), "[^a-z0-9]+", "-")

  ggsave(
    file.path(
      results_dir,
      sprintf("Figure-stacked-diverging-fdr-fdp-cutoffs-%s.pdf", file_stub)
    ),
    gg,
    width = 8,
    height = 8
  )

  ggsave(
    file.path(
      results_dir,
      sprintf("Figure-stacked-diverging-fdr-fdp-cutoffs-%s.svg", file_stub)
    ),
    gg,
    width = 8,
    height = 8
  )
})

# Convenience handles if you want them individually, e.g.:
# gg_astral <- plots_by_level[["Astral"]]
# gg_qexhf  <- plots_by_level[["QExHF"]]

# All in one plot
# gg <- ggplot(plot_dat, aes(x = dataset_name, y = ratio, fill = library_type)) +
#   geom_col(width = 0.65) +
#   geom_text(
#     aes(
#       label = paste0(
#         fdp_at_1pct_fdr,
#         "\n",
#         scales::percent(ratio, accuracy = 0.1),
#         "\n",
#         original_discoveries,
#         " \u2192 ",
#         fdp_discoveries
#       )
#     ),
#     hjust = -0.05,
#     size = 3,
#     lineheight = 0.9
#   ) +
#   facet_grid(level ~ instrument, scales = "free_y", space = "free_y") +
#   coord_flip() +
#   scale_y_continuous(
#     labels = scales::percent,
#     expand = expansion(mult = c(0, 0.25))
#   ) +
#   theme_pubr(base_size = 12, border = TRUE) +
#   theme(
#     axis.text.x = element_text(angle = 0, hjust = 0.5),
#     axis.text.y = element_text(angle = 0, hjust = 1),
#     legend.position = "bottom",
#     plot.margin = unit(c(0.1, 0.5, 0.2, 0.1), "cm"),
#     strip.text = element_text(size = 12)
#   ) +
#   labs(
#     x = "Dataset",
#     y = "Discovery inflation at 1% paired FDP",
#     fill = "Library type"
#   )

# print(gg)

# dir.create("output", showWarnings = FALSE)

# ggsave(
#   "output/Figure-diann_fdp-fdr-comparison-all-runs.pdf",
#   gg,
#   width = 12,
#   height = 12
# )

# ggsave(
#   "output/Figure-diann_fdp-fdr-comparison-all-runs.svg",
#   gg,
#   width = 12,
#   height = 12
# )
