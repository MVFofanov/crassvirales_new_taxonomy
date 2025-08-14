#!/usr/bin/env Rscript

# ---- deps (simple) ----
pkgs_cran <- c("ape","ggplot2","dplyr","readr","tibble","tidyr","stringr")
for (p in pkgs_cran) if (!requireNamespace(p, quietly=TRUE)) install.packages(p, repos="https://cloud.r-project.org")
if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager", repos="https://cloud.r-project.org")
for (p in c("ggtree","treeio")) if (!requireNamespace(p, quietly=TRUE)) BiocManager::install(p, ask=FALSE, update=FALSE)

library(ape)
library(ggplot2)
library(ggtree)
library(treeio)
library(dplyr)
library(readr)
library(tibble)
library(tidyr)
library(stringr)

# ---- palettes ----
DEFAULT_COLOR <- "#bfbfbf"
  crass_family_colors <- c(
    "Intestiviridae"="#EE3B3B","Crevaviridae"="#EE9A00","Suoliviridae"="#4169E1",
    "Steigviridae"="#00CED1","Epsilon"="#CD2990","Zeta"="#006400"
  )
  non_cv_colors <- c("Unknown"="#BDBDBD","outgroup"="#808080","NA"="#BDBDBD")
  family_palette <- c(crass_family_colors, non_cv_colors)
  origin_palette <- c("MAG"="#1f77b4","Isolate"="#d62728")
  NOT_CV_TAGS    <- c("outgroup","NA")
  
  # ---- paths ----
  proj_dir  <- "C:/crassvirales/crassvirales_new_taxonomy/crassvirales_prophages/blast_prophages_vs_ncbi_and_gtdb"
  crassus_results_dir <- "C:/crassvirales/CrassUS_old/CrassUS/results/prophage_analysis"
  
  tree_file <- file.path(crassus_results_dir, "5_phylogenies/3_iToL/Terl_iToL_renamed.nwk")
  ann_file  <- file.path(proj_dir, "Crassphage_prophage_analysis_annotation.tsv")
  itol_file <- file.path(crassus_results_dir, "5_phylogenies/3_iToL/TerL_iToL_simplified.txt")
  itol_simplified_file <- file.path(dirname(itol_file), "TerL_iToL_simplified.txt")
  
  out_png <- file.path(proj_dir, "TerL_Crassvirales_pruned_collapsed.png")
  out_svg <- file.path(proj_dir, "TerL_Crassvirales_pruned_collapsed.svg")
  out_pdf <- file.path(proj_dir, "TerL_Crassvirales_pruned_collapsed.pdf")
  
  cat("Tree file: ", normalizePath(tree_file, winslash = "/"), "\n")
  cat("Prophage annotation: ", normalizePath(ann_file, winslash = "/"), "\n")
  cat("iTOL simplified: ", normalizePath(itol_file, winslash = "/"), "\n")
  
  # ---- 1) load ----
  tr <- read.tree(tree_file)
  
  # Use ONLY genome_origin from the annotation file
  ann <- suppressMessages(readr::read_tsv(ann_file, show_col_types = FALSE)) %>%
    transmute(
      label           = as.character(protein_id),
      prophage_family = as.character(crassvirales_family),
      origin          = ifelse(genome_origin %in% c("MAG","Isolate"), genome_origin, NA_character_)
    )
  
  # iTOL mapping (for non-prophage references)
  itol <- suppressMessages(readr::read_tsv(itol_simplified_file, show_col_types = FALSE)) %>%
    transmute(
      label       = as.character(genome_id),
      itol_family = as.character(crassvirales_family)
    )
  
  # ---- 2) tip metadata ----
  tips_tbl <- tibble(label = tr$tip.label) %>%
    left_join(itol, by="label") %>%
    left_join(ann,  by="label") %>%
    mutate(
      is_prophage     = !is.na(prophage_family),
      family          = dplyr::coalesce(prophage_family, itol_family),
      family          = ifelse(family %in% names(family_palette), family, "Unknown"),
      in_crassvirales = !is.na(family) & !(family %in% NOT_CV_TAGS) & family != "Unknown"
    )
  
  # ---- 3) prune to MRCA of all Crassvirales leaves ----
  cv_tips <- tips_tbl %>% filter(in_crassvirales) %>% pull(label) %>% unique()
  if (length(cv_tips) < 2) stop("Not enough Crassvirales tips to define an MRCA (need >= 2).")
  
  mrca_node <- ape::getMRCA(tr, cv_tips)
  if (is.null(mrca_node) || is.na(mrca_node)) stop("MRCA not found for the provided Crassvirales tips.")
  
  cv_clade          <- ape::extract.clade(tr, mrca_node)
  cv_tip_labels     <- cv_clade$tip.label
  tr_pruned         <- ape::keep.tip(tr, cv_tip_labels)
  
  tips_tbl_pruned <- tips_tbl %>%
    filter(label %in% tr_pruned$tip.label) %>%
    mutate(family = factor(family, levels = names(family_palette)))  # for stable legend
  
  # ---- 4) select “pure non-prophage” nodes to collapse ----
  internal_nodes <- (Ntip(tr_pruned) + 1):(Ntip(tr_pruned) + tr_pruned$Nnode)
  desc_tip_count <- desc_prophage_count <- integer(length(internal_nodes))
  for (i in seq_along(internal_nodes)) {
    nd <- internal_nodes[i]
    tips_under <- ape::extract.clade(tr_pruned, nd)$tip.label
    desc_tip_count[i]      <- length(tips_under)
    desc_prophage_count[i] <- sum(tips_tbl_pruned$is_prophage[match(tips_under, tips_tbl_pruned$label)], na.rm = TRUE)
  }
  collapse_nodes <- internal_nodes[desc_tip_count > 1 & desc_prophage_count == 0]
  collapsed_info <- tibble(node = collapse_nodes,
                           collapsed_n = desc_tip_count[match(collapse_nodes, internal_nodes)])
  
  message(sprintf("Will collapse %d clades (pure non-prophage); total collapsed tips across them = %d",
                  length(collapse_nodes), sum(collapsed_info$collapsed_n)))
  
  # ---- 5) side table (ONLY for prophages) ----
  side_df <- tips_tbl_pruned %>%
    transmute(
      label,
      `Crassvirales family` = ifelse(is_prophage, as.character(family), NA_character_),
      `Genome origin`       = ifelse(is_prophage, origin,            NA_character_)
    ) %>%
    # keep EXACT tree tip order to avoid misalignment
    right_join(tibble(label = tr_pruned$tip.label), by = "label") %>%
    arrange(match(label, tr_pruned$tip.label)) %>%
    tibble::column_to_rownames("label")
  
  # Build the two data.frames gheatmap expects (1 column each)
  family_df <- side_df[, "Crassvirales family", drop = FALSE]
  origin_df <- side_df[, "Genome origin",       drop = FALSE]
  
  # ---- 6) layout numbers so labels never overlap panels ----
  # measure tree span *after* collapse to compute offsets in data units
  p0 <- ggtree(tr_pruned, layout = "rectangular") %<+% tips_tbl_pruned
  for (nd in collapse_nodes) p0 <- collapse(p0, node = nd)
  x_span <- diff(range(p0$data$x, na.rm = TRUE))
  
  # Keep labels close; keep panels narrow
  lab_offset    <- 0.03 * x_span
  gap_after_lab <- 0.012 * x_span
  width_family  <- 0.06 * x_span
  gap_between   <- 0.015 * x_span
  width_origin  <- 0.04 * x_span
  
  # Push both panels further to the right (tweak this to taste)
  panel_shift   <- 0.06 * x_span
  
  # Final offsets (do NOT call xlim_tree yet)
  offset_family <- lab_offset + gap_after_lab + panel_shift
  offset_origin <- offset_family + width_family + gap_between
  extra_right   <- 0.10 * x_span
  
  # ---- 7) plot ----
  p <- ggtree(tr_pruned, layout = "rectangular") %<+% tips_tbl_pruned +
    theme_tree2() +
    geom_tippoint(aes(color = family), size = 1.2, alpha = 0.9) +
    ggtitle("TerL — Crassvirales clade (pruned), non-prophage clades collapsed",
            subtitle = "Triangles = collapsed clades; numbers = # tips inside")
  
  # collapse and add triangle counts
  for (nd in collapse_nodes) p <- collapse(p, node = nd)
  pd <- p$data
  label_pos <- pd %>% select(node, x, y) %>% inner_join(collapsed_info, by = "node")
  p <- p + geom_label(data = label_pos, aes(x = x, y = y, label = collapsed_n),
                      size = 3, label.size = 0.2, label.padding = unit(0.12, "lines"),
                      fill = "white")
  
  # ring prophage tips (no legend)
  p <- p + geom_tippoint(
    data = subset(pd, isTip & !is.na(is_prophage) & is_prophage),
    shape = 21, stroke = 0.2, size = 2, fill = NA, color = "black",
    show.legend = FALSE
  )
  
  # aligned tip labels, kept close to tips
  p <- p + geom_tiplab(size = 2, align = TRUE, offset = lab_offset, linesize = 0)
  
  # Reserve canvas space ONCE (include panel_shift and both panel widths)
  p <- p + xlim_tree(
    max(p$data$x, na.rm = TRUE) +
      lab_offset + gap_after_lab + panel_shift +
      width_family + gap_between + width_origin +
      extra_right
  )
  
  # ---- SIDE PANELS (simple, non-overlapping, scalable) ----
  
  # Build side tables (ONLY for prophages), keep exact tip order
  side_df <- tips_tbl_pruned %>%
    transmute(
      label,
      `Crassvirales family` = ifelse(is_prophage, as.character(family), NA_character_),
      `Genome origin`       = ifelse(is_prophage, origin,            NA_character_)
    ) %>%
    right_join(tibble(label = tr_pruned$tip.label), by = "label") %>%
    arrange(match(label, tr_pruned$tip.label)) %>%
    tibble::column_to_rownames("label")
  
  # One column per gheatmap call
  family_df <- side_df[, "Crassvirales family", drop = FALSE]
  origin_df <- side_df[, "Genome origin",       drop = FALSE]
  
  # Measure span AFTER collapsing to set widths/offsets in data units
  p0 <- ggtree(tr_pruned, layout = "rectangular") %<+% tips_tbl_pruned
  for (nd in collapse_nodes) p0 <- collapse(p0, node = nd)
  x_span <- diff(range(p0$data$x, na.rm = TRUE))
  
  # Keep labels close to tips
  lab_offset  <- 0.025 * x_span
  
  # Panel sizes + consistent gap (tweak to taste)
  panel_widths <- c(family = 0.045 * x_span,
                    origin = 0.035 * x_span)
  panel_gap    <- 0.012 * x_span
  
  # Push the whole block of panels to the right of labels
  panel_shift  <- 0.06  * x_span
  
  # Compute offsets for each panel (no overlap)
  panel_names  <- names(panel_widths)
  cum_left     <- c(0, cumsum(head(panel_widths, -1) + panel_gap))
  panel_offsets <- lab_offset + panel_shift + cum_left
  names(panel_offsets) <- panel_names
  
  # Draw panels in order
  p <- gheatmap(
    p, family_df,
    offset   = panel_offsets["family"],
    width    = panel_widths["family"],
    colnames = TRUE, colnames_position = "top",
    font.size = 3, hjust = 0, color = NA
  )
  
  p <- gheatmap(
    p, origin_df,
    offset   = panel_offsets["origin"],
    width    = panel_widths["origin"],
    colnames = TRUE, colnames_position = "top",
    font.size = 3, hjust = 0, color = NA
  )
  
  # One simple legend for both panels
  p <- p + scale_fill_manual(
    name     = "Side panels",
    values   = c(family_palette, origin_palette),
    na.value = "#FFFFFF",   # blank for non-prophage rows
    drop     = FALSE
  )
  
  # Reserve canvas space ONCE (labels + shift + widths + gaps)
  total_panels_width <- sum(panel_widths) + panel_gap * (length(panel_widths) - 1)
  extra_right        <- 0.10 * x_span
  
  p <- p + xlim_tree(
    max(p$data$x, na.rm = TRUE) +
      lab_offset + panel_shift + total_panels_width + extra_right
  )
  
  
  # ---- save ----
  ggsave(out_png, p, width = 18, height = 14, dpi = 300)
  ggsave(out_svg, p, width = 18, height = 14)
  ggsave(out_pdf, p, width = 18, height = 14)
  message(sprintf("Saved plots:\n  %s\n  %s", out_png, out_pdf))
  