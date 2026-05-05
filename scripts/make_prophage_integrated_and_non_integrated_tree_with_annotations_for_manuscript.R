library(ggtree)
library(ggplot2)
library(ape)
library(dplyr)
library(readr)
library(tibble)
library(ggnewscale)

# =========================================================
# Paths and files
# =========================================================

# Main working directory
#WD <- "C:/crassvirales/crassvirales_new_taxonomy/crassvirales_prophages/mafft_iqtree_TerL_analysis_reference_best_TerLs_without_KC821624_1_100418/iqtree_trees"
# WD <- "C:/crassvirales/crassvirales_new_taxonomy/crassvirales_prophages/tree_plots_for_manuscript"
WD <- "C:/crassvirales/crassvirales_new_taxonomy/crassvirales_prophages/tree_plots_for_manuscript/MRCA_Crassvirales_iqtree_trees"
# WD <- "C:/crassvirales/crassvirales_new_taxonomy/crassvirales_prophages/tree_plots_for_manuscript/MRCA_integrated_prophage_candidates_iqtree_trees"

DATA_DIR <- file.path(WD, "data")

# Output directory
OUTPUT_DIR <- file.path(WD, "figures")

# Shared input files
SOURCE_FILE <- file.path(DATA_DIR, "all_genomad_contigs_mags_or_isolates.txt")
PAIRS_TSV   <- file.path(DATA_DIR, "pairs_prophage_with_flanks_vs_bacterial.tsv")
CHECKV_FILE <- file.path(DATA_DIR, "checkv_quality_summary.tsv")
CANDIDATE_LIST_FILE <- file.path(DATA_DIR, "Crassvirales_integrated_prophage_candidate_list.txt")

# # Tree files
# TREE_FILES <- c(
#   file.path(DATA_DIR, "TerL_gappy0.7.treefile")
# )

TREE_FILES <- c(
  file.path(DATA_DIR, "TerL_gappy0.7.treefile"),
  file.path(DATA_DIR, "TerL_gappy0.8.treefile"),
  file.path(DATA_DIR, "TerL_gappy0.9.treefile"),
  file.path(DATA_DIR, "TerL_kpic.treefile"),
  file.path(DATA_DIR, "TerL_smart-gap.treefile"),
  file.path(DATA_DIR, "TerL_untrimmed.treefile")
)


# Create output directory if needed
if (!dir.exists(OUTPUT_DIR)) {
  dir.create(OUTPUT_DIR, recursive = TRUE)
}

# =========================================================
# Settings
# =========================================================

# ---- Ring spacing settings ----
# ---- Ring spacing settings ----
BAR_LINEWIDTH <- 0.3

PROPHAGE_LINEWIDTH     <- 0.25
PROPHAGE_RATIO_LINEWIDTH     <- 0.25
PHYLUM_LINEWIDTH       <- BAR_LINEWIDTH
CLASS_LINEWIDTH        <- BAR_LINEWIDTH
SOURCE_LINEWIDTH       <- BAR_LINEWIDTH
COMPLETENESS_LINEWIDTH <- BAR_LINEWIDTH

TREE_TO_GENOMAD_GAP <- 0.005 #0.03
#TREE_TO_GENOMAD_GAP <- 0.01 #0.03
GENOMAD_TO_RATIO_GAP <- 0.03
RATIO_TO_PHYLUM_GAP <- 0.03
PHYLUM_TO_CLASS_GAP <- 0.02
CLASS_TO_DOT_GAP <- 0.03
TEXT_EXTRA_GAP <- 0.03

GENOMAD_DIV <- 500000
GENOMAD_LEVELS <- c(50000, 100000, 150000)

#RATIO_LEVELS <- c(0.25, 0.50, 0.75, 1.00)
RATIO_LEVELS <- c(0.50, 1.00)
RATIO_MULT <- max(GENOMAD_LEVELS) / GENOMAD_DIV

# 50% of max prophage-size bar width (= 0.3 / 2 = 0.15)
PHYLUM_RING_WIDTH <- 0.15
CLASS_RING_WIDTH  <- 0.15

DOT_SIZE    <- 0.4
DOT_STROKE  <- 0.15

SOURCE_RING_WIDTH <- CLASS_RING_WIDTH
SOURCE_RING_GAP <- CLASS_TO_DOT_GAP

SHORT_ID_TO_LABEL_GAP <- 0.02
SHORT_ID_TEXT_SIZE <- 0.3
SHORT_ID_TEXT_SIZE_COLLAPSED <- 1
SHORT_ID_TEXT_COLOR <- "black"

SOURCE_COLORS <- c(
  "isolate" = "black",
  "MAG"     = "grey80"
)

SOURCE_TO_COMPLETENESS_GAP <- 0.02
COMPLETENESS_RING_WIDTH <- 0.15

# ROOTING_MODES <- c(
#   "outgroup",
#   "crass_mrca",
#   "crass_plus_integrated_mrca"
# )

ROOTING_MODES <- c(
  "outgroup"
)

CRASSVIRALES_COLOR_SCHEME <- c(
  "Intestiviridae" = "#EE3B3B",
  "Crevaviridae"   = "#EE9A00",
  "Suoliviridae"   = "#4169E1",
  "Steigviridae"   = "#00CED1",
  "Epsilon"        = "#CD2990",
  "Zeta"           = "#006400",
  "Outgroup"       = "violet",
  "Other"          = "black"
  #"Other"          = "grey60"
)

CRASS_CLADE_FILL <- c(
  "Intestiviridae" = scales::alpha("#EE3B3B", 0.2),
  "Crevaviridae"   = scales::alpha("#EE9A00", 0.2),
  "Suoliviridae"   = scales::alpha("#4169E1", 0.2),
  "Steigviridae"   = scales::alpha("#00CED1", 0.2),
  "Epsilon"        = scales::alpha("#CD2990", 0.2),
  "Zeta"           = scales::alpha("#006400", 0.2)
)

# ---- Phylum color scheme ----
PHYLUM_COLORS <- c(
  "Bacteroidota"   = "#33a02c",
  "Bacillota"      = "#1f78b4",
  "Pseudomonadota" = "#ff7f00",
  "Actinomycetota" = "#6A3D9A",
  "Bacteria"       = "gold",
  "Other"          = "grey70"
)

# ---- Class color scheme ----
CLASS_COLORS <- c(
  "Flavobacteriia"      = "#b2df8a",
  "Bacteroidia"         = "#a6cee3",
  "Cytophagia"          = "#ffff99",
  "Saprospiria"         = "#cab2d6",
  "Chitinophagia"       = "#6a3d9a",
  "Clostridia"          = "#b15928",
  "Alphaproteobacteria" = "#E31A1C",
  "Gammaproteobacteria" = "#FB9A99",
  "Betaproteobacteria"  = "#fdbf6f",
  "Other"               = "grey70"
)

COMPLETENESS_COLORS <- c(
  "50-74"  = "#FDBF6F",  # orange
  "75-100" = "#33A02C"   # green
)

FAMILIES_TO_COLLAPSE <- c("Intestiviridae", "Crevaviridae", "Suoliviridae")

COLLAPSED_TRIANGLE_LENGTH <- 0.18
COLLAPSED_TRIANGLE_HALF_HEIGHT <- 8
COLLAPSED_TRIANGLE_ALPHA <- 0.95
COLLAPSED_TRIANGLE_BORDER <- 0.2

crass_families <- setdiff(names(CRASSVIRALES_COLOR_SCHEME), c("Outgroup", "Other"))

# ---------- helper functions ----------
get_parent_node <- function(tree, node) {
  p <- tree$edge[tree$edge[, 2] == node, 1]
  if (length(p) == 0) return(NA_integer_)
  as.integer(p[1])
}

get_descendants <- function(tree, node) {
  children <- tree$edge[tree$edge[, 1] == node, 2]
  if (length(children) == 0) return(integer(0))
  out <- children
  for (ch in children) {
    out <- c(out, get_descendants(tree, ch))
  }
  out
}

get_tip_descendants <- function(tree, node) {
  Ntip <- length(tree$tip.label)
  d <- get_descendants(tree, node)
  d[d <= Ntip]
}

leaf_to_accession <- function(x) {
  x %>%
    stringr::str_replace("\\|.*$", "") %>%
    stringr::str_replace("_[A-Z]{8}_CDS_[0-9]+$", "") %>%
    stringr::str_replace("\\r$", "") %>%
    stringr::str_trim()
}

prophage_fragment_to_accession <- function(x) {
  x %>%
    stringr::str_replace("_[0-9]+-[0-9]+$", "") %>%
    stringr::str_replace("\\r$", "") %>%
    stringr::str_trim()
}

find_family_pure_clades <- function(tree, annot_all, crass_families) {
  Ntip <- length(tree$tip.label)
  internal_nodes <- (Ntip + 1):(Ntip + tree$Nnode)

  node_info <- lapply(internal_nodes, function(node) {
    tip_ids <- get_tip_descendants(tree, node)
    tip_labels <- tree$tip.label[tip_ids]

    desc_df <- annot_all %>%
      dplyr::filter(label %in% tip_labels)

    crass_desc <- desc_df %>%
      dplyr::filter(family %in% crass_families)

    fams <- sort(unique(crass_desc$family))

    tibble(
      node = node,
      n_tips = length(tip_labels),
      n_crass = nrow(crass_desc),
      family = if (length(fams) == 1) fams else NA_character_,
      is_pure = length(fams) == 1 && nrow(crass_desc) > 0
    )
  }) %>%
    dplyr::bind_rows()

  pure_nodes <- node_info %>%
    dplyr::filter(is_pure)

  if (nrow(pure_nodes) == 0) {
    return(tibble(node = integer(), family = character(), n_tips = integer()))
  }

  # keep only maximal pure clades:
  # remove a node if its parent has the same pure family
  pure_nodes$parent <- vapply(pure_nodes$node, function(n) get_parent_node(tree, n), integer(1))

  parent_family_map <- pure_nodes %>%
    dplyr::select(node, family) %>%
    tibble::deframe()

  pure_nodes <- pure_nodes %>%
    dplyr::mutate(
      parent_family = unname(parent_family_map[as.character(parent)]),
      keep = is.na(parent_family) | parent_family != family
    ) %>%
    dplyr::filter(keep) %>%
    dplyr::select(node, family, n_tips, n_crass)

  pure_nodes
}

leaf_to_checkv_id <- function(x) {
  x %>%
    stringr::str_replace("_[A-Z]{8}_CDS_[0-9]+$", "") %>%
    stringr::str_replace("\\r$", "") %>%
    stringr::str_trim()
}

build_base_tree_plot <- function(tree_grouped, annot_all, open_angle = 180) {
  ggtree(
    tree_grouped,
    layout = "fan",
    open.angle = open_angle,
    aes(color = group, fill = group),
    size = 0.05
  ) %<+% annot_all +
    scale_color_manual(
      values = CRASSVIRALES_COLOR_SCHEME,
      na.value = CRASSVIRALES_COLOR_SCHEME["Other"]
    ) +
    scale_fill_manual(
      values = CRASSVIRALES_COLOR_SCHEME,
      na.value = CRASSVIRALES_COLOR_SCHEME["Other"]
    ) +
    theme_tree() +
    theme(
      legend.position = "none"
      # plot.margin = margin(20, 80, 20, 20)
      # plot.margin = margin(1, 3, 1, 1)
    )
}

build_collapsed_tree_from_pure_clades <- function(tree,
                                                  family_pure_clades,
                                                  families_to_collapse) {
  
  collapsed_node_family_df <- family_pure_clades %>%
    dplyr::filter(family %in% families_to_collapse) %>%
    dplyr::mutate(collapse_id = dplyr::row_number()) %>%
    dplyr::select(node, family, collapse_id)
  
  if (nrow(collapsed_node_family_df) == 0) {
    return(list(
      tree_collapsed = tree,
      marker_tip_df = tibble(
        marker_label = character(),
        family = character(),
        representative_tip = character(),
        removed_tip = character()
      )
    ))
  }
  
  tips_to_drop <- character(0)
  marker_tip_df_list <- list()
  
  for (i in seq_len(nrow(collapsed_node_family_df))) {
    node <- collapsed_node_family_df$node[i]
    fam  <- collapsed_node_family_df$family[i]
    cid  <- collapsed_node_family_df$collapse_id[i]
    
    clade_tips <- tree$tip.label[get_tip_descendants(tree, node)]
    
    if (length(clade_tips) == 0) {
      next
    }
    
    representative_tip <- clade_tips[1]
    removed_tips <- setdiff(clade_tips, representative_tip)
    marker_label <- paste0("collapsed_", fam, "_", cid)
    
    tips_to_drop <- c(tips_to_drop, removed_tips)
    
    marker_tip_df_list[[length(marker_tip_df_list) + 1]] <- tibble(
      marker_label = marker_label,
      family = fam,
      representative_tip = representative_tip
    )
  }
  
  tips_to_drop <- unique(tips_to_drop)
  
  tree_collapsed <- tree
  if (length(tips_to_drop) > 0) {
    tree_collapsed <- ape::drop.tip(tree_collapsed, tips_to_drop)
  }
  
  marker_tip_df <- dplyr::bind_rows(marker_tip_df_list)
  
  if (nrow(marker_tip_df) > 0) {
    idx <- match(marker_tip_df$representative_tip, tree_collapsed$tip.label)
    tree_collapsed$tip.label[idx] <- marker_tip_df$marker_label
  }
  
  list(
    tree_collapsed = tree_collapsed,
    marker_tip_df = marker_tip_df
  )
}

build_collapsed_annotation <- function(tree_collapsed,
                                       annot_all,
                                       marker_tip_df) {
  
  real_tip_labels <- setdiff(tree_collapsed$tip.label, marker_tip_df$marker_label)
  
  annot_real <- annot_all %>%
    dplyr::filter(label %in% real_tip_labels) %>%
    dplyr::mutate(
      is_collapsed_marker = FALSE,
      collapsed_family = NA_character_
    )
  
  annot_markers <- tibble(label = marker_tip_df$marker_label) %>%
    dplyr::left_join(marker_tip_df, by = c("label" = "marker_label")) %>%
    dplyr::mutate(
      genome_id = "collapsed_clade",
      family = family,
      is_crassvirales = TRUE,
      is_prophage = FALSE,
      integration_status = NA_character_,
      source_type = NA_character_,
      bact_phylum2 = NA_character_,
      bact_class2 = NA_character_,
      bact_genomad_length_num = NA_real_,
      bact_prophage_ratio_num = NA_real_,
      completeness = NA_real_,
      completeness_group = NA_character_,
      completeness_group_50plus = NA_character_,
      is_selected_candidate = FALSE,
      short_prophage_id = NA_integer_,
      short_prophage_label = NA_character_,
      is_collapsed_marker = TRUE,
      collapsed_family = family
    )
  
  dplyr::bind_rows(annot_real, annot_markers)
}

add_collapsed_tip_markers <- function(p,
                                      annot_all_collapsed,
                                      triangle_size = 2.5,
                                      triangle_stroke = 0.25) {
  
  marker_df <- p$data %>%
    dplyr::filter(isTip, label %in% annot_all_collapsed$label) %>%
    dplyr::select(label, x, y) %>%
    dplyr::left_join(
      annot_all_collapsed %>%
        dplyr::select(label, is_collapsed_marker, collapsed_family),
      by = "label"
    ) %>%
    dplyr::filter(is_collapsed_marker)
  
  if (nrow(marker_df) == 0) {
    return(p)
  }
  
  p +
    ggnewscale::new_scale_fill() +
    geom_point(
      data = marker_df,
      aes(x = x, y = y, fill = collapsed_family),
      inherit.aes = FALSE,
      shape = 24,
      size = triangle_size,
      stroke = triangle_stroke,
      colour = "black"
    ) +
    scale_fill_manual(
      values = CRASSVIRALES_COLOR_SCHEME,
      na.value = CRASSVIRALES_COLOR_SCHEME["Other"]
    )
}

build_and_save_short_prophage_ids <- function(p_for_order,
                                              annot_all,
                                              output_dir,
                                              stem,
                                              mode) {

  # ---- Build short prophage IDs from visual tree order ----
  tip_order_df <- p_for_order$data %>%
    dplyr::filter(
      isTip,
      label %in% annot_all$label
    ) %>%
    dplyr::select(label, x, y) %>%
    dplyr::left_join(
      annot_all %>%
        dplyr::select(
          label,
          is_prophage,
          integration_status,
          completeness,
          completeness_group,
          completeness_group_50plus
        ),
      by = "label"
    ) %>%
    dplyr::filter(is_prophage) %>%
    dplyr::arrange(desc(y)) %>%   # if numbering looks reversed, change to arrange(y)
    dplyr::mutate(
      short_prophage_id = dplyr::row_number()
    ) %>%
    dplyr::select(label, short_prophage_id)

  annot_all_mode <- annot_all %>%
    dplyr::left_join(tip_order_df, by = "label") %>%
    dplyr::mutate(
      short_prophage_label = paste0(
        short_prophage_id,
        dplyr::case_when(
          integration_status == "integrated"     ~ "i",
          integration_status == "non_integrated" ~ "n",
          TRUE                                   ~ ""
        ),
        dplyr::if_else(is_selected_candidate, "*", "")
      )
    )

  short_id_table <- annot_all_mode %>%
    dplyr::filter(is_prophage, !is.na(short_prophage_id)) %>%
    dplyr::arrange(short_prophage_id) %>%
    dplyr::select(
      short_prophage_id,
      short_prophage_label,
      is_selected_candidate,
      label,
      genome_id,
      family,
      bact_accession,
      bact_topology,
      integration_status,
      source_type,
      bact_genomad_length_num,
      bact_prophage_ratio_num,
      bact_phylum,
      bact_class,
      completeness,
      contamination,
      checkv_quality,
      miuvig_quality,
      completeness_group,
      completeness_group_50plus,
      everything()
    )

  out_short_id_tsv <- file.path(
    output_dir,
    paste0(stem, "_", mode, "_prophage_short_ids.tsv")
  )

  readr::write_tsv(short_id_table, out_short_id_tsv)

  cat("[OK] Saved short ID table:", out_short_id_tsv, "\n")

  return(list(
    tip_order_df = tip_order_df,
    annot_all_mode = annot_all_mode,
    short_id_table = short_id_table,
    out_short_id_tsv = out_short_id_tsv
  ))
}

annotate_genomad_length <- function(p, visible_labels = NULL,
                                    genomad_div = GENOMAD_DIV,
                                    genomad_levels = GENOMAD_LEVELS,
                                    tree_to_genomad_gap = TREE_TO_GENOMAD_GAP,
                                    text_offset = 0.03) {

  tree_data <- p$data %>%
    dplyr::filter(isTip, !is.na(y), !is.na(label)) %>%
    dplyr::select(
      label, x, y, is_prophage,
      bact_genomad_length_num,
      bact_topology,
      integration_status
    )

  if (!is.null(visible_labels)) {
    tree_data <- tree_data %>%
      dplyr::filter(label %in% visible_labels)
  }

  if (nrow(tree_data) == 0) {
    tree_outer_r <- max(p$data$x, na.rm = TRUE)
    genomad_base_r <- tree_outer_r + tree_to_genomad_gap

    genomad_ring_positions <- tibble(
      genomad_len = genomad_levels,
      ring_x = genomad_base_r + (genomad_levels / genomad_div)
    )

    return(list(
      tree_outer_r = tree_outer_r,
      genomad_base_r = genomad_base_r,
      genomad_ring_positions = genomad_ring_positions,
      genomad_ring_segments_df = tibble(),
      genomad_label_df = tibble(),
      genomad_ring_df = tibble(),
      guide_df = tibble()
    ))
  }

  tree_outer_r <- max(p$data$x, na.rm = TRUE)
  genomad_base_r <- tree_outer_r + tree_to_genomad_gap

  genomad_ring_positions <- tibble(
    genomad_len = genomad_levels,
    ring_x = genomad_base_r + (genomad_levels / genomad_div)
  )

  y_range <- range(tree_data$y, na.rm = TRUE)

  genomad_ring_segments_df <- genomad_ring_positions %>%
    dplyr::mutate(
      y = y_range[1],
      yend = y_range[2]
    )

  genomad_label_df <- genomad_ring_positions %>%
    dplyr::mutate(
      y = y_range[2],
      label = as.character(genomad_len / 1000),
      x = ring_x + text_offset
    )

  genomad_ring_df <- tree_data %>%
    dplyr::filter(
      is_prophage,
      !is.na(bact_genomad_length_num)
    ) %>%
    dplyr::mutate(
      geno_group = dplyr::if_else(
        !is.na(bact_topology) & bact_topology == "Provirus",
        "integrated_provirus",
        "non_integrated_prophage"
      ),
      geno_x = genomad_base_r,
      geno_xend = genomad_base_r + bact_genomad_length_num / genomad_div
    )

  guide_df <- genomad_ring_df %>%
    dplyr::filter(
      !is.na(x),
      !is.na(y),
      !is.na(geno_x),
      !is.na(integration_status),
      integration_status == "integrated"
    ) %>%
    dplyr::mutate(
      guide_x = x,
      guide_xend = geno_x
    )

  list(
    tree_outer_r = tree_outer_r,
    genomad_base_r = genomad_base_r,
    genomad_ring_positions = genomad_ring_positions,
    genomad_ring_segments_df = genomad_ring_segments_df,
    genomad_label_df = genomad_label_df,
    genomad_ring_df = genomad_ring_df,
    guide_df = guide_df
  )
}


draw_genomad_length <- function(p,
                                genomad_ann,
                                prophage_linewidth = PROPHAGE_LINEWIDTH) {

  p +
    geom_segment(
      data = genomad_ann$guide_df,
      aes(
        x = guide_x,
        xend = guide_xend,
        y = y,
        yend = y
      ),
      inherit.aes = FALSE,
      colour = "red",
      linetype = "dotted",
      linewidth = 0.1,
      lineend = "round"
    ) +
    geom_segment(
      data = genomad_ann$genomad_ring_segments_df,
      aes(
        x = ring_x,
        xend = ring_x,
        y = y,
        yend = yend
      ),
      inherit.aes = FALSE,
      colour = "grey50",
      linetype = "dotted",
      linewidth = 0.15,
      lineend = "round"
    ) +
    geom_text(
      data = genomad_ann$genomad_label_df,
      aes(
        x = x,
        y = y,
        label = label
      ),
      inherit.aes = FALSE,
      size = 1,
      vjust = -0.2,
      angle = 45,
      colour = "grey30"
    ) +
    ggnewscale::new_scale_color() +
    geom_segment(
      data = genomad_ann$genomad_ring_df,
      aes(
        x = geno_x,
        xend = geno_xend,
        y = y,
        yend = y,
        color = geno_group
      ),
      linewidth = prophage_linewidth,
      lineend = "butt",
      inherit.aes = FALSE
    ) +
    scale_color_manual(
      values = c(
        "integrated_provirus" = "red",
        "non_integrated_prophage" = "blue"
      )
    )
}

annotate_ratio <- function(p,
                           genomad_ann,
                           visible_labels = NULL,
                           ratio_levels = RATIO_LEVELS,
                           ratio_mult = RATIO_MULT,
                           genomad_to_ratio_gap = GENOMAD_TO_RATIO_GAP,
                           text_extra_gap = TEXT_EXTRA_GAP) {

  tree_data <- p$data %>%
    dplyr::filter(isTip, !is.na(y), !is.na(label)) %>%
    dplyr::select(
      label, x, y, is_prophage,
      bact_prophage_ratio_num,
      bact_topology
    )

  if (!is.null(visible_labels)) {
    tree_data <- tree_data %>%
      dplyr::filter(label %in% visible_labels)
  }

  if (nrow(tree_data) == 0) {
    genomad_outer_r <- if (nrow(genomad_ann$genomad_ring_df) > 0) {
      max(genomad_ann$genomad_ring_df$geno_xend, na.rm = TRUE)
    } else {
      genomad_ann$genomad_base_r
    }

    ratio_base_r <- genomad_outer_r + genomad_to_ratio_gap

    ratio_ring_positions <- tibble(
      ratio_value = ratio_levels,
      ring_x = ratio_base_r + ratio_levels * ratio_mult
    )

    return(list(
      genomad_outer_r = genomad_outer_r,
      ratio_base_r = ratio_base_r,
      ratio_ring_positions = ratio_ring_positions,
      ratio_ring_segments_df = tibble(),
      ratio_label_df = tibble(),
      ratio_df = tibble(),
      ratio_outer_r = ratio_base_r
    ))
  }

  y_range <- range(tree_data$y, na.rm = TRUE)

  genomad_outer_r <- if (nrow(genomad_ann$genomad_ring_df) > 0) {
    max(genomad_ann$genomad_ring_df$geno_xend, na.rm = TRUE)
  } else {
    genomad_ann$genomad_base_r
  }

  ratio_base_r <- genomad_outer_r + genomad_to_ratio_gap

  ratio_ring_positions <- tibble(
    ratio_value = ratio_levels,
    ring_x = ratio_base_r + ratio_levels * ratio_mult
  )

  ratio_ring_segments_df <- ratio_ring_positions %>%
    dplyr::mutate(
      y = y_range[1],
      yend = y_range[2]
    )

  ratio_label_df <- ratio_ring_positions %>%
    dplyr::mutate(
      y = y_range[2],
      label = sprintf("%.2f", ratio_value),
      x = ring_x + text_extra_gap
    )

  ratio_df <- tree_data %>%
    dplyr::filter(
      is_prophage,
      !is.na(bact_prophage_ratio_num)
    ) %>%
    dplyr::mutate(
      ratio_group = dplyr::if_else(
        !is.na(bact_topology) & bact_topology == "Provirus",
        "integrated_provirus",
        "non_integrated_prophage"
      ),
      ratio_x = ratio_base_r,
      ratio_xend = ratio_base_r + bact_prophage_ratio_num * ratio_mult
    ) %>%
    dplyr::filter(!is.na(ratio_xend))

  ratio_outer_r <- if (nrow(ratio_df) > 0) {
    max(ratio_df$ratio_xend, na.rm = TRUE)
  } else {
    ratio_base_r
  }

  list(
    genomad_outer_r = genomad_outer_r,
    ratio_base_r = ratio_base_r,
    ratio_ring_positions = ratio_ring_positions,
    ratio_ring_segments_df = ratio_ring_segments_df,
    ratio_label_df = ratio_label_df,
    ratio_df = ratio_df,
    ratio_outer_r = ratio_outer_r
  )
}

draw_ratio <- function(p,
                       ratio_ann,
                       prophage_ratio_linewidth = PROPHAGE_RATIO_LINEWIDTH) {

  p +
    geom_segment(
      data = ratio_ann$ratio_ring_segments_df,
      aes(
        x = ring_x,
        xend = ring_x,
        y = y,
        yend = yend
      ),
      inherit.aes = FALSE,
      colour = "grey50",
      linetype = "dotted",
      linewidth = 0.15,
      lineend = "round"
    ) +
    geom_text(
      data = ratio_ann$ratio_label_df,
      aes(
        x = x,
        y = y,
        label = label
      ),
      inherit.aes = FALSE,
      size = 1,
      vjust = -0.2,
      angle = 45,
      colour = "grey30"
    ) +
    ggnewscale::new_scale_color() +
    geom_segment(
      data = ratio_ann$ratio_df,
      aes(
        x = ratio_x,
        xend = ratio_xend,
        y = y,
        yend = y,
        color = ratio_group
      ),
      linewidth = prophage_ratio_linewidth,
      lineend = "butt",
      inherit.aes = FALSE
    ) +
    scale_color_manual(
      values = c(
        "integrated_provirus"     = "red",
        "non_integrated_prophage" = "blue"
      )
    )
}

annotate_phylum_class <- function(p,
                                  ratio_ann,
                                  visible_labels = NULL,
                                  ratio_to_phylum_gap = RATIO_TO_PHYLUM_GAP,
                                  phylum_ring_width = PHYLUM_RING_WIDTH,
                                  phylum_to_class_gap = PHYLUM_TO_CLASS_GAP,
                                  class_ring_width = CLASS_RING_WIDTH) {

  tree_data <- p$data %>%
    dplyr::filter(isTip, !is.na(y), !is.na(label)) %>%
    dplyr::select(
      label, y, is_prophage,
      bact_phylum2,
      bact_class2
    )

  if (!is.null(visible_labels)) {
    tree_data <- tree_data %>%
      dplyr::filter(label %in% visible_labels)
  }

  ratio_outer_r <- ratio_ann$ratio_outer_r

  phylum_base_r <- ratio_outer_r + ratio_to_phylum_gap

  phylum_ring_df <- tree_data %>%
    dplyr::filter(
      is_prophage,
      !is.na(bact_phylum2)
    ) %>%
    dplyr::mutate(
      phyl_x = phylum_base_r,
      phyl_xend = phylum_base_r + phylum_ring_width
    )

  class_base_r <- phylum_base_r + phylum_ring_width + phylum_to_class_gap

  class_ring_df <- tree_data %>%
    dplyr::filter(
      is_prophage,
      !is.na(bact_class2)
    ) %>%
    dplyr::mutate(
      class_x = class_base_r,
      class_xend = class_base_r + class_ring_width
    )

  class_outer_r <- if (nrow(class_ring_df) > 0) {
    max(class_ring_df$class_xend, na.rm = TRUE)
  } else {
    class_base_r
  }

  list(
    ratio_outer_r = ratio_outer_r,
    phylum_base_r = phylum_base_r,
    phylum_ring_df = phylum_ring_df,
    class_base_r = class_base_r,
    class_ring_df = class_ring_df,
    class_outer_r = class_outer_r
  )
}

draw_phylum_class <- function(p,
                              phylum_class_ann,
                              bar_linewidth = BAR_LINEWIDTH,
                              phylum_colors = PHYLUM_COLORS,
                              class_colors = CLASS_COLORS) {

  p +
    ggnewscale::new_scale_color() +
    geom_segment(
      data = phylum_class_ann$phylum_ring_df,
      aes(
        x = phyl_x,
        xend = phyl_xend,
        y = y,
        yend = y,
        colour = bact_phylum2
      ),
      linewidth = bar_linewidth,
      lineend = "butt",
      inherit.aes = FALSE,
      show.legend = FALSE
    ) +
    scale_color_manual(
      values = phylum_colors,
      na.value = phylum_colors["Other"]
    ) +
    ggnewscale::new_scale_color() +
    geom_segment(
      data = phylum_class_ann$class_ring_df,
      aes(
        x = class_x,
        xend = class_xend,
        y = y,
        yend = y,
        colour = bact_class2
      ),
      linewidth = bar_linewidth,
      lineend = "butt",
      inherit.aes = FALSE,
      show.legend = FALSE
    ) +
    scale_color_manual(
      values = class_colors,
      na.value = class_colors["Other"]
    )
}

annotate_source_ring <- function(p,
                                 phylum_class_ann,
                                 visible_labels = NULL,
                                 source_ring_gap = SOURCE_RING_GAP,
                                 source_ring_width = SOURCE_RING_WIDTH) {

  tree_data <- p$data %>%
    dplyr::filter(isTip, !is.na(y), !is.na(label)) %>%
    dplyr::select(
      label, y, is_prophage, source_type
    )

  if (!is.null(visible_labels)) {
    tree_data <- tree_data %>%
      dplyr::filter(label %in% visible_labels)
  }

  source_base_r <- phylum_class_ann$class_outer_r + source_ring_gap

  source_ring_df <- tree_data %>%
    dplyr::filter(
      is_prophage,
      !is.na(source_type)
    ) %>%
    dplyr::mutate(
      source_x = source_base_r,
      source_xend = source_base_r + source_ring_width
    )

  source_outer_r <- if (nrow(source_ring_df) > 0) {
    max(source_ring_df$source_xend, na.rm = TRUE)
  } else {
    source_base_r
  }

  list(
    source_base_r = source_base_r,
    source_ring_df = source_ring_df,
    source_outer_r = source_outer_r
  )
}

draw_source_ring <- function(p,
                             source_ann,
                             bar_linewidth = BAR_LINEWIDTH,
                             source_colors = SOURCE_COLORS) {

  p +
    ggnewscale::new_scale_color() +
    geom_segment(
      data = source_ann$source_ring_df,
      aes(
        x = source_x,
        xend = source_xend,
        y = y,
        yend = y,
        colour = source_type
      ),
      linewidth = bar_linewidth,
      lineend = "butt",
      inherit.aes = FALSE,
      show.legend = FALSE
    ) +
    scale_color_manual(
      values = source_colors
    )
}

annotate_completeness_ring <- function(p,
                                       source_ann,
                                       visible_labels = NULL,
                                       source_to_completeness_gap = SOURCE_TO_COMPLETENESS_GAP,
                                       completeness_ring_width = COMPLETENESS_RING_WIDTH) {

  tree_data <- p$data %>%
    dplyr::filter(isTip, !is.na(y), !is.na(label)) %>%
    dplyr::select(
      label, y, is_prophage, completeness_group_50plus
    )

  if (!is.null(visible_labels)) {
    tree_data <- tree_data %>%
      dplyr::filter(label %in% visible_labels)
  }

  if (nrow(tree_data) == 0) {
    completeness_base_r <- source_ann$source_outer_r + source_to_completeness_gap

    plot_outer_limit <- max(
      c(p$data$x, completeness_base_r),
      na.rm = TRUE
    )

    return(list(
      completeness_base_r = completeness_base_r,
      completeness_ring_df = tibble(),
      completeness_outer_r = completeness_base_r,
      plot_outer_limit = plot_outer_limit
    ))
  }

  completeness_base_r <- source_ann$source_outer_r + source_to_completeness_gap

  completeness_ring_df <- tree_data %>%
    dplyr::filter(
      is_prophage,
      !is.na(completeness_group_50plus)
    ) %>%
    dplyr::mutate(
      comp_x = completeness_base_r,
      comp_xend = completeness_base_r + completeness_ring_width
    )

  completeness_outer_r <- if (nrow(completeness_ring_df) > 0) {
    max(completeness_ring_df$comp_xend, na.rm = TRUE)
  } else {
    completeness_base_r
  }

  plot_outer_limit <- max(
    c(p$data$x, completeness_outer_r),
    na.rm = TRUE
  )

  list(
    completeness_base_r = completeness_base_r,
    completeness_ring_df = completeness_ring_df,
    completeness_outer_r = completeness_outer_r,
    plot_outer_limit = plot_outer_limit
  )
}

draw_completeness_ring <- function(p,
                                   completeness_ann,
                                   bar_linewidth = BAR_LINEWIDTH,
                                   completeness_colors = COMPLETENESS_COLORS,
                                   x_padding = 0.05) {

  p +
    ggnewscale::new_scale_color() +
    geom_segment(
      data = completeness_ann$completeness_ring_df,
      aes(
        x = comp_x,
        xend = comp_xend,
        y = y,
        yend = y,
        colour = completeness_group_50plus
      ),
      linewidth = bar_linewidth,
      lineend = "butt",
      inherit.aes = FALSE,
      show.legend = FALSE
    ) +
    scale_color_manual(
      values = completeness_colors
    ) +
    xlim(NA, completeness_ann$plot_outer_limit + x_padding)
}


annotate_short_id_labels <- function(p,
                                     annot_all_mode,
                                     outer_limit,
                                     visible_labels = NULL,
                                     short_id_to_label_gap = SHORT_ID_TO_LABEL_GAP) {

  tree_data <- p$data %>%
    dplyr::filter(isTip, !is.na(y), !is.na(label)) %>%
    dplyr::select(label, x, y, angle)

  if (!is.null(visible_labels)) {
    tree_data <- tree_data %>%
      dplyr::filter(label %in% visible_labels)
  }

  if (nrow(tree_data) == 0) {
    return(list(
      label_df = tibble(),
      label_outer_limit = outer_limit
    ))
  }

  label_df <- tree_data %>%
    dplyr::left_join(
      annot_all_mode %>%
        dplyr::select(
          label,
          short_prophage_id,
          short_prophage_label,
          is_prophage,
          integration_status,
          completeness
        ),
      by = "label"
    ) %>%
    dplyr::filter(
      is_prophage,
      !is.na(short_prophage_id),
      integration_status == "integrated" |
        (integration_status == "non_integrated" & !is.na(completeness) & completeness >= 50)
    ) %>%
    dplyr::mutate(
      label_x = outer_limit + short_id_to_label_gap,
      text_angle = ifelse(angle > 90 & angle < 270, angle + 180, angle),
      text_hjust = ifelse(angle > 90 & angle < 270, 1, 0)
    )

  label_outer_limit <- if (nrow(label_df) > 0) {
    max(label_df$label_x, na.rm = TRUE)
  } else {
    outer_limit
  }

  list(
    label_df = label_df,
    label_outer_limit = label_outer_limit
  )
}

draw_short_id_labels <- function(p,
                                 label_ann,
                                 short_id_text_size = SHORT_ID_TEXT_SIZE,
                                 short_id_text_color = SHORT_ID_TEXT_COLOR,
                                 x_padding = 0.05) {

  p +
    geom_text(
      data = label_ann$label_df,
      aes(
        x = label_x,
        y = y,
        label = short_prophage_label,
        angle = text_angle,
        hjust = text_hjust
      ),
      inherit.aes = FALSE,
      size = short_id_text_size,
      colour = short_id_text_color
    ) +
    xlim(NA, label_ann$label_outer_limit + x_padding)
}

save_plot_pair <- function(p_full,
                           p_collapsed,
                           output_dir,
                           stem,
                           mode,
                           width = 21.0,
                           height = 29.7,
                           units = "cm",
                           dpi = 1200) {

  out_png_full <- file.path(
    output_dir,
    paste0(stem, "_", mode, "_bare_fan180.png")
  )

  out_png_collapsed <- file.path(
    output_dir,
    paste0(stem, "_", mode, "_bare_fan180_collapsed.png")
  )

  ggsave(
    out_png_full,
    p_full,
    width = width,
    height = height,
    units = units,
    dpi = dpi
  )

  ggsave(
    out_png_collapsed,
    p_collapsed,
    width = width,
    height = height,
    units = units,
    dpi = dpi
  )

  cat("Saved:", out_png_full, "\n")
  cat("Saved:", out_png_collapsed, "\n")

  list(
    out_png_full = out_png_full,
    out_png_collapsed = out_png_collapsed
  )
}

root_tree_by_mode <- function(tree, mode,
                              annot_all,
                              crass_families,
                              outgroup_tips,
                              stem) {

  if (mode == "unrooted") {
    message("[ROOT] Unrooted tree")
    return(tree)
  }

  if (mode == "outgroup") {
    if (length(outgroup_tips) != 1 || !(outgroup_tips %in% tree$tip.label)) {
      stop("[ROOT] Outgroup rooting requested but invalid outgroup tip for: ", stem)
    }

    message("[ROOT] Rooting at outgroup: ", outgroup_tips)
    tree <- ape::root(tree, outgroup = outgroup_tips, resolve.root = TRUE)
    tree <- ape::reorder.phylo(tree, "cladewise")
    return(tree)
  }

  if (mode == "crass_mrca") {
    crass_tips <- annot_all %>%
      filter(family %in% crass_families) %>%
      pull(label) %>%
      intersect(tree$tip.label)

    if (length(crass_tips) < 2) {
      stop("[ROOT] <2 Crass tips for MRCA rooting in: ", stem)
    }

    noncrass_tips <- setdiff(tree$tip.label, crass_tips)
    if (length(noncrass_tips) < 1) {
      stop("[ROOT] No non-Crass tips available for MRCA rooting in: ", stem)
    }

    if (length(outgroup_tips) == 1 && outgroup_tips %in% noncrass_tips) {
      out_tip <- outgroup_tips
    } else {
      out_tip <- noncrass_tips[1]
    }

    message("[ROOT] Rooting at Crass MRCA using outgroup tip: ", out_tip)
    tree <- ape::root(tree, outgroup = out_tip, resolve.root = TRUE)
    tree <- ape::reorder.phylo(tree, "cladewise")
    return(tree)
  }

  if (mode == "crass_plus_integrated_mrca") {
    crass_tips <- annot_all %>%
      dplyr::filter(family %in% crass_families) %>%
      dplyr::pull(label) %>%
      intersect(tree$tip.label)

    integ_tips <- annot_all %>%
      dplyr::filter(!is.na(integrated_prophage_label), integrated_prophage_label != "") %>%
      dplyr::pull(label) %>%
      intersect(tree$tip.label)

    target_tips <- unique(c(crass_tips, integ_tips))

    if (length(target_tips) < 2) {
      stop("[ROOT] <2 tips for crass+integrated MRCA rooting in: ", stem)
    }

    non_target_tips <- setdiff(tree$tip.label, target_tips)
    if (length(non_target_tips) < 1) {
      stop("[ROOT] No tips outside (crass+integrated) set for rooting in: ", stem)
    }

    if (length(outgroup_tips) == 1 && outgroup_tips %in% non_target_tips) {
      out_tip <- outgroup_tips
    } else {
      out_tip <- non_target_tips[1]
    }

    message("[ROOT] Rooting at MRCA(crass refs + integrated-labeled) using outgroup tip: ", out_tip)
    tree <- ape::root(tree, outgroup = out_tip, resolve.root = TRUE)
    tree <- ape::reorder.phylo(tree, "cladewise")
    return(tree)
  }

  stop("Unknown rooting mode: ", mode)
}

source_df <- read_tsv(SOURCE_FILE, show_col_types = FALSE) %>%
  mutate(contig = as.character(contig)) %>%
  distinct(contig, .keep_all = TRUE)

cat("[INFO] source_df rows:", nrow(source_df), "\n")

# ---- Build mapping: accession -> prophage_label ----
#PAIRS_TSV <- file.path(WD, "pairs_prophage_with_flanks_vs_bacterial.tsv")

pairs <- readr::read_tsv(PAIRS_TSV, show_col_types = FALSE) %>%
  transmute(
    prophage = as.character(prophage),
    prophage_label = as.character(prophage_label),
    bact_accession = prophage_fragment_to_accession(prophage)
  ) %>%
  filter(!is.na(bact_accession), bact_accession != "",
         !is.na(prophage_label), prophage_label != "")

pairs_map <- pairs %>%
  distinct(bact_accession, .keep_all = TRUE) %>%
  transmute(
    bact_accession,
    integrated_prophage_label = prophage_label
  )

cat("[INFO] pairs_map rows:", nrow(pairs_map), "\n")

checkv_df <- read_tsv(CHECKV_FILE, show_col_types = FALSE) %>%
  mutate(
    checkv_id = as.character(contig_id),
    completeness = suppressWarnings(as.numeric(completeness)),
    contamination = suppressWarnings(as.numeric(contamination)),
    completeness_group = case_when(
      is.na(completeness) ~ NA_character_,
      completeness < 25 ~ "0-24",
      completeness < 50 ~ "25-49",
      completeness < 75 ~ "50-74",
      TRUE ~ "75-100"
    ),
    completeness_group_50plus = case_when(
      is.na(completeness) ~ NA_character_,
      completeness >= 75 ~ "75-100",
      completeness >= 50 ~ "50-74",
      TRUE ~ NA_character_
    )
  ) %>%
  distinct(checkv_id, .keep_all = TRUE)

cat("[INFO] checkv_df rows:", nrow(checkv_df), "\n")

candidate_ids <- readr::read_lines(CANDIDATE_LIST_FILE) %>%
  stringr::str_replace("\\r$", "") %>%
  stringr::str_trim() %>%
  .[. != ""] %>%
  unique()

cat("[INFO] candidate_ids loaded:", length(candidate_ids), "\n")

for (tree_file in TREE_FILES) {

  message("==== Processing tree file: ", tree_file)

  stem <- tools::file_path_sans_ext(basename(tree_file))

  ANNOT_FILE      <- file.path(DATA_DIR, paste0(stem, "_tree_protein_taxonomy.tsv"))
  BACT_ANNOT_FILE <- file.path(DATA_DIR, paste0(stem, "_tree_bacterial_taxonomy_with_genomad_and_lengths.tsv"))

  tree_raw <- read.tree(tree_file)
  annot_raw <- read_tsv(ANNOT_FILE, show_col_types = FALSE)
  bact_annot_raw <- read_tsv(BACT_ANNOT_FILE, show_col_types = FALSE)

  annot_all <- tibble(label = tree_raw$tip.label) %>%
    left_join(annot_raw %>% rename(label = leaf_label), by = "label") %>%
    mutate(
      genome_id = if_else(is.na(genome_id), "unknown", genome_id),
      family    = if_else(is.na(family), "unknown", family)
    )

  outgroup_genomes <- c("NC_021803")
  outgroup_tips <- annot_all %>%
    filter(genome_id %in% outgroup_genomes) %>%
    pull(label)

  bact_annot <- bact_annot_raw %>%
    rename(label = leaf_label) %>%
    rename_with(~ ifelse(startsWith(.x, "bact_"), .x, paste0("bact_", .x)), -label)

  print(colnames(bact_annot))

  annot_all <- annot_all %>%
    left_join(bact_annot, by = "label") %>%
    mutate(
      is_crassvirales = genome_id != "unknown" & !(genome_id %in% outgroup_genomes),
      is_prophage     = !is.na(bact_domain) & bact_domain == "Bacteria",

      bact_genomad_length_num = suppressWarnings(as.numeric(bact_genomad_length)),
      bact_contig_length_num  = suppressWarnings(as.numeric(bact_bacterial_length)),
      bact_prophage_ratio_num = suppressWarnings(as.numeric(bact_prophage_ratio)),

      prophage_ratio = bact_genomad_length_num / bact_contig_length_num,
      prophage_ratio = pmin(prophage_ratio, 1)
    ) %>%
    mutate(
      contig_id = label %>%
        stringr::str_replace("\\|.*$", "") %>%
        stringr::str_replace("_[A-Z]{8}_CDS_[0-9]+$", "")
    ) %>%
    left_join(source_df, by = c("contig_id" = "contig")) %>%
    mutate(bact_accession = leaf_to_accession(label)) %>%
    left_join(pairs_map, by = "bact_accession") %>%
    mutate(checkv_id = leaf_to_checkv_id(label)) %>%
    left_join(
      checkv_df %>%
        select(
          checkv_id,
          completeness,
          contamination,
          checkv_quality,
          miuvig_quality,
          completeness_group,
          completeness_group_50plus
        ),
      by = "checkv_id"
    ) %>%
    mutate(
      integration_status = case_when(
        is_prophage & !is.na(bact_topology) & bact_topology == "Provirus" ~ "integrated",
        is_prophage ~ "non_integrated",
        TRUE ~ NA_character_
      ),
      source_type = case_when(
        sample_origin == "isolate" ~ "isolate",
        sample_origin == "MAG"     ~ "MAG",
        TRUE                       ~ NA_character_
      )
    )

  annot_all <- annot_all %>%
    mutate(
      bact_phylum2 = case_when(
        is.na(bact_phylum) ~ NA_character_,
        bact_phylum %in% names(PHYLUM_COLORS) ~ bact_phylum,
        TRUE ~ "Other"
      ),
      bact_class2 = case_when(
        is.na(bact_class) ~ NA_character_,
        bact_class %in% names(CLASS_COLORS) ~ bact_class,
        TRUE ~ "Other"
      ),
      bact_seq_name = bact_seq_name %>%
        as.character() %>%
        stringr::str_replace("\\r$", "") %>%
        stringr::str_trim(),
      is_selected_candidate = !is.na(bact_seq_name) &
        bact_seq_name %in% candidate_ids
    )

  cat("[INFO] Tips with bact_seq_name in candidate list:",
      sum(annot_all$is_selected_candidate, na.rm = TRUE), "\n")

  cat("[INFO] Tips with integrated_prophage_label:",
      sum(!is.na(annot_all$integrated_prophage_label) & annot_all$integrated_prophage_label != ""),
      "\n")

  for (mode in ROOTING_MODES) {

    message("---- Rooting mode: ", mode)

    tree <- root_tree_by_mode(
      tree = tree_raw,
      mode = mode,
      annot_all = annot_all,
      crass_families = crass_families,
      outgroup_tips = outgroup_tips,
      stem = stem
    )

    family_pure_clades <- find_family_pure_clades(
      tree = tree,
      annot_all = annot_all,
      crass_families = crass_families
    )

    print(family_pure_clades)

    # ---- MRCA of all Crassvirales tips ----
    crass_tips <- annot_all %>%
      dplyr::filter(family %in% crass_families) %>%
      dplyr::pull(label) %>%
      intersect(tree$tip.label)

    if (length(crass_tips) < 2) {
      stop("Cannot compute Crassvirales MRCA: <2 Crassvirales tips in tree for: ", stem)
    }

    crass_mrca_node <- ape::getMRCA(tree, crass_tips)

    # ---- Get all descendant tips of Crass MRCA ----
    crass_mrca_tip_ids <- get_tip_descendants(tree, crass_mrca_node)
    crass_mrca_tip_labels <- tree$tip.label[crass_mrca_tip_ids]

    # ---- Filter prophage tips ----
    prophage_labels <- annot_all %>%
      dplyr::filter(
        label %in% crass_mrca_tip_labels,
        is_prophage
      ) %>%
      dplyr::pull(label)

    # ---- Output file ----
    out_txt <- file.path(
      OUTPUT_DIR,
      paste0(stem, "_", mode, "_prophage_in_crass_mrca.txt")
    )

    # ---- Save (one label per line) ----
    writeLines(prophage_labels, out_txt)

    cat("[OK] Saved prophage labels in Crass MRCA:",
        out_txt,
        " (n=", length(prophage_labels), ")\n", sep = "")

    cat("[DEBUG] Total tips in MRCA:", length(crass_mrca_tip_labels), "\n")
    cat("[DEBUG] Prophage tips in MRCA:", length(prophage_labels), "\n")

    if (is.null(crass_mrca_node) || is.na(crass_mrca_node)) {
      stop("Could not compute Crassvirales MRCA for: ", stem, " in mode: ", mode)
    }

    cat("[INFO] Crassvirales MRCA node =", crass_mrca_node, "\n")

    # attach group annotations to tree
    groups_list <- lapply(crass_families, function(fam) {
      annot_all %>%
        filter(family == fam) %>%
        pull(label)
    })
    names(groups_list) <- crass_families
    groups_list[["Outgroup"]] <- outgroup_tips

    tree_grouped_full <- groupOTU(tree, groups_list)
    
    collapsed_res <- build_collapsed_tree_from_pure_clades(
      tree = tree,
      family_pure_clades = family_pure_clades,
      families_to_collapse = FAMILIES_TO_COLLAPSE
    )
    
    tree_collapsed <- collapsed_res$tree_collapsed
    marker_tip_df <- collapsed_res$marker_tip_df
    
    groups_list_collapsed <- lapply(crass_families, function(fam) {
      tips <- annot_all %>%
        dplyr::filter(family == fam) %>%
        dplyr::pull(label)
      
      marker_tips <- marker_tip_df %>%
        dplyr::filter(family == fam) %>%
        dplyr::pull(marker_label)
      
      intersect(c(tips, marker_tips), tree_collapsed$tip.label)
    })
    names(groups_list_collapsed) <- crass_families
    
    outgroup_tips_collapsed <- intersect(outgroup_tips, tree_collapsed$tip.label)
    groups_list_collapsed[["Outgroup"]] <- outgroup_tips_collapsed
    
    tree_grouped_collapsed <- groupOTU(tree_collapsed, groups_list_collapsed)
    
    annot_all_collapsed <- build_collapsed_annotation(
      tree_collapsed = tree_collapsed,
      annot_all = annot_all,
      marker_tip_df = marker_tip_df
    )
    
    family_pure_clades_collapsed <- find_family_pure_clades(
      tree = tree_collapsed,
      annot_all = annot_all_collapsed,
      crass_families = crass_families
    )
    
    OPEN_ANGLE_FULL <- 180
    OPEN_ANGLE_COLLAPSED <- 180
    
    p_base_full <- build_base_tree_plot(
      tree_grouped = tree_grouped_full,
      annot_all = annot_all,
      open_angle = OPEN_ANGLE_FULL
    )
    
    p_base_collapsed <- build_base_tree_plot(
      tree_grouped = tree_grouped_collapsed,
      annot_all = annot_all_collapsed,
      open_angle = OPEN_ANGLE_COLLAPSED
    )
    
    p_full <- p_base_full
    p_collapsed <- p_base_collapsed
    
    p_collapsed <- add_collapsed_tip_markers(
      p = p_collapsed,
      annot_all_collapsed = annot_all_collapsed
    )

    if (nrow(family_pure_clades) > 0) {
      for (i in seq_len(nrow(family_pure_clades))) {
        fam <- family_pure_clades$family[i]
        node <- family_pure_clades$node[i]
        
        p_full <- p_full +
          geom_hilight(
            node = node,
            fill = CRASS_CLADE_FILL[[fam]],
            alpha = 0.25,
            extend = 0.002
          )
      }
    }
    
    if (nrow(family_pure_clades_collapsed) > 0) {
      for (i in seq_len(nrow(family_pure_clades_collapsed))) {
        fam <- family_pure_clades_collapsed$family[i]
        node <- family_pure_clades_collapsed$node[i]
        
        if (!(fam %in% FAMILIES_TO_COLLAPSE)) {
          p_collapsed <- p_collapsed +
            geom_hilight(
              node = node,
              fill = CRASS_CLADE_FILL[[fam]],
              alpha = 0.25,
              extend = 0.002
            )
        }
      }
    }
    
    short_id_res <- build_and_save_short_prophage_ids(
      p_for_order = p_base_full,
      annot_all = annot_all,
      output_dir = OUTPUT_DIR,
      stem = stem,
      mode = mode
    )

    annot_all_mode <- short_id_res$annot_all_mode
    
    annot_all_collapsed_mode <- annot_all_collapsed %>%
      dplyr::select(-short_prophage_id, -short_prophage_label) %>%
      dplyr::left_join(
        annot_all_mode %>%
          dplyr::select(
            label,
            short_prophage_id,
            short_prophage_label
          ),
        by = "label"
      )

    plot_variants <- list(
      full = list(
        p = p_full,
        visible_labels = NULL,
        annot = annot_all_mode,
        short_id_text_size = SHORT_ID_TEXT_SIZE
      ),
      collapsed = list(
        p = p_collapsed,
        visible_labels = NULL,
        annot = annot_all_collapsed_mode,
        short_id_text_size = SHORT_ID_TEXT_SIZE_COLLAPSED
      )
    )

    for (variant_name in names(plot_variants)) {

      p_current <- plot_variants[[variant_name]]$p
      visible_labels <- plot_variants[[variant_name]]$visible_labels
      short_id_text_size_current <- plot_variants[[variant_name]]$short_id_text_size

      genomad_ann <- annotate_genomad_length(
        p = p_current,
        visible_labels = visible_labels
      )
      p_current <- draw_genomad_length(p_current, genomad_ann)

      ratio_ann <- annotate_ratio(
        p = p_current,
        genomad_ann = genomad_ann,
        visible_labels = visible_labels
      )
      p_current <- draw_ratio(p_current, ratio_ann)

      phylum_class_ann <- annotate_phylum_class(
        p = p_current,
        ratio_ann = ratio_ann,
        visible_labels = visible_labels
      )
      p_current <- draw_phylum_class(p_current, phylum_class_ann)

      source_ann <- annotate_source_ring(
        p = p_current,
        phylum_class_ann = phylum_class_ann,
        visible_labels = visible_labels
      )
      p_current <- draw_source_ring(p_current, source_ann)

      completeness_ann <- annotate_completeness_ring(
        p = p_current,
        source_ann = source_ann,
        visible_labels = visible_labels
      )
      p_current <- draw_completeness_ring(p_current, completeness_ann)

      label_ann <- annotate_short_id_labels(
        p = p_current,
        annot_all_mode = plot_variants[[variant_name]]$annot,
        outer_limit = completeness_ann$plot_outer_limit,
        visible_labels = visible_labels
      )
      p_current <- draw_short_id_labels(
        p = p_current,
        label_ann = label_ann,
        short_id_text_size = short_id_text_size_current
      )

      plot_variants[[variant_name]]$p <- p_current
    }

    p_full <- plot_variants$full$p
    p_collapsed <- plot_variants$collapsed$p

    save_plot_pair(
      p_full = p_full,
      p_collapsed = p_collapsed,
      output_dir = OUTPUT_DIR,
      stem = stem,
      mode = mode
    )
  }
}