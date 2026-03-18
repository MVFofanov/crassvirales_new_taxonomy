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
WD <- "C:/crassvirales/crassvirales_new_taxonomy/crassvirales_prophages/tree_plots_for_manuscript"
DATA_DIR <- file.path(WD, "data")

# Output directory
OUTPUT_DIR <- file.path(WD, "figures")

# Shared input files
SOURCE_FILE <- file.path(DATA_DIR, "all_genomad_contigs_mags_or_isolates.txt")
PAIRS_TSV   <- file.path(DATA_DIR, "pairs_prophage_with_flanks_vs_bacterial.tsv")

# Tree files
TREE_FILES <- c(
  file.path(DATA_DIR, "TerL_gappy0.8.treefile")
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
TREE_TO_GENOMAD_GAP <- 0.03
GENOMAD_TO_RATIO_GAP <- 0.03
RATIO_TO_PHYLUM_GAP <- 0.03
PHYLUM_TO_CLASS_GAP <- 0.02
CLASS_TO_DOT_GAP <- 0.03
TEXT_EXTRA_GAP <- 0.03

GENOMAD_DIV <- 500000
GENOMAD_LEVELS <- c(50000, 100000, 150000)

RATIO_LEVELS <- c(0.25, 0.50, 0.75, 1.00)
RATIO_MULT <- max(GENOMAD_LEVELS) / GENOMAD_DIV

# 50% of max prophage-size bar width (= 0.3 / 2 = 0.15)
PHYLUM_RING_WIDTH <- 0.15
CLASS_RING_WIDTH  <- 0.15

DOT_SIZE    <- 0.4
DOT_STROKE  <- 0.15

ROOTING_MODES <- c(
  "outgroup",
  "crass_mrca",
  "crass_plus_integrated_mrca"
)

CRASSVIRALES_COLOR_SCHEME <- c(
  "Intestiviridae" = "#EE3B3B",
  "Crevaviridae"   = "#EE9A00",
  "Suoliviridae"   = "#4169E1",
  "Steigviridae"   = "#00CED1",
  "Epsilon"        = "#CD2990",
  "Zeta"           = "#006400",
  "Outgroup"       = "violet",
  "Other"          = "grey60"
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
      )
    )
  
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
    
    out_png <- file.path(
      OUTPUT_DIR,
      paste0(stem, "_", mode, "_bare_fan180.png")
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
    
    tree_grouped <- groupOTU(tree, groups_list)
    
    OPEN_ANGLE <- 180
    
    p_bare <- ggtree(
      tree_grouped,
      layout = "fan",
      open.angle = OPEN_ANGLE,
      aes(color = group),
      size = 0.2
    ) %<+% annot_all +
      scale_color_manual(
        values = CRASSVIRALES_COLOR_SCHEME,
        na.value = CRASSVIRALES_COLOR_SCHEME["Other"]
      )  +
      geom_point2(
        aes(subset = (node == crass_mrca_node)),
        colour = "black",
        size = 2
      ) +
      geom_text2(
        aes(subset = (node == crass_mrca_node), label = paste0("MRCA ", node)),
        colour = "black",
        size = 2,
        nudge_x = 0.02
      ) +
      theme_tree() +
      theme(
        legend.position = "none",
        plot.margin = margin(20, 80, 20, 20)
      )
    
    if (nrow(family_pure_clades) > 0) {
      for (i in seq_len(nrow(family_pure_clades))) {
        fam <- family_pure_clades$family[i]
        node <- family_pure_clades$node[i]
        
        p_bare <- p_bare +
          geom_hilight(
            node = node,
            fill = CRASS_CLADE_FILL[[fam]],
            alpha = 0.25,
            extend = 0.002
          )
      }
    }
    
    # ---- Prophage length bars for fan tree ----
    # ---- Outer annotation rings ----
    # ---- Prophage length bars for fan tree ----
    tree_outer_r <- max(p_bare$data$x, na.rm = TRUE)
    
    genomad_base_r <- tree_outer_r + TREE_TO_GENOMAD_GAP
    genomad_div    <- GENOMAD_DIV
    
    # ---- Prophage length ring starts after class ring ----
    genomad_base_r <- class_base_r + CLASS_RING_WIDTH + CLASS_TO_GENOMAD_GAP
    genomad_div    <- GENOMAD_DIV
    
    # ---- Dotted reference lines for prophage length ----
    genomad_levels <- GENOMAD_LEVELS
    
    genomad_ring_positions <- tibble(
      genomad_len = genomad_levels,
      ring_x      = genomad_base_r + (genomad_levels / genomad_div)
    )
    
    y_range <- range(p_bare$data$y, na.rm = TRUE)
    
    genomad_ring_segments_df <- genomad_ring_positions %>%
      mutate(
        y    = y_range[1],
        yend = y_range[2]
      )
    
    genomad_label_df <- genomad_ring_positions %>%
      mutate(
        y     = y_range[2],
        label = as.character(genomad_len / 1000),   # 50, 100, 150
        x     = ring_x + 0.03
      )
    
    genomad_ring_df <- p_bare$data %>%
      dplyr::filter(
        isTip,
        is_prophage,
        !is.na(bact_genomad_length_num),
        !is.na(y)
      ) %>%
      dplyr::mutate(
        geno_group = if_else(
          !is.na(bact_topology) & bact_topology == "Provirus",
          "integrated_provirus",
          "non_integrated_prophage"
        ),
        geno_x    = genomad_base_r,
        geno_xend = genomad_base_r + bact_genomad_length_num / genomad_div
      )
    
    p_bare <- p_bare +
      # dotted reference lines
      geom_segment(
        data = genomad_ring_segments_df,
        aes(x = ring_x, xend = ring_x, y = y, yend = yend),
        inherit.aes = FALSE,
        colour = "grey50",
        linetype = "dotted",
        linewidth = 0.15,
        lineend = "round"
      ) +
      geom_text(
        data = genomad_label_df,
        aes(x = x, y = y, label = label),
        inherit.aes = FALSE,
        size = 1,
        vjust = -0.2,
        angle = 45,
        colour = "grey30"
      ) +
      ggnewscale::new_scale_color() +
      geom_segment(
        data = genomad_ring_df,
        aes(
          x = geno_x,
          xend = geno_xend,
          y = y,
          yend = y,
          color = geno_group
        ),
        #linewidth = 0.5,
        linewidth = 0.2,
        lineend = "butt",
        inherit.aes = FALSE
      ) +
      scale_color_manual(
        values = c(
          "integrated_provirus"     = "red",
          "non_integrated_prophage" = "blue"
        )
      )
    
    # Prophage ratio annotation
    
    # ---- Prophage ratio annotation ----
    genomad_outer_r <- max(genomad_ring_df$geno_xend, na.rm = TRUE)
    ratio_base_r <- genomad_outer_r + GENOMAD_TO_RATIO_GAP
    
    ratio_levels <- RATIO_LEVELS
    
    ratio_ring_positions <- tibble(
      ratio_value = ratio_levels,
      ring_x      = ratio_base_r + ratio_levels * RATIO_MULT
    )
    
    ratio_ring_segments_df <- ratio_ring_positions %>%
      mutate(
        y    = y_range[1],
        yend = y_range[2]
      )
    
    ratio_label_df <- ratio_ring_positions %>%
      mutate(
        y     = y_range[2],
        label = sprintf("%.2f", ratio_value),
        x     = ring_x + TEXT_EXTRA_GAP
      )
    
    ratio_df <- p_bare$data %>%
      dplyr::filter(
        isTip,
        is_prophage,
        !is.na(bact_prophage_ratio_num),
        !is.na(y)
      ) %>%
      dplyr::mutate(
        ratio_group = if_else(
          !is.na(bact_topology) & bact_topology == "Provirus",
          "integrated_provirus",
          "non_integrated_prophage"
        ),
        ratio_x    = ratio_base_r,
        ratio_xend = ratio_base_r + bact_prophage_ratio_num * RATIO_MULT
      ) %>%
      dplyr::filter(!is.na(ratio_xend))
    
    p_bare <- p_bare +
      # dotted reference lines for ratio
      geom_segment(
        data = ratio_ring_segments_df,
        aes(x = ring_x, xend = ring_x, y = y, yend = yend),
        inherit.aes = FALSE,
        colour = "grey50",
        linetype = "dotted",
        linewidth = 0.15,
        lineend = "round"
      ) +
      geom_text(
        data = ratio_label_df,
        aes(x = x, y = y, label = label),
        inherit.aes = FALSE,
        size = 1,
        vjust = -0.2,
        angle = 45,
        colour = "grey30"
      ) +
      ggnewscale::new_scale_color() +
      geom_segment(
        data = ratio_df,
        aes(
          x = ratio_x,
          xend = ratio_xend,
          y = y,
          yend = y,
          color = ratio_group
        ),
        linewidth = 0.2,
        lineend = "butt",
        inherit.aes = FALSE
      ) +
      scale_color_manual(
        values = c(
          "integrated_provirus"     = "red",
          "non_integrated_prophage" = "blue"
        )
      )
    
    # ---- Dot layer: integration + source ----
    ratio_outer_r <- max(ratio_df$ratio_xend, na.rm = TRUE)
    
    # ---- Phylum ring after ratio ----
    phylum_base_r <- ratio_outer_r + RATIO_TO_PHYLUM_GAP
    
    phylum_ring_df <- p_bare$data %>%
      dplyr::filter(
        isTip,
        is_prophage,
        !is.na(bact_phylum2),
        !is.na(y)
      ) %>%
      dplyr::mutate(
        phyl_x    = phylum_base_r,
        phyl_xend = phylum_base_r + PHYLUM_RING_WIDTH
      )
    
    # ---- Class ring after phylum ----
    class_base_r <- phylum_base_r + PHYLUM_RING_WIDTH + PHYLUM_TO_CLASS_GAP
    
    class_ring_df <- p_bare$data %>%
      dplyr::filter(
        isTip,
        is_prophage,
        !is.na(bact_class2),
        !is.na(y)
      ) %>%
      dplyr::mutate(
        class_x    = class_base_r,
        class_xend = class_base_r + CLASS_RING_WIDTH
      )
    
    p_bare <- p_bare +
      # ---- Phylum ring ----
    ggnewscale::new_scale_color() +
      geom_segment(
        data = phylum_ring_df,
        aes(
          x = phyl_x,
          xend = phyl_xend,
          y = y,
          yend = y,
          colour = bact_phylum2
        ),
        linewidth = 0.5,
        lineend = "butt",
        inherit.aes = FALSE,
        show.legend = FALSE
      ) +
      scale_color_manual(
        values = PHYLUM_COLORS,
        na.value = PHYLUM_COLORS["Other"]
      ) +
      
      # ---- Class ring ----
    ggnewscale::new_scale_color() +
      geom_segment(
        data = class_ring_df,
        aes(
          x = class_x,
          xend = class_xend,
          y = y,
          yend = y,
          colour = bact_class2
        ),
        linewidth = 0.5,
        lineend = "butt",
        inherit.aes = FALSE,
        show.legend = FALSE
      ) +
      scale_color_manual(
        values = CLASS_COLORS,
        na.value = CLASS_COLORS["Other"]
      )
    
    class_outer_r <- max(class_ring_df$class_xend, na.rm = TRUE)
    
    dot_df <- p_bare$data %>%
      dplyr::filter(
        isTip,
        is_prophage,
        !is.na(integration_status)
      ) %>%
      dplyr::mutate(
        dot_x = class_outer_r + CLASS_TO_DOT_GAP,
        dot_fill = if_else(source_type == "isolate", integration_status, "white")
      )
    
    outer_limit <- max(
      c(
        p_bare$data$x,
        genomad_ring_df$geno_xend,
        ratio_df$ratio_xend,
        phylum_ring_df$phyl_xend,
        class_ring_df$class_xend,
        dot_df$dot_x
      ),
      na.rm = TRUE
    )
    
    p_bare <- p_bare + xlim(NA, outer_limit + 0.05)
    
    p_bare <- p_bare +
      ggnewscale::new_scale_color() +
      geom_point(
        data = dot_df,
        aes(
          x = dot_x,
          y = y,
          color = integration_status,
          fill  = dot_fill
        ),
        shape = 21,
        size = 0.1,
        stroke = 0.1,
        inherit.aes = FALSE
      ) +
      scale_color_manual(
        values = c(
          "integrated"     = "red",
          "non_integrated" = "blue"
        )
      ) +
      scale_fill_manual(
        values = c(
          "integrated"     = "red",
          "non_integrated" = "blue",
          "white"          = "white"
        )
      )
    
    outer_limit <- max(
      c(
        p_bare$data$x,
        phylum_ring_df$phyl_xend,
        class_ring_df$class_xend,
        genomad_ring_df$geno_xend,
        ratio_df$ratio_xend,
        dot_df$dot_x
      ),
      na.rm = TRUE
    )
    
    p_bare <- p_bare + xlim(NA, outer_limit + 0.05)
    
    
    cat("[DEBUG] With source info:",
        sum(!is.na(annot_all$source_type)), "\n")
    
    out_png <- file.path(
      OUTPUT_DIR,
      paste0(stem, "_", mode, "_bare_fan180.png")
    )
    
    ggsave(
      out_png,
      p_bare,
      width = 20,
      height = 20,
      units = "cm",
      dpi = 1200
    )
    
    cat("Saved:", out_png, "\n")
  }
}