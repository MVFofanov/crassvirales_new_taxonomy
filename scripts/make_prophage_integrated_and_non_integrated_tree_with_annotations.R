library(ggtree)
library(ggplot2)
library(ape)
library(dplyr)
library(readr)
library(tibble)
library(ggnewscale)  # for multiple colour scales
library(cowplot)

# ---- Paths ----
# WD <- "C:/crassvirales/crassvirales_new_taxonomy/crassvirales_prophages/main_figures/Figure_2_phylogenetic_tree"
# WD <- "C:/crassvirales/crassvirales_new_taxonomy/crassvirales_prophages/maft_iqtree_TerL_analysis/iqtree_trees"
# WD <- "C:/crassvirales/crassvirales_new_taxonomy/crassvirales_prophages/mafft_iqtree_TerL_analysis_reference_best_TerLs/iqtree_trees"
WD <- "C:/crassvirales/crassvirales_new_taxonomy/crassvirales_prophages/mafft_iqtree_TerL_analysis_reference_best_TerLs_without_KC821624_1_100418/iqtree_trees"

# tree_file <- file.path(WD, "TerL.treefile")
# tree_file <- file.path(WD, "TerL_gappy0.5.treefile")

# tree_files <- c(
#   file.path(WD, "TerL_gappy0.5.treefile"),
#   file.path(WD, "TerL_gappy0.7__supports_ufb1000_nm3000.treefile"),
#   file.path(WD, "TerL_gappy0.9__supports_ufb1000_nm3000.treefile"),
#   file.path(WD, "TerL_gappy0.95.treefile"),
#   file.path(WD, "TerL_gappy0.99.treefile"),
#   file.path(WD, "TerL_kpic__supports_ufb1000_nm3000.treefile"),
#   file.path(WD, "TerL_smart-gap__supports_ufb1000_nm3000.treefile")
# )

# tree_files <- c(
#   file.path(WD, "TerL_untrimmed.treefile"),
#   file.path(WD, "TerL_gappy0.5.treefile"),
#   file.path(WD, "TerL_gappy0.6.treefile"),
#   file.path(WD, "TerL_gappy0.7.treefile"),
#   file.path(WD, "TerL_gappy0.8.treefile"),
#   file.path(WD, "TerL_gappy0.9.treefile"),
#   file.path(WD, "TerL_kpic.treefile"),
#   file.path(WD, "TerL_smart-gap.treefile")
# )

# tree_files <- c(
#   file.path(WD, "TerL_gappy0.5.treefile"),
#   file.path(WD, "TerL_gappy0.6.treefile"),
#   file.path(WD, "TerL_gappy0.7.treefile"),
#   file.path(WD, "TerL_gappy0.8.treefile"),
#   file.path(WD, "TerL_gappy0.9.treefile"),
#   file.path(WD, "TerL_kpic.treefile"),
#   file.path(WD, "TerL_smart-gap.treefile")
# )

tree_files <- c(
  file.path(WD, "TerL_gappy0.8.treefile")
)

# annot_file <- file.path(WD, "TerL_tree_protein_taxonomy.tsv")
# annot_file <- file.path(WD, "TerL_gappy0.5_tree_protein_taxonomy.tsv")
# # bact_annot_file <- file.path(WD, "TerL_tree_bacterial_taxonomy_with_genomad_and_lengths.tsv")
# bact_annot_file <- file.path(WD, "TerL_gappy0.5_tree_bacterial_taxonomy_with_genomad_and_lengths.tsv")
# # out_png_circ <- file.path(WD, "Crassvirales_prophages_Fig.2_TerL_tree_circular_annotated.png")
# out_png_circ <- file.path(WD, "TerL_gappy0.5_tree_circular_annotated.png")

# ---- User options ----
# SHOW_TIP_LABELS <- TRUE   # set to FALSE to hide tip labels
SHOW_CRASS_LABELS <- FALSE
#SHOW_PROPHAGE_LABELS  <- FALSE   # <- NEW FLAG
ROOT_AT_OUTGROUP <- FALSE
ROOT_AT_CRASS_MRCA <- TRUE   # root at MRCA of all Crassvirales tips
SHOW_PROPHAGE_LABELS <- TRUE

SHOW_CRASS_MRCA_NODE <- TRUE

# ROOTING_MODES <- c(
#   "unrooted",
#   "outgroup",
#   "crass_mrca",
#   "crass_mrca_parent",
#   "crass_plus_integrated_mrca"   # NEW
# )

ROOTING_MODES <- c(
  "outgroup"
)

# ---- Color scheme for viral families + outgroup ----
CRASSVIRALES_COLOR_SCHEME <- c(
  "Intestiviridae" = "#EE3B3B",
  "Crevaviridae"   = "#EE9A00",
  "Suoliviridae"   = "#4169E1",
  "Steigviridae"   = "#00CED1",
  "Epsilon"        = "#CD2990",
  "Zeta"           = "#006400",
  "Outgroup"       = "violet"
)
crass_families <- setdiff(names(CRASSVIRALES_COLOR_SCHEME), "Outgroup")

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

# ---------- helper functions ----------
get_nodes_in_clade <- function(tree, node) {
  if (is.null(node) || is.na(node)) return(integer(0))
  unique(c(node, get_descendants(tree, node)))
}


get_parent_node <- function(tree, node) {
  p <- tree$edge[tree$edge[, 2] == node, 1]
  if (length(p) == 0) return(NA_integer_)  # node is root
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

get_ancestors <- function(tree, node) {
  parents <- integer(0)
  cur <- node
  repeat {
    p <- tree$edge[tree$edge[, 2] == cur, 1]
    if (length(p) == 0) break
    parents <- c(parents, p)
    cur <- p
  }
  parents
}

# Return maximal (non-nested) monophyletic clade nodes for a given family
get_family_clade_nodes <- function(tree, tip_family, fam, min_tips = 2) {
  Ntip <- length(tree$tip.label)
  internal_nodes <- (Ntip + 1):(Ntip + tree$Nnode)
  
  # candidate nodes: all descendant tips are fam, and clade has >= min_tips tips
  cand <- integer(0)
  for (nd in internal_nodes) {
    desc_tip_ids <- get_tip_descendants(tree, nd)
    labs <- tree$tip.label[desc_tip_ids]
    fams <- unique(tip_family[labs])
    
    if (length(labs) >= min_tips && length(fams) == 1 && fams == fam) {
      cand <- c(cand, nd)
    }
  }
  
  if (length(cand) == 0) return(integer(0))
  
  # keep only maximal clades (drop those nested inside another candidate)
  keep <- rep(TRUE, length(cand))
  for (i in seq_along(cand)) {
    if (any(get_ancestors(tree, cand[i]) %in% cand)) keep[i] <- FALSE
  }
  cand[keep]
}

# Maximal (non-nested) clades where:
# - descendant tips contain Crassvirales members from exactly ONE family (fam)
# - BUT can also contain non-Crass tips (unknown / prophage / outgroup etc.)
# "Big" is controlled by min_crass (number of Crass tips of fam) and/or min_total.
get_family_clade_nodes_allow_noncrass <- function(
    tree, tip_family, crass_families, fam,
    min_crass = 2,      # <- tweak (e.g., 2, 5, 10)
    min_total = 0       # <- optional: require total tips too
) {
  Ntip <- length(tree$tip.label)
  internal_nodes <- (Ntip + 1):(Ntip + tree$Nnode)
  
  cand <- integer(0)
  
  for (nd in internal_nodes) {
    desc_tip_ids <- get_tip_descendants(tree, nd)
    labs <- tree$tip.label[desc_tip_ids]
    fams_all <- tip_family[labs]
    
    # which Crass families are present in this node?
    crass_present <- unique(fams_all[fams_all %in% crass_families])
    
    # must contain exactly one Crass family, and it must be "fam"
    if (length(crass_present) == 1 && crass_present == fam) {
      
      # "big" criteria: number of Crass tips of that fam, and optional total tips
      n_crass_fam <- sum(fams_all == fam, na.rm = TRUE)
      n_total <- length(labs)
      
      if (n_crass_fam >= min_crass && n_total >= min_total) {
        cand <- c(cand, nd)
      }
    }
  }
  
  if (length(cand) == 0) return(integer(0))
  
  # keep only maximal nodes: drop those nested inside another candidate of same family
  keep <- rep(TRUE, length(cand))
  for (i in seq_along(cand)) {
    if (any(get_ancestors(tree, cand[i]) %in% cand)) keep[i] <- FALSE
  }
  cand[keep]
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
    
    Ntip <- length(tree$tip.label)
    root_node <- Ntip + 1
    
    cat("is.rooted:", ape::is.rooted(tree), "\n")
    cat("root node:", root_node, "\n")
    cat("crass_mrca_node:", crass_mrca_node, "\n")
    cat("MRCA equals root?", crass_mrca_node == root_node, "\n")
    
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
    
    # Prefer real outgroup if available
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
  
  if (mode == "crass_mrca_parent") {
    
    crass_tips <- annot_all %>%
      filter(family %in% crass_families) %>%
      pull(label) %>%
      intersect(tree$tip.label)
    
    if (length(crass_tips) < 2) {
      stop("[ROOT] <2 Crass tips for MRCA-parent rooting in: ", stem)
    }
    
    mrca <- ape::getMRCA(tree, crass_tips)
    if (is.null(mrca) || is.na(mrca)) {
      stop("[ROOT] Could not compute Crass MRCA in: ", stem)
    }
    
    parent <- get_parent_node(tree, mrca)
    if (is.na(parent)) {
      warning("[ROOT] Crass MRCA is already the root; returning rooted-as-is for: ", stem)
      return(tree)
    }
    
    message("[ROOT] Rooting at parent of Crass MRCA. MRCA=", mrca, " parent=", parent)
    
    tree <- ape::root(tree, node = parent, resolve.root = TRUE)
    tree <- ape::reorder.phylo(tree, "cladewise")
    return(tree)
  }
  
  if (mode == "crass_plus_integrated_mrca") {
    
    # 1) reference Crass tips
    crass_tips <- annot_all %>%
      dplyr::filter(family %in% crass_families) %>%
      dplyr::pull(label) %>%
      intersect(tree$tip.label)
    
    # 2) tips that have integrated prophage label
    integ_tips <- annot_all %>%
      dplyr::filter(!is.na(integrated_prophage_label), integrated_prophage_label != "") %>%
      dplyr::pull(label) %>%
      intersect(tree$tip.label)
    
    target_tips <- unique(c(crass_tips, integ_tips))
    
    if (length(target_tips) < 2) {
      stop("[ROOT] <2 tips for crass+integrated MRCA rooting in: ", stem,
           " (crass=", length(crass_tips), ", integrated=", length(integ_tips), ")")
    }
    
    # pick an outgroup tip OUTSIDE the target clade
    non_target_tips <- setdiff(tree$tip.label, target_tips)
    if (length(non_target_tips) < 1) {
      stop("[ROOT] No tips outside (crass+integrated) set for rooting in: ", stem)
    }
    
    # Prefer your known outgroup tip if it is outside target_tips
    if (length(outgroup_tips) == 1 && outgroup_tips %in% non_target_tips) {
      out_tip <- outgroup_tips
    } else {
      out_tip <- non_target_tips[1]
    }
    
    message("[ROOT] Rooting at MRCA(crass refs + integrated-labeled) using outgroup tip: ", out_tip,
            " | n_crass=", length(crass_tips), " n_integrated=", length(integ_tips),
            " n_target=", length(target_tips))
    
    tree <- ape::root(tree, outgroup = out_tip, resolve.root = TRUE)
    tree <- ape::reorder.phylo(tree, "cladewise")
    return(tree)
  }
  
  stop("Unknown rooting mode: ", mode)
}

find_biggest_noncrass_clade <- function(tree, crass_tips, forbidden_nodes = integer(0)) {
  Ntip <- length(tree$tip.label)
  internal_nodes <- (Ntip + 1):(Ntip + tree$Nnode)
  crass_set <- unique(crass_tips)
  forbidden_nodes <- unique(forbidden_nodes)
  
  best_node <- NA_integer_
  best_size <- -1L
  
  for (nd in internal_nodes) {
    if (nd %in% forbidden_nodes) next
    
    tip_ids <- get_tip_descendants(tree, nd)
    labs <- tree$tip.label[tip_ids]
    
    # must contain NO Crass tips
    if (any(labs %in% crass_set)) next
    
    n <- length(labs)
    if (n > best_size) {
      best_size <- n
      best_node <- nd
    }
  }
  
  list(node = best_node, n_tips = best_size)
}

# leaf tip label -> accession (handles both "|" and "_XXXX_CDS_####" suffixes)
leaf_to_accession <- function(x) {
  x %>%
    stringr::str_replace("\\|.*$", "") %>%                  # drop everything after "|"
    stringr::str_replace("_[A-Z]{8}_CDS_[0-9]+$", "") %>%   # drop CDS suffix if present
    stringr::str_replace("\\r$", "") %>%                    # in case of CRLF
    stringr::str_trim()
}

# pairs$prophage value "{accession}_{start}-{end}" -> accession
prophage_fragment_to_accession <- function(x) {
  x %>%
    stringr::str_replace("_[0-9]+-[0-9]+$", "") %>%         # drop _start-end
    stringr::str_replace("\\r$", "") %>%
    stringr::str_trim()
}

# Pick an explicit outgroup tip that is NOT in target_tips
pick_outgroup_tip <- function(tree, target_tips, outgroup_tips) {
  non_target <- setdiff(tree$tip.label, target_tips)
  if (length(non_target) < 1) stop("[MRCA] No tip outside target set; cannot pick outgroup.")
  
  if (length(outgroup_tips) == 1 && outgroup_tips %in% non_target) {
    return(outgroup_tips)
  }
  non_target[1]
}

# Return a tibble of prophage tips (integrated+non-integrated) inside MRCA(target_tips)
get_prophage_in_target_mrca <- function(tree, annot_all, target_tips) {
  target_tips <- intersect(target_tips, tree$tip.label)
  if (length(target_tips) < 2) {
    return(list(mrca_node = NA_integer_, tips = character(0), df = tibble()))
  }
  
  mrca_node <- ape::getMRCA(tree, target_tips)
  if (is.null(mrca_node) || is.na(mrca_node)) {
    return(list(mrca_node = NA_integer_, tips = character(0), df = tibble()))
  }
  
  tip_ids <- get_tip_descendants(tree, mrca_node)
  tip_labels <- tree$tip.label[tip_ids]
  
  df <- annot_all %>%
    dplyr::filter(label %in% tip_labels, is_prophage) %>%   # <- integrated + non-integrated
    dplyr::arrange(label)
  
  list(mrca_node = as.integer(mrca_node), tips = tip_labels, df = df)
}

# Save both TSV (full rows) + TXT (labels only)
save_prophage_lists <- function(df, out_tsv, out_txt) {
  readr::write_tsv(df, out_tsv)
  writeLines(df$label, out_txt)
  cat("[OK] Saved:", out_tsv, " (n=", nrow(df), ")\n", sep = "")
  cat("[OK] Saved:", out_txt, "\n")
}

extract_prophage_id_from_label <- function(x) {
  x %>%
    stringr::str_replace("_[A-Z0-9]+_CDS_[0-9]+$", "") %>%  # remove protein suffix
    stringr::str_replace("\\r$", "") %>%
    stringr::str_trim()
}



# # ---- Build mapping: accession -> prophage_label (ONCE, before loops) ----
# PAIRS_TSV <- file.path(WD, "pairs_prophage_with_flanks_vs_bacterial.tsv")
# 
# pairs <- readr::read_tsv(PAIRS_TSV, show_col_types = FALSE) %>%
#   transmute(
#     prophage = as.character(prophage),
#     prophage_label = as.character(prophage_label)
#   ) %>%
#   filter(!is.na(prophage), prophage != "")
# 
# dups <- pairs %>% count(prophage) %>% filter(n > 1)
# if (nrow(dups) > 0) {
#   warning("Duplicate prophage values in pairs.tsv: ",
#           paste(dups$prophage, collapse = ", "))
# }
# 
# cat("[INFO] pairs rows:", nrow(pairs), "\n")
# cat("[INFO] unique prophage values in pairs:", dplyr::n_distinct(pairs$prophage), "\n")

# ---- Read integrated prophage candidate list ----
CANDIDATE_LIST_TXT <- file.path(WD, "Crassvirales_integrated_prophage_candidate_list.txt")

candidate_ids <- readLines(CANDIDATE_LIST_TXT, warn = FALSE) %>%
  stringr::str_replace("\\r$", "") %>%
  stringr::str_trim()

candidate_ids <- candidate_ids[candidate_ids != ""]
candidate_ids <- unique(candidate_ids)

cat("[INFO] candidate prophage ids in list:", length(candidate_ids), "\n")



for (tree_file in tree_files) {
  
  message("==== Processing tree file: ", tree_file)
  stem <- tools::file_path_sans_ext(basename(tree_file))
  
  tree_raw <- read.tree(tree_file)
  stem <- tools::file_path_sans_ext(basename(tree_file))
  
  annot_file      <- file.path(WD, paste0(stem, "_tree_protein_taxonomy.tsv"))
  bact_annot_file <- file.path(WD, paste0(stem, "_tree_bacterial_taxonomy_with_genomad_and_lengths.tsv"))
  
  annot_raw <- read_tsv(annot_file, show_col_types = FALSE)
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
  
  # --- JOIN BACTERIAL TABLE ONCE (IMPORTANT) ---
  bact_annot_raw <- read_tsv(bact_annot_file, show_col_types = FALSE)
  bact_annot <- bact_annot_raw %>%
    rename(label = leaf_label) %>%
    rename_with(~ ifelse(startsWith(.x, "bact_"), .x, paste0("bact_", .x)), -label)
  
  
  annot_all <- annot_all %>% left_join(bact_annot, by = "label")
  
  # compute flags ONCE
  annot_all <- annot_all %>%
    mutate(
      is_crassvirales = genome_id != "unknown" & !(genome_id %in% outgroup_genomes),
      is_prophage     = !is.na(bact_domain) & bact_domain == "Bacteria",
      
      bact_genomad_length_num = suppressWarnings(as.numeric(bact_genomad_length)),
      bact_prophage_ratio_num = suppressWarnings(as.numeric(bact_prophage_ratio))
    )
  
  annot_all <- annot_all %>%
    mutate(
      prophage_id = dplyr::if_else(
        is_prophage,
        extract_prophage_id_from_label(label),
        NA_character_
      )
    )
  

  
  # mark whether extracted prophage_id is in candidate list
  annot_all <- annot_all %>%
    mutate(
      integrated_prophage_label = dplyr::if_else(
        !is.na(prophage_id) & prophage_id %in% candidate_ids,
        prophage_id,
        NA_character_
      )
    )
  
  n_prophage_leaves_with_id <- annot_all %>%
    dplyr::filter(is_prophage, !is.na(prophage_id), prophage_id != "") %>%
    nrow()
  
  unique_tree_prophage_ids <- annot_all %>%
    dplyr::filter(is_prophage, !is.na(prophage_id), prophage_id != "") %>%
    dplyr::pull(prophage_id) %>%
    unique()
  
  found_candidate_ids <- intersect(candidate_ids, unique_tree_prophage_ids)
  not_found_candidate_ids <- setdiff(candidate_ids, unique_tree_prophage_ids)
  
  cat("[INFO] prophage leaves given prophage_id:", n_prophage_leaves_with_id, "\n")
  cat("[INFO] unique prophage_ids extracted from tree leaves:", length(unique_tree_prophage_ids), "\n")
  cat("[INFO] candidate prophage ids found in tree:", length(found_candidate_ids), " / ", length(candidate_ids), "\n", sep = "")
  cat("[INFO] candidate prophage ids NOT found in tree:", length(not_found_candidate_ids), "\n")
  
  if (length(not_found_candidate_ids) > 0) {
    cat("[INFO] Candidate prophage ids not found in tree:\n")
    print(not_found_candidate_ids)
  }
  

  cat("[INFO] Tips with integrated_prophage_label:",
      sum(!is.na(annot_all$integrated_prophage_label) & annot_all$integrated_prophage_label != ""),
      "\n")
  
  cat("Columns in annot_raw:\n")
  print(colnames(annot_raw))

  cat("First 20 tree labels:\n")
  print(head(annot_all$label, 20))
  
  # annot_all %>%
  #   dplyr::filter(grepl("NZ_JAHXKQ010000025", label) | grepl("NZ_JAHXKQ010000025", genome_id)) %>%
  #   dplyr::select(label, genome_id, protein_id) %>%
  #   print(n = 20)
  # 
  # annot_all %>%
  #   dplyr::filter(genome_id %in% pairs$prophage) %>%
  #   dplyr::select(label, genome_id, protein_id, family) %>%
  #   print(n = 50)

  
  
  # # safety: ensure expected columns exist
  # needed <- c("bact_domain","bact_topology","bact_genomad_length","bact_prophage_ratio","bact_phylum","bact_class")
  # missing <- setdiff(needed, names(annot_all))
  # if (length(missing) > 0) {
  #   warning("Missing bacterial columns in ", stem, ": ", paste(missing, collapse=", "))
  #   for (nm in missing) annot_all[[nm]] <- NA
  # }
  # 
  # # compute these once (same for all modes)
  # annot_all <- annot_all %>%
  #   mutate(
  #     is_crassvirales = genome_id != "unknown" & !(genome_id %in% outgroup_genomes),
  #     is_prophage     = !is.na(bact_domain) & bact_domain == "Bacteria",
  #     bact_genomad_length_num = suppressWarnings(as.numeric(bact_genomad_length))
  #   )
  
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
    
    # ---------------------------------------------------------
    # MRCA of the 16 integrated prophage candidate tips only
    # ---------------------------------------------------------
    candidate_prophage_tips <- annot_all %>%
      dplyr::filter(
        is_prophage,
        !is.na(prophage_id),
        prophage_id %in% candidate_ids
      ) %>%
      dplyr::pull(label) %>%
      unique() %>%
      intersect(tree$tip.label)
    
    cat("[INFO] candidate_prophage_tips found in tree:", length(candidate_prophage_tips), "\n")
    
    if (length(candidate_prophage_tips) < 2) {
      warning("[MRCA] Fewer than 2 candidate prophage tips found in tree for: ", stem)
    } else {
      candidate_mrca_node <- ape::getMRCA(tree, candidate_prophage_tips)
      
      if (is.null(candidate_mrca_node) || is.na(candidate_mrca_node)) {
        warning("[MRCA] Could not compute MRCA of candidate prophage tips for: ", stem)
      } else {
        candidate_mrca_tip_ids <- get_tip_descendants(tree, candidate_mrca_node)
        candidate_mrca_tip_labels <- tree$tip.label[candidate_mrca_tip_ids]
        
        # keep ONLY is_prophage leaves from that MRCA
        candidate_mrca_prophage_df <- annot_all %>%
          dplyr::filter(
            label %in% candidate_mrca_tip_labels,
            is_prophage
          ) %>%
          dplyr::arrange(label)
        
        candidate_mrca_prophage_labels <- candidate_mrca_prophage_df$label
        
        out_candidate_mrca_txt <- file.path(
          WD,
          paste0(stem, "_", mode, "_is_prophage_leaves_in_MRCA_of_16_integrated_candidates.txt")
        )
        writeLines(candidate_mrca_prophage_labels, out_candidate_mrca_txt)
        
        out_candidate_mrca_tsv <- file.path(
          WD,
          paste0(stem, "_", mode, "_is_prophage_leaves_in_MRCA_of_16_integrated_candidates.tsv")
        )
        readr::write_tsv(candidate_mrca_prophage_df, out_candidate_mrca_tsv)
        
        out_candidate_tip_txt <- file.path(
          WD,
          paste0(stem, "_", mode, "_16_integrated_candidate_tip_labels.txt")
        )
        writeLines(candidate_prophage_tips, out_candidate_tip_txt)
        
        cat("[OK] candidate_mrca_node =", candidate_mrca_node,
            " | descendant prophage leaves =", nrow(candidate_mrca_prophage_df), "\n")
        cat("[OK] Saved:", out_candidate_mrca_txt, "\n")
        cat("[OK] Saved:", out_candidate_mrca_tsv, "\n")
        cat("[OK] Saved:", out_candidate_tip_txt, "\n")
      }
    }
    
    # -----------------------------
    # MRCA(crass refs + integrated-labeled) prophage tips
    # Save BEFORE and AFTER rooting
    # -----------------------------
    
    # define target set once
    crass_ref_tips <- annot_all %>%
      dplyr::filter(family %in% crass_families) %>%
      dplyr::pull(label) %>%
      intersect(tree_raw$tip.label)
    
    integrated_labeled_tips <- annot_all %>%
      dplyr::filter(!is.na(integrated_prophage_label), integrated_prophage_label != "") %>%
      dplyr::pull(label) %>%
      intersect(tree_raw$tip.label)
    
    target_tips <- unique(c(crass_ref_tips, integrated_labeled_tips))
    
    if (length(target_tips) >= 2) {
      
      # Choose ONE explicit outgroup tip outside the target set (consistent definition)
      out_tip <- pick_outgroup_tip(tree_raw, target_tips, outgroup_tips)
      
      # ---- "BEFORE rooting": use a TEMP rooted copy of tree_raw (do not reorder labels etc.)
      tree_pre <- ape::root(tree_raw, outgroup = out_tip, resolve.root = TRUE)
      tree_pre <- ape::reorder.phylo(tree_pre, "cladewise")
      
      pre_res <- get_prophage_in_target_mrca(tree_pre, annot_all, target_tips)
      
      out_pre_tsv <- file.path(WD, paste0(stem, "_", mode, "_crassPlusInt_MRCA_prophages_PRE.tsv"))
      out_pre_txt <- file.path(WD, paste0(stem, "_", mode, "_crassPlusInt_MRCA_prophages_PRE.txt"))
      save_prophage_lists(pre_res$df, out_pre_tsv, out_pre_txt)
      
      cat("[INFO] PRE MRCA node =", pre_res$mrca_node,
          " | out_tip =", out_tip,
          " | n_target =", length(target_tips), "\n", sep = "")
      
      # ---- "AFTER rooting": use the final rooted tree for this mode (your `tree`)
      post_res <- get_prophage_in_target_mrca(tree, annot_all, target_tips)
      
      out_post_tsv <- file.path(WD, paste0(stem, "_", mode, "_crassPlusInt_MRCA_prophages_POST.tsv"))
      out_post_txt <- file.path(WD, paste0(stem, "_", mode, "_crassPlusInt_MRCA_prophages_POST.txt"))
      save_prophage_lists(post_res$df, out_post_tsv, out_post_txt)
      
      cat("[INFO] POST MRCA node =", post_res$mrca_node,
          " | n_target =", length(target_tips), "\n", sep = "")
      
    } else {
      cat("[WARN] Not enough target tips for Crass+Integrated MRCA (n=", length(target_tips), ")\n", sep = "")
    }
    
    
    out_png_circ <- file.path(WD, paste0(stem, "_tree_", mode, "_circular_annotated.png"))
    
  
  # # ---- Optional rooting at outgroup tip ----
  # if (ROOT_AT_OUTGROUP && !ROOT_AT_CRASS_MRCA) {
  #   if (length(outgroup_tips) == 0) stop("ROOT_AT_OUTGROUP=TRUE but no outgroup tip found for: ", stem)
  #   if (length(outgroup_tips) > 1) stop("Multiple outgroup tips found for: ", stem, " -> ", paste(outgroup_tips, collapse=", "))
  #   if (!(outgroup_tips %in% tree$tip.label)) stop("Outgroup tip not in tree for: ", stem)
  #   
  #   tree <- root(tree, outgroup = outgroup_tips, resolve.root = TRUE)
  #   tree <- reorder.phylo(tree, order = "cladewise")
  # }
  # 
  # # ---- Optional rooting at MRCA of all Crassvirales tips ----
  # # ---- Optional rooting at MRCA of all Crassvirales tips (robust) ----
  # if (ROOT_AT_CRASS_MRCA) {
  #   
  #   crass_tips <- annot_all %>%
  #     dplyr::filter(family %in% crass_families) %>%
  #     dplyr::pull(label) %>%
  #     intersect(tree$tip.label)
  #   
  #   if (length(crass_tips) < 2) {
  #     stop("ROOT_AT_CRASS_MRCA=TRUE but fewer than 2 Crassvirales tips found in tree for: ", stem)
  #   }
  #   
  #   noncrass_tips <- setdiff(tree$tip.label, crass_tips)
  #   if (length(noncrass_tips) < 1) {
  #     stop("ROOT_AT_CRASS_MRCA=TRUE but no non-Crass tips exist in tree for: ", stem)
  #   }
  #   
  #   # Compute MRCA on current tree
  #   crass_mrca_node <- ape::getMRCA(tree, crass_tips)
  #   if (is.null(crass_mrca_node) || is.na(crass_mrca_node)) {
  #     stop("Could not compute Crassvirales MRCA for: ", stem)
  #   }
  #   cat("[INFO] Crassvirales MRCA node id (pre-root) =", crass_mrca_node, "\n")
  #   
  #   # --- Choose ONE non-Crass tip as explicit outgroup (always monophyletic) ---
  #   # Prefer your known outgroup genome if present; otherwise take the first non-Crass tip
  #   if (length(outgroup_tips) == 1 && outgroup_tips %in% noncrass_tips) {
  #     out_tip <- outgroup_tips
  #   } else {
  #     out_tip <- noncrass_tips[1]
  #   }
  #   cat("[INFO] Using explicit outgroup tip =", out_tip, "\n")
  #   
  #   # Root on that single tip
  #   tree <- ape::root(tree, outgroup = out_tip, resolve.root = TRUE)
  #   tree <- ape::reorder.phylo(tree, order = "cladewise")
  #   
  #   # --- Now rotate so that Crass clade is the "main" side (optional but makes consistent plots) ---
  #   # After rooting, Crass MRCA should be definable; we just recompute it for plotting.
  #   crass_mrca_node <- ape::getMRCA(tree, intersect(crass_tips, tree$tip.label))
  #   cat("[INFO] Crassvirales MRCA node id (post-root) =", crass_mrca_node, "\n")
  # }
  
  
  
  # ---- MRCA of all Crassvirales tips (node id) ----
  crass_tips <- annot_all %>%
    filter(family %in% crass_families) %>%
    pull(label) %>%
    intersect(tree$tip.label)

  crass_mrca_node <- ape::getMRCA(tree, crass_tips)
  
  # ---- Save prophage tips inside Crass MRCA clade ----
  
  # all descendant tips (labels) of the Crass MRCA node
  crass_mrca_tip_ids <- get_tip_descendants(tree, crass_mrca_node)
  crass_mrca_tip_labels <- tree$tip.label[crass_mrca_tip_ids]
  
  # keep only prophage tips (is_prophage == TRUE)
  prophage_in_crass_mrca <- annot_all %>%
    dplyr::filter(
      label %in% crass_mrca_tip_labels,
      is_prophage
    ) %>%
    dplyr::arrange(label)
  
  out_prophage_list_tsv <- file.path(
    WD,
    paste0(stem, "_", mode, "_prophage_tips_in_crass_mrca.tsv")
  )
  
  readr::write_tsv(prophage_in_crass_mrca, out_prophage_list_tsv)
  cat("[OK] Saved prophage tips in Crass MRCA clade:", out_prophage_list_tsv,
      " (n=", nrow(prophage_in_crass_mrca), ")\n", sep = "")
  
  out_prophage_list_txt <- file.path(
    WD,
    paste0(stem, "_", mode, "_prophage_tips_in_crass_mrca.txt")
  )
  
  writeLines(prophage_in_crass_mrca$label, out_prophage_list_txt)
  cat("[OK] Saved labels:", out_prophage_list_txt, "\n")

  crass_nodes <- get_nodes_in_clade(tree, crass_mrca_node)
  
  biggest_nc <- find_biggest_noncrass_clade(tree, crass_tips, forbidden_nodes = crass_nodes)
  biggest_noncrass_node <- biggest_nc$node
  biggest_noncrass_size <- biggest_nc$n_tips
  
  crass_mrca_parent_node <- get_parent_node(tree, crass_mrca_node)

  
  # label to show on the tree
  crass_parent_label <- if (is.na(crass_parent_node)) {
    "Parent of Crass MRCA: ROOT"
  } else {
    paste0("Parent of Crass MRCA: ", crass_parent_node)
  }
  
  cat("[INFO] crass_mrca_node =", crass_mrca_node,
      " parent =", crass_mrca_parent_node, "\n")
  
  
  if (length(crass_tips) < 2) {
    stop("Cannot compute Crassvirales MRCA: <2 Crassvirales tips in tree for: ", stem)
  }
  
  # ---- MRCA of (Crass references + integrated-labeled tips) ----
  crass_ref_tips <- annot_all %>%
    dplyr::filter(family %in% crass_families) %>%
    dplyr::pull(label) %>%
    intersect(tree$tip.label)
  
  integrated_labeled_tips <- annot_all %>%
    dplyr::filter(!is.na(integrated_prophage_label), integrated_prophage_label != "") %>%
    dplyr::pull(label) %>%
    intersect(tree$tip.label)
  
  target_tips <- unique(c(crass_ref_tips, integrated_labeled_tips))
  
  crass_plus_integrated_mrca_node <- if (length(target_tips) >= 2) {
    ape::getMRCA(tree, target_tips)
  } else {
    NA_integer_
  }
  
  candidate_only_mrca_node <- if (length(candidate_prophage_tips) >= 2) {
    ape::getMRCA(tree, candidate_prophage_tips)
  } else {
    NA_integer_
  }
  
  cat("[INFO] candidate_only_mrca_node =", candidate_only_mrca_node,
      " (n_candidate_tips=", length(candidate_prophage_tips), ")\n")
  
  cat("[INFO] crass_plus_integrated_mrca_node =", crass_plus_integrated_mrca_node,
      " (n_target=", length(target_tips), ")\n")
  
  
  # crass_mrca_node <- ape::getMRCA(tree, crass_tips)
  # 
  # if (is.null(crass_mrca_node) || is.na(crass_mrca_node)) {
  #   stop("Could not compute Crassvirales MRCA for: ", stem)
  # }
  
  # cat("[INFO] Crassvirales MRCA node id =", crass_mrca_node, "\n")
  
  
  # ---- Bacterial annotation ----
  # bact_annot_raw <- read_tsv(bact_annot_file, show_col_types = FALSE)
  
  # annot_file      <- file.path(WD, paste0(stem, "_tree_protein_taxonomy.tsv"))
  # bact_annot_file <- file.path(WD, paste0(stem, "_tree_bacterial_taxonomy_with_genomad_and_lengths.tsv"))
  #   
  # if (!file.exists(annot_file)) stop("Missing annot_file: ", annot_file)
  # if (!file.exists(bact_annot_file)) stop("Missing bact_annot_file: ", bact_annot_file)
  
  # after: annot_all <- annot_all %>% left_join(bact_annot, by="label")
  # needed <- c("bact_domain","bact_topology","bact_genomad_length","bact_prophage_ratio","bact_phylum","bact_class")
  # missing <- setdiff(needed, names(annot_all))
  # if (length(missing) > 0) {
  #   warning("Missing bacterial columns in ", stem, ": ", paste(missing, collapse=", "))
  #   for (nm in missing) annot_all[[nm]] <- NA
  # }

  
  # bact_annot <- bact_annot_raw %>%
  #   rename(label = leaf_label) %>%
  #   rename_with(~ paste0("bact_", .x), -label)
  # 
  # annot_all <- annot_all %>%
  #   left_join(bact_annot, by = "label") %>%
  #   mutate(
  #     is_crassvirales = genome_id != "unknown" & !(genome_id %in% outgroup_genomes),
  #     is_prophage     = !is.na(bact_domain) & bact_domain == "Bacteria",
  #     bact_genomad_length_num = suppressWarnings(as.numeric(bact_genomad_length))
  #   )
  
  # =========================
  # ---- EXTRA FIGURE: unrooted tree with Crassvirales reference clades in red ----
  # Place this block right AFTER you create annot_all with is_crassvirales / is_prophage
  # =========================
  
  # out_png_unroot <- file.path(WD, "TerL_tree_unrooted_crass_reference_branches.png")
  # 
  # # Tip labels that are Crassvirales references (not unknown, not outgroup)
  # ref_tips <- annot_all %>%
  #   filter(is_crassvirales) %>%
  #   pull(label)
  # 
  # # Convert labels -> tip indices (ape numbering: tips are 1..Ntip)
  # ref_tip_ids <- match(ref_tips, tree$tip.label)
  # ref_tip_ids <- ref_tip_ids[!is.na(ref_tip_ids)]
  # 
  # # All nodes on paths from each reference tip to the root
  # nodes_to_red <- unique(c(
  #   ref_tip_ids,
  #   unlist(lapply(ref_tip_ids, function(tid) get_ancestors(tree, tid)))
  # ))
  # 
  # # Create node-wise colour annotation for ggtree
  # Ntip <- length(tree$tip.label)
  # n_nodes_total <- Ntip + tree$Nnode
  # 
  # branch_df <- tibble(
  #   node = 1:n_nodes_total,
  #   branch_col = if_else(node %in% nodes_to_red, "Crassvirales_ref_path", "Other")
  # )
  # 
  # p_unroot <- ggtree(tree, layout = "unrooted") %<+% branch_df +
  #   geom_tree(aes(color = branch_col), linewidth = 0.25) +
  #   scale_color_manual(
  #     values = c(
  #       "Crassvirales_ref_path" = "red",
  #       "Other" = "black"
  #     ),
  #     guide = "none"
  #   ) +
  #   theme_tree() +
  #   theme(plot.margin = margin(5, 5, 5, 5))
  # 
  # ggsave(out_png_unroot, p_unroot, width = 16, height = 16, units = "cm", dpi = 1200)
  # cat("Saved:", out_png_unroot, "\n")
  
  
  # ---- Phylum + Class categories ----
  annot_all <- annot_all %>%
    mutate(
      bact_phylum2 = case_when(
        bact_phylum %in% names(PHYLUM_COLORS) ~ bact_phylum,
        TRUE                                   ~ "Other"
      ),
      bact_class2 = case_when(
        bact_class %in% names(CLASS_COLORS) ~ bact_class,
        TRUE                                ~ "Other"
      )
    )
  
  # ---- Mark families in tree ----
  groups_list <- lapply(crass_families, function(fam) {
    annot_all %>%
      filter(family == fam) %>%
      pull(label)
  })
  names(groups_list) <- crass_families
  groups_list[["Outgroup"]] <- outgroup_tips
  
  tree_grouped <- groupOTU(tree, groups_list)
  
  # ---- Collapse Crassvirales-only clades ----
  tip_family <- annot_all$family
  names(tip_family) <- annot_all$label
  tip_family <- tip_family[tree$tip.label]  # align to tree tips
  
  # ---- Find highlight clades per family (allow non-Crass tips inside) ----
  MIN_CRASS_TIPS <- 5   # <-- define what "big" means (change to 2/10/etc.)
  MIN_TOTAL_TIPS <- 0   # <-- optional (keep 0 if you don't care)
  
  fam_clades_list <- lapply(crass_families, function(fam) {
    nodes <- get_family_clade_nodes_allow_noncrass(
      tree = tree,
      tip_family = tip_family,
      crass_families = crass_families,
      fam = fam,
      min_crass = MIN_CRASS_TIPS,
      min_total = MIN_TOTAL_TIPS
    )
    if (length(nodes) == 0) return(NULL)
    tibble(node = nodes, fam = fam)
  })
  
  fam_clades_df <- dplyr::bind_rows(fam_clades_list)
  
  
  
  Ntip <- length(tree$tip.label)
  all_nodes <- (Ntip + 1):(Ntip + tree$Nnode)
  
  is_coll <- logical(Ntip + tree$Nnode)
  names(is_coll) <- as.character(1:(Ntip + tree$Nnode))
  
  for (nd in all_nodes) {
    desc <- get_tip_descendants(tree, nd)
    labs <- tree$tip.label[desc]
    fams <- unique(tip_family[labs])
    if (all(fams %in% crass_families) && length(fams) == 1) {
      is_coll[as.character(nd)] <- TRUE
    }
  }
  
  cand <- as.integer(names(is_coll)[is_coll])
  to_drop <- logical(length(cand))
  for (i in seq_along(cand)) {
    if (any(get_ancestors(tree, cand[i]) %in% cand)) {
      to_drop[i] <- TRUE
    }
  }
  collapse_nodes <- cand[!to_drop]
  
  cat("Collapsing", length(collapse_nodes), "Clades\n")
  
  # ---- Compute node leaf counts ----
  n_nodes_total <- Ntip + tree$Nnode
  n_leaves <- integer(n_nodes_total)
  names(n_leaves) <- as.character(1:n_nodes_total)
  
  for (nd in 1:n_nodes_total) {
    if (nd <= Ntip) {
      n_leaves[nd] <- 1L
    } else {
      n_leaves[nd] <- length(get_tip_descendants(tree, nd))
    }
  }
  
  node_sizes_df <- tibble(
    node    = as.integer(names(n_leaves)),
    bar_len = log10(n_leaves) + 0.1
  )
  
  # ---- Build circular tree ----
  
  # p_circ <- ggtree(
  #   tree_grouped,
  #   layout = "circular",
  #   aes(color = group),
  #   size = 0.2
  # ) %<+% annot_all
  # 
  
  OPEN_ANGLE <- 90   # 360 - 340 = 20 degrees gap
  
  p_circ <- ggtree(
    tree_grouped,
    layout = "circular",
    open.angle = OPEN_ANGLE,
    aes(color = group),
    size = 0.2
  ) %<+% annot_all
  
  collapse_nodes <- setdiff(collapse_nodes, crass_nodes)
  
  for (nd in collapse_nodes) {
    p_circ <- collapse(p_circ, node = nd)
  }
  
  # ---- NEW: collapse biggest non-Crass clade (outside Crass MRCA) ----
  collapse_noncrass <- TRUE  # add this flag near your other options if you want
  if (collapse_noncrass && !is.na(biggest_noncrass_node)) {
    message("[COLLAPSE] Biggest non-Crass clade node=", biggest_noncrass_node,
            " tips=", biggest_noncrass_size)
    p_circ <- collapse(p_circ, node = biggest_noncrass_node)
  }
  
  # ---- Show MRCA node of all Crassvirales tips on the plot ----
  if (SHOW_CRASS_MRCA_NODE) {
    
    # MRCA point
    p_circ <- p_circ +
      geom_point2(aes(subset = (node == crass_mrca_node)), size = 2, colour = "black") +
      geom_text2(
        aes(subset = (node == crass_mrca_node), label = paste0("MRCA Crass: ", node)),
        size = 2, nudge_x = 0.02
      )
    
    if (!is.na(candidate_only_mrca_node)) {
      p_circ <- p_circ +
        geom_point2(
          aes(subset = (node == candidate_only_mrca_node)),
          size = 2,
          colour = "darkgreen"
        ) +
        geom_text2(
          aes(
            subset = (node == candidate_only_mrca_node),
            label = paste0("MRCA 16 candidates: ", node)
          ),
          size = 2,
          nudge_x = 0.02,
          colour = "darkgreen"
        )
    }
    
    # Parent point (if exists)
    if (!is.na(crass_mrca_parent_node)) {
      p_circ <- p_circ +
        geom_point2(aes(subset = (node == crass_mrca_parent_node)), size = 2, colour = "orange") +
        geom_text2(
          aes(subset = (node == crass_mrca_parent_node), label = paste0("Parent: ", node)),
          size = 2, nudge_x = 0.02
        )
    }
    
    # MRCA of Crass refs + integrated-labeled tips
    if (!is.na(crass_plus_integrated_mrca_node)) {
      p_circ <- p_circ +
        geom_point2(aes(subset = (node == crass_plus_integrated_mrca_node)),
                    size = 2, colour = "purple") +
        geom_text2(
          aes(subset = (node == crass_plus_integrated_mrca_node),
              label = paste0("MRCA Crass+Int: ", node)),
          size = 2, nudge_x = 0.02
        )
    }
  }
  
  
  
  
  # ---- Label only prophage tips with size > 50 kb ----
  # ---- Optional prophage labels (>50 kb) ----
  # if (SHOW_PROPHAGE_LABELS) {
  #   
  #   tip_label_df <- p_circ$data %>%
  #     dplyr::filter(
  #       isTip,
  #       is_prophage,
  #       !is.na(bact_genomad_length_num),
  #       bact_genomad_length_num > 50000
  #     )
  #   
  #   if (nrow(tip_label_df) > 0) {
  #     p_circ <- p_circ +
  #       geom_tiplab(
  #         data = tip_label_df,
  #         aes(label = label),
  #         size = 1.3,
  #         offset = 0.02,
  #         align = FALSE,
  #         linetype = "dotted",
  #         linewidth = 0.2
  #       )
  #   }
  # }
  
  SHOW_PROPHAGE_LABELS <- TRUE  # or keep your flag
  
  if (SHOW_PROPHAGE_LABELS) {
    tip_label_df <- p_circ$data %>%
      dplyr::filter(
        isTip, is_prophage,
        !is.na(integrated_prophage_label),
        integrated_prophage_label != ""
      )
    
    if (nrow(tip_label_df) > 0) {
      p_circ <- p_circ +
        geom_tiplab(
          data = tip_label_df,
          aes(label = integrated_prophage_label),
          size = 1.6, offset = 0.02,
          align = FALSE, linetype = "dotted", linewidth = 0.2
        )
    }
  }
  
  
  
  
  # ---- Clade highlights (behind branches) ----
  # ---- Highlight clades (family-colored) ----
  if (nrow(fam_clades_df) > 0) {
    p_circ <- p_circ +
      geom_hilight(
        data = fam_clades_df,
        aes(node = node, fill = fam),
        alpha = 0.5,     # transparency
        extend = 0.4,     # how far outward the highlight goes
        show.legend = FALSE
      ) +
      scale_fill_manual(values = CRASSVIRALES_COLOR_SCHEME[crass_families])
  }
  
  
  
  # if (SHOW_TIP_LABELS) {
  #   p_circ <- p_circ +
  #     geom_tiplab(
  #       aes(label = label),
  #       size = 1.5,
  #       offset = 0.01,
  #       align = FALSE,
  #       linetype = "dotted"
  #     )
  # }
  
  
  p_circ$data <- p_circ$data %>%
    left_join(node_sizes_df, by = "node")
  
  max_x <- max(p_circ$data$x, na.rm = TRUE)
  
  status_ring_r <- max_x * 1.02
  base_radius   <- max_x * 1.08
  scale_factor  <- 0.3
  
  # ---- Bar nodes (leaf count) ----
  bar_nodes <- p_circ$data %>%
    filter(isTip | node %in% collapse_nodes) %>%
    filter(!is.na(bar_len), is.finite(bar_len)) %>%
    mutate(
      bar_group = if_else(group %in% crass_families, as.character(group), "Other")
    )
  
  # ---- Dotted reference rings for clade sizes 10, 100, 1000 ----
  # ---- Dotted reference rings for clade sizes 1, 10, 100 ----
  clade_sizes <- c(1, 10, 100)               # <- this is the key change
  log_vals    <- log10(clade_sizes)
  
  ring_positions <- tibble(
    clade_size = clade_sizes,
    ring_x     = base_radius + (log_vals + 0.1) * scale_factor
  )
  
  y_range <- range(p_circ$data$y, na.rm = TRUE)
  
  ring_segments_df <- ring_positions %>%
    mutate(
      y    = y_range[1],
      yend = y_range[2]
    )
  
  # ---- Labels for dotted clade-size rings (show only numbers) ----
  ring_label_df <- ring_positions %>%
    mutate(
      y = y_range[2],                 # place near top; change to y_range[1] for bottom
      label = as.character(clade_size),
      x = ring_x + 0.03               # small offset to the right; tweak if needed
    )
  
  
  # ---- Status ring ----
  status_ring_df <- p_circ$data %>%
    filter(isTip, !is.na(y)) %>%
    mutate(
      ring_group = case_when(
        is_crassvirales ~ "ref_crass",
        !is.na(bact_topology) & bact_topology == "Provirus" ~ "integrated_provirus",
        !is.na(bact_topology) & bact_topology != "Provirus" ~ "non_integrated_prophage",
        TRUE                                                ~ "other"
      ),
      ring_x    = status_ring_r,
      ring_xend = status_ring_r + 0.03
    )
  
  # ---- genomad_length ring ----
  bar_outer_r    <- base_radius + max(bar_nodes$bar_len) * scale_factor
  genomad_base_r <- bar_outer_r + 0.4
  genomad_div    <- 500000  # controls bar height for genomad_length
  
  # ---- Dotted reference rings for genomad_length (50kb, 100kb, 150kb) ----
  genomad_levels <- c(50000, 100000, 150000)
  
  genomad_ring_positions <- tibble(
    genomad_len = genomad_levels,
    ring_x      = genomad_base_r + (genomad_levels / genomad_div)
  )
  
  y_range <- range(p_circ$data$y, na.rm = TRUE)
  
  genomad_ring_segments_df <- genomad_ring_positions %>%
    mutate(
      y    = y_range[1],
      yend = y_range[2]
    )
  
  genomad_ring_df <- p_circ$data %>%
    filter(isTip, is_prophage, !is.na(bact_genomad_length), !is.na(y)) %>%
    mutate(
      geno_group = if_else(bact_topology == "Provirus",
                           "integrated_provirus",
                           "non_integrated_prophage"),
      geno_x    = genomad_base_r,
      geno_xend = genomad_base_r + bact_genomad_length / genomad_div
    )
  
  # ---- Labels for dotted genomad_length rings (show only numbers) ----
  genomad_label_df <- genomad_ring_positions %>%
    mutate(
      y     = y_range[2],                         # top; use y_range[1] for bottom
      label = as.character(genomad_len / 1000),   # show 50/100/150 (kb)
      x     = ring_x + 0.03                       # small offset to the right
    )
  
  
  # ---- prophage_ratio ring ----
  prophage_base_r <- genomad_base_r + 0.4
  prophage_scale  <- 0.3
  
  # ---- Dotted reference rings for prophage_ratio (0.25, 0.5, 0.75) ----
  prophage_levels <- c(0.25, 0.5, 0.75, 1)
  
  prophage_ring_positions <- tibble(
    prophage_ratio = prophage_levels,
    ring_x         = prophage_base_r + prophage_levels * prophage_scale
  )
  
  prophage_ring_segments_df <- prophage_ring_positions %>%
    mutate(
      y    = y_range[1],
      yend = y_range[2]
    )
  
  
  prophage_ring_df <- p_circ$data %>%
    filter(isTip, is_prophage, !is.na(bact_prophage_ratio), !is.na(y)) %>%
    mutate(
      geno_group = if_else(bact_topology == "Provirus",
                           "integrated_provirus",
                           "non_integrated_prophage"),
      proph_x    = prophage_base_r,
      proph_xend = prophage_base_r + bact_prophage_ratio * prophage_scale
    )
  
  # ---- Labels for dotted prophage_ratio rings (show only numbers) ----
  prophage_label_df <- prophage_ring_positions %>%
    mutate(
      y     = y_range[2],                         # top; use y_range[1] for bottom
      label = format(prophage_ratio, trim = TRUE, scientific = FALSE),
      x     = ring_x + 0.03                       # small offset to the right
    )
  
  
  # ---- Phylum ring ----
  phylum_base_r <- prophage_base_r + 0.4
  
  phylum_ring_df <- p_circ$data %>%
    filter(isTip, is_prophage, !is.na(bact_phylum2), !is.na(y)) %>%
    mutate(
      phyl_x    = phylum_base_r,
      phyl_xend = phylum_base_r + 0.3
    )
  
  # ---- Class ring ----
  class_base_r <- phylum_base_r + 0.4
  
  class_ring_df <- p_circ$data %>%
    filter(isTip, is_prophage, !is.na(bact_class2), !is.na(y)) %>%
    mutate(
      class_x    = class_base_r,
      class_xend = class_base_r + 0.3
    )
  
  # ---- Color scale for viral/status/leaf ----
  COLOR_SCALE_ALL <- c(
    CRASSVIRALES_COLOR_SCHEME,
    "Other"                   = "grey60",
    "ref_crass"               = "black",
    "integrated_provirus"     = "red",
    "non_integrated_prophage" = "blue",
    # "non_integrated_prophage" = "#1F77b499",
    "other"                   = "grey60"
  )
  
  # --- Crassvirales families legend grob (TOP-LEFT) ---
  CRASS_LEGEND_COLORS <- CRASSVIRALES_COLOR_SCHEME[
    names(CRASSVIRALES_COLOR_SCHEME) != "Outgroup"
  ]
  
  p_crass_leg <- ggplot(
    tibble(cat = names(CRASS_LEGEND_COLORS), x = 1, y = seq_along(CRASS_LEGEND_COLORS)),
    aes(x, y, colour = cat)
  ) +
    geom_point(size = 3) +
    scale_color_manual(values = CRASS_LEGEND_COLORS) +
    guides(colour = guide_legend(title = "Crassvirales family")) +
    theme_void() +
    theme(legend.position = "left")
  
  crass_legend_grob <- cowplot::get_legend(p_crass_leg)
  
  # --- Prophage integration status legend grob (TOP-RIGHT) ---
  STATUS_COLORS <- c(
    "Integrated prophage"     = "red",
    "Non-integrated prophage" = "blue"
    # "Non-integrated prophage" =  "#1F77b499"
  )
  
  p_status_leg <- ggplot(
    tibble(cat = names(STATUS_COLORS), x = 1, y = seq_along(STATUS_COLORS)),
    aes(x, y, colour = cat)
  ) +
    geom_point(size = 3) +
    scale_color_manual(values = STATUS_COLORS) +
    guides(colour = guide_legend(title = "Prophage status")) +
    theme_void() +
    theme(legend.position = "right")
  
  status_legend_grob <- cowplot::get_legend(p_status_leg)
  
  
  
  # --- Phylum legend grob ---
  p_phylum_leg <- ggplot(
    tibble(cat = names(PHYLUM_COLORS), x = 1, y = seq_along(PHYLUM_COLORS)),
    aes(x, y, colour = cat)
  ) +
    geom_point(size = 3) +
    scale_color_manual(values = PHYLUM_COLORS) +
    guides(colour = guide_legend(title = "Bacterial phylum")) +
    theme_void() +
    theme(legend.position = "left")
  
  phylum_legend_grob <- cowplot::get_legend(p_phylum_leg)
  
  # --- Class legend grob ---
  p_class_leg <- ggplot(
    tibble(cat = names(CLASS_COLORS), x = 1, y = seq_along(CLASS_COLORS)),
    aes(x, y, colour = cat)
  ) +
    geom_point(size = 3) +
    scale_color_manual(values = CLASS_COLORS) +
    guides(colour = guide_legend(title = "Bacterial class")) +
    theme_void() +
    theme(legend.position = "right")
  
  class_legend_grob <- cowplot::get_legend(p_class_leg)
  
  
  p_circ <- p_circ +
    # 1) Viral families, status, leaf-count bars
    scale_color_manual(values = COLOR_SCALE_ALL, na.value = "black") +
    theme_tree() +
    theme(legend.position = "none") +   # <- IMPORTANT: turn off default legends (we'll place Phylum/Class manually)
    
    # Dotted guide rings for clade sizes 1, 10, 100
    geom_segment(
      data = ring_segments_df,
      aes(x = ring_x, xend = ring_x, y = y, yend = yend),
      inherit.aes = FALSE,
      colour = "grey50",
      linetype = "dotted",
      linewidth = 0.1,
      lineend = "round"
    ) +
    
    geom_text(
      data = ring_label_df,
      aes(x = x, y = y, label = label),
      inherit.aes = FALSE,
      size = 1,
      vjust = -0.2,
      angle = 45,
      colour = "grey30"
    ) +
    
    # Status ring
    geom_segment(
      data = status_ring_df,
      aes(x = ring_x, xend = ring_xend, y = y, yend = y, colour = ring_group),
      linewidth = 0.5,
      inherit.aes = FALSE
    ) +
    
    # Leaf count bars
    geom_segment(
      data = bar_nodes,
      aes(
        x    = base_radius,
        xend = base_radius + bar_len * scale_factor,
        y    = y,
        yend = y,
        colour = bar_group
      ),
      linewidth = 0.5,
      inherit.aes = FALSE
    ) +
    
    # Genomad dotted guide rings (50kb, 100kb, 150kb)
    geom_segment(
      data = genomad_ring_segments_df,
      aes(x = ring_x, xend = ring_x, y = y, yend = yend),
      inherit.aes = FALSE,
      colour = "grey50",
      linetype = "dotted",
      linewidth = 0.1,
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
    
    # 2) New colour scale for genomad/prophage rings
    ggnewscale::new_scale_color() +
    
    # Genomad length bars
    geom_segment(
      data = genomad_ring_df,
      aes(x = geno_x, xend = geno_xend, y = y, yend = y, colour = geno_group),
      linewidth = 0.5,
      inherit.aes = FALSE
    ) +
    
    # Prophage_ratio dotted guide rings (0.25, 0.5, 0.75)
    geom_segment(
      data = prophage_ring_segments_df,
      aes(x = ring_x, xend = ring_x, y = y, yend = yend),
      inherit.aes = FALSE,
      colour = "grey50",
      linetype = "dotted",
      linewidth = 0.1,
      lineend = "round"
    ) +
    
    geom_text(
      data = prophage_label_df,
      aes(x = x, y = y, label = label),
      inherit.aes = FALSE,
      size = 0.65,
      vjust = -0.2,
      angle = 45,
      colour = "grey30"
    ) +
    
    
    # Prophage_ratio bars
    geom_segment(
      data = prophage_ring_df,
      aes(x = proph_x, xend = proph_xend, y = y, yend = y, colour = geno_group),
      linewidth = 0.5,
      inherit.aes = FALSE
    ) +
    
    scale_color_manual(
      values = c(
        "integrated_provirus"     = "red",
        "non_integrated_prophage" = "blue"
        # "non_integrated_prophage" = "#1F77b499"
      )
    ) +
    
    # 3) New colour scale for Phylum ring (NO legend here; legend added manually below)
    ggnewscale::new_scale_color() +
    
    geom_segment(
      data = phylum_ring_df,
      aes(x = phyl_x, xend = phyl_xend, y = y, yend = y, colour = bact_phylum2),
      linewidth = 0.5,
      inherit.aes = FALSE,
      show.legend = FALSE
    ) +
    
    scale_color_manual(values = PHYLUM_COLORS) +
    
    # 4) New colour scale for Class ring (NO legend here; legend added manually below)
    ggnewscale::new_scale_color() +
    
    geom_segment(
      data = class_ring_df,
      aes(x = class_x, xend = class_xend, y = y, yend = y, colour = bact_class2),
      linewidth = 0.5,
      inherit.aes = FALSE,
      show.legend = FALSE
    ) +
    
    scale_color_manual(values = CLASS_COLORS)
  
  # =========================
  # ---- Manual Phylum + Class legends in bottom corners ----
  # =========================
  # (Requires: library(cowplot))
  
  # Phylum legend grob
  p_phylum_leg <- ggplot(
    tibble(cat = names(PHYLUM_COLORS), x = 1, y = seq_along(PHYLUM_COLORS)),
    aes(x, y, colour = cat)
  ) +
    geom_point(size = 3) +
    scale_color_manual(values = PHYLUM_COLORS) +
    guides(colour = guide_legend(title = "Phylum")) +
    theme_void() +
    theme(legend.position = "left")
  
  phylum_legend_grob <- cowplot::get_legend(p_phylum_leg)
  
  # Class legend grob
  p_class_leg <- ggplot(
    tibble(cat = names(CLASS_COLORS), x = 1, y = seq_along(CLASS_COLORS)),
    aes(x, y, colour = cat)
  ) +
    geom_point(size = 3) +
    scale_color_manual(values = CLASS_COLORS) +
    guides(colour = guide_legend(title = "Class", ncol = 3, byrow = TRUE)) +  # <-- HERE
    theme_void() +
    theme(legend.position = "right")
  
  class_legend_grob <- cowplot::get_legend(p_class_leg)
  
  if (SHOW_CRASS_LABELS) {
    
    tip_label_df <- p_circ$data %>%
      filter(
        isTip,
        is_crassvirales
      )
    
    p_circ <- p_circ +
      geom_tiplab(
        data = tip_label_df,
        aes(label = label),
        size = 1.5,
        offset = 0.01,
        align = FALSE,
        linetype = "dotted"
      )
  }
  
  # after you compute all ring bases (class_base_r etc.)
  # outer_r <- class_base_r + 0.35
  # 
  # p_circ <- p_circ +
  #   coord_equal(
  #     xlim = c(-outer_r, outer_r),
  #     ylim = c(-outer_r, outer_r),
  #     clip = "off"
  #   )
  
  
  # Compose final plot with legends in corners
  p_final <- cowplot::ggdraw() +
    # main plot
    cowplot::draw_plot(
      p_circ,
      x = 0, y = 0, width = 1, height = 1
    ) +
    
    # top-left: Crassvirales families
    cowplot::draw_grob(
      crass_legend_grob,
      x = 0, y = 0.65, width = 0.30, height = 0.18
    ) +
    
    # top-right: Prophage status
    cowplot::draw_grob(
      status_legend_grob,
      x = 0.62, y = 0.75, width = 0.30, height = 0.12
    ) +
    
    # bottom-left: Phylum
    cowplot::draw_grob(
      phylum_legend_grob,
      x = 0, y = 0.2, width = 0.30, height = 0.08
    ) +
    
    # bottom-right: Class
    cowplot::draw_grob(
      class_legend_grob,
      x = 0.4, y = 0.07, width = 0.36, height = 0.08
    )
  
  
  
  # =========================
  # ---- Save ----
  # =========================
  ggsave(out_png_circ, p_final,
         width = 18, height = 22, units = "cm", dpi = 1200)
  
  cat("Saved:", out_png_circ, "\n")
  }
}