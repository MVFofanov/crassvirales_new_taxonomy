#!/usr/bin/env Rscript

# ======================= deps =======================
pkgs_cran <- c("ape","ggplot2","gggenes","dplyr","readr","tibble","tidyr","stringr")
for (pkg in pkgs_cran) if (!requireNamespace(pkg, quietly=TRUE)) install.packages(pkg, repos="https://cloud.r-project.org")
if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager", repos="https://cloud.r-project.org")
for (pkg in c("ggtree","treeio")) if (!requireNamespace(pkg, quietly=TRUE)) BiocManager::install(pkg, ask=FALSE, update=FALSE)
if (!requireNamespace("ggnewscale", quietly=TRUE)) install.packages("ggnewscale", repos="https://cloud.r-project.org")

library(ape)
library(ggplot2)
library(ggtree)
library(gggenes)
library(treeio)
library(ggnewscale)
library(dplyr)
library(readr)
library(tibble)
library(tidyr)
library(stringr)

# ======================= palettes =======================
DEFAULT_COLOR <- "#bfbfbf"
  
crass_family_colors <- c(
  "Intestiviridae"="#EE3B3B","Crevaviridae"="#EE9A00","Suoliviridae"="#4169E1",
  "Steigviridae"="#00CED1","Epsilon"="#CD2990","Zeta"="#006400"
)
non_cv_colors <- c("Unknown"="#BDBDBD","outgroup"="#808080","NA"="#BDBDBD")
family_palette <- c(crass_family_colors, non_cv_colors)

origin_palette <- c("MAG"="#1f77b4", "Isolate"="#d62728")

phylum_color_map <- c(
  "Bacillota"="#1f77b4","Firmicutes"="#1f77b4","Bdellovibrionota"="#bfef45",
  "Bacteroidota"="#2ca02c","Bacteroidetes"="#2ca02c","Proteobacteria"="#9A6324",
  "Pseudomonadota"="#9A6324","Mycoplasmatota"="#911eb4","Chloroflexota"="#ffe119",
  "Ignavibacteriota"="#fabed4","Candidatus Melainabacteria"="#42d4f4",
  "Candidatus Borrarchaeota"="#000000","Candidatus Pacearchaeota"="#000000",
  "Methanobacteriota"="#000000","Unknown"=DEFAULT_COLOR
)

class_color_map <- c(
  "Bacilli"="#e6194B","Clostridia"="#f58231","Erysipelotrichia"="#808000",
  "Chitinophagia"="#f032e6","Cytophagia"="#000075","Ignavibacteria"="#fabed4",
  "Bdellovibrionia"="#bfef45","Mollicutes"="#911eb4","Vampirovibriophyceae"="#42d4f4",
  "Candidatus Borrarchaeota"="#000000","Candidatus Pacearchaeota"="#000000",
  "Halobacteria"="#000000","Methanobacteria"="#000000",
  "Unknown"="#999999"
)

NOT_CV_TAGS <- c("outgroup","NA")

# ======================= paths =======================
proj_dir  <- "C:/crassvirales/crassvirales_new_taxonomy/crassvirales_prophages/blast_prophages_vs_ncbi_and_gtdb"
crassus_results_dir <- "C:/crassvirales/CrassUS_old/CrassUS/results/prophage_analysis"

tree_file <- file.path(crassus_results_dir, "5_phylogenies/3_iToL/Terl_iToL_renamed.nwk")
ann_file  <- file.path(proj_dir, "Crassphage_prophage_analysis_annotation.tsv")
itol_file <- file.path(crassus_results_dir, "5_phylogenies/3_iToL/TerL_iToL_simplified.txt")
itol_simplified_file <- file.path(dirname(itol_file), "TerL_iToL_simplified.txt")

genes_file <- file.path(crassus_results_dir, "4_ORF/3_functional_annot_table_renamed.tsv")

out_png <- file.path(proj_dir, "Crassphage_prophage_analysis_with_tree_plot.png")
out_svg <- file.path(proj_dir, "Crassphage_prophage_analysis_with_tree_plot.svg")
out_pdf <- file.path(proj_dir, "Crassphage_prophage_analysis_with_tree_plot.pdf")

cat("Tree file: ", normalizePath(tree_file, winslash="/"), "\n")
cat("Prophage annotation: ", normalizePath(ann_file, winslash="/"), "\n")
cat("iTOL simplified: ", normalizePath(itol_file, winslash="/"), "\n")

# ======================= 1) load tree =======================
tr <- read.tree(tree_file)

# ======================= 2) load + normalize annotations =======================
ann_raw <- suppressMessages(readr::read_tsv(ann_file, show_col_types = FALSE))

# pick best join key (protein_id vs genome_id)
possible_keys <- c("genome_id","protein_id")
present_keys  <- possible_keys[possible_keys %in% names(ann_raw)]
if (!length(present_keys)) stop("Neither genome_id nor protein_id present in annotation table.")
overlap  <- vapply(present_keys, function(k) sum(ann_raw[[k]] %in% tr$tip.label), integer(1))
best_key <- present_keys[which.max(overlap)]
cat(sprintf("Joining annotations on '%s' (matches %d/%d tips)\n",
            best_key, max(overlap), length(tr$tip.label)))

ann <- ann_raw %>%
  mutate(
    label = as.character(.data[[best_key]]),
    
    # normalize origin
    origin = ifelse(genome_origin %in% c("MAG","Isolate"),
                    as.character(genome_origin),
                    NA_character_),
    
    # normalize family/phyla/class and map unknowns to "Unknown"
    prophage_family = str_trim(crassvirales_family),
    prophage_family = ifelse(prophage_family %in% names(family_palette), prophage_family, "Unknown"),
    
    bacterial_phylum = str_trim(bacterial_phylum),
    bacterial_phylum = ifelse(bacterial_phylum %in% names(phylum_color_map), bacterial_phylum, "Unknown"),
    
    bacterial_class = str_trim(bacterial_class),
    bacterial_class = ifelse(bacterial_class %in% names(class_color_map), bacterial_class, "Unknown"),
    
    bacterial_order_family  = as.character(bacterial_order_and_family),
    bacterial_tax_lineage   = as.character(bacterial_taxonomy_lineage),
    bacterial_id            = as.character(bacterial_id),
    prophage_length         = suppressWarnings(as.numeric(prophage_length)),
    bacterial_contig_length = suppressWarnings(as.numeric(bacterial_contig_length)),
    prophage_start         = suppressWarnings(as.numeric(prophage_start)),
    prophage_end           = suppressWarnings(as.numeric(prophage_end)),
    genome_id               = as.character(genome_id),
    protein_id              = as.character(protein_id)
  ) %>%
  select(label, prophage_family, origin, bacterial_phylum, bacterial_class,
         bacterial_order_family, bacterial_tax_lineage, bacterial_id,
         prophage_length, bacterial_contig_length, prophage_start, prophage_end,
         genome_id, protein_id)

# ======================= 3) build tip metadata =======================
itol <- suppressMessages(readr::read_tsv(itol_simplified_file, show_col_types = FALSE)) %>%
  transmute(label = as.character(genome_id),
            itol_family = as.character(crassvirales_family))

tips_tbl <- tibble(label = tr$tip.label) %>%
  left_join(itol, by="label") %>%
  left_join(ann,  by="label") %>%
  mutate(
    is_prophage     = !is.na(prophage_family),
    family          = dplyr::coalesce(prophage_family, itol_family),
    family          = ifelse(family %in% names(family_palette), family, "Unknown"),
    in_crassvirales = !is.na(family) & !(family %in% NOT_CV_TAGS) & family != "Unknown"
  )
cat("Matched prophage rows: ", sum(tips_tbl$is_prophage, na.rm=TRUE), " of ", nrow(tips_tbl), "\n", sep="")

# ======================= 4) prune to MRCA of Crassvirales =======================
cv_tips <- tips_tbl %>% filter(in_crassvirales) %>% pull(label) %>% unique()
if (length(cv_tips) < 2) stop("Not enough Crassvirales tips to define an MRCA (need >= 2).")
mrca_node <- ape::getMRCA(tr, cv_tips)
if (is.null(mrca_node) || is.na(mrca_node)) stop("MRCA not found for the provided Crassvirales tips.")

cv_clade      <- ape::extract.clade(tr, mrca_node)
cv_tip_labels <- cv_clade$tip.label
tr_pruned     <- ape::keep.tip(tr, cv_tip_labels)

tips_tbl_pruned <- tips_tbl %>%
  filter(label %in% tr_pruned$tip.label) %>%
  mutate(family = factor(family, levels = names(family_palette)))

# ======================= 5) choose non-prophage nodes to collapse =======================
internal_nodes <- (Ntip(tr_pruned) + 1):(Ntip(tr_pruned) + tr_pruned$Nnode)
desc_tip_count <- desc_prophage_count <- integer(length(internal_nodes))
for (i in seq_along(internal_nodes)) {
  nd <- internal_nodes[i]
  tips_under <- ape::extract.clade(tr_pruned, nd)$tip.label
  desc_tip_count[i]      <- length(tips_under)
  desc_prophage_count[i] <- sum(tips_tbl_pruned$is_prophage[match(tips_under, tips_tbl_pruned$label)], na.rm=TRUE)
}
collapse_nodes <- internal_nodes[desc_tip_count > 1 & desc_prophage_count == 0]
collapsed_info <- tibble(node = collapse_nodes,
                         collapsed_n = desc_tip_count[match(collapse_nodes, internal_nodes)])
message(sprintf("Will collapse %d clades (pure non-prophage); total collapsed tips across them = %d",
                length(collapse_nodes), sum(collapsed_info$collapsed_n)))

# ======================= 5b) fill families for non-prophage (incl. collapsed clades) =======================
# Start from current family calls
tips_tbl_filled <- tips_tbl_pruned %>%
  mutate(family_all = family)

# For each collapsed pure non-prophage clade, pick the majority (or sole) family among its descendant tips
for (nd in collapse_nodes) {
  desc <- ape::extract.clade(tr_pruned, nd)$tip.label
  fams <- tips_tbl_pruned$family[match(desc, tips_tbl_pruned$label)]
  fams_clean <- fams[!is.na(fams) & fams != "Unknown" & !(fams %in% NOT_CV_TAGS)]
  if (length(fams_clean) > 0) {
    fam_mode <- names(sort(table(fams_clean), decreasing = TRUE))[1]
    # assign to all descendant tips for the "family_all" column
    idx <- match(desc, tips_tbl_filled$label)
    tips_tbl_filled$family_all[idx] <- factor(fam_mode, levels = levels(tips_tbl_pruned$family))
  }
}

# Family call per collapsed clade (mode among descendants)
clade_fam_tbl <- lapply(collapse_nodes, function(nd) {
  desc <- ape::extract.clade(tr_pruned, nd)$tip.label
  fams <- tips_tbl_pruned$family[match(desc, tips_tbl_pruned$label)]
  fams_clean <- fams[!is.na(fams) & fams != "Unknown" & !(fams %in% NOT_CV_TAGS)]
  fam_mode <- if (length(fams_clean) > 0) {
    as.character(names(sort(table(fams_clean), decreasing = TRUE))[1])
  } else {
    "Unknown"
  }
  tibble(node = nd, fam = fam_mode)
}) %>% bind_rows()


# ======================= 6) geometry for labels & panels =======================
# preview (collapsed) to measure in data units
p0 <- ggtree(tr_pruned, layout="rectangular") %<+% tips_tbl_pruned
for (nd in collapse_nodes) p0 <- collapse(p0, node=nd)
x_span <- diff(range(p0$data$x, na.rm=TRUE))

lab_offset  <- 0.03 * x_span   # distance between tips and labels
panel_shift <- 1.7 * x_span   # push block of panels to the right of labels, 0.060 * x_span
panel_gap   <- 0.15 * x_span   # uniform gap between panels, 0.012 * x_span
panel_w     <- 0.045 * x_span   # SAME width for all panels (change this one number), 0.045 * x_span 

## ---- panel layout (now with a 'len' panel) ----
# widths: keep existing four panels at panel_w; give length panel a bit more room
# len_width <- 0.06 * x_span  # tweak to taste

# ---- panel layout: 4 heatmap panels only ----
## ---- panel layout (left → right) ----

# 4 heatmap panels
panel_names  <- c("family","origin","phylum","class")
panel_widths <- setNames(rep(panel_w, length(panel_names)), panel_names)
cum_left       <- c(0, cumsum(head(panel_widths, -1) + panel_gap))
panel_offsets  <- lab_offset + panel_shift + cum_left
names(panel_offsets) <- panel_names
panels_total_width <- sum(panel_widths) + panel_gap * (length(panel_widths) - 1)

# contig panel (gray bar + red prophage segment)
contig_gap    <- 1.5 * x_span
contig_width  <- 0.6 * x_span
contig_offset <- lab_offset + panel_shift + panels_total_width + contig_gap

# gene map panel (per-ORF rectangles), to the RIGHT of contig
gene_gap    <- 0.05 * x_span
gene_width  <- 2  * x_span
gene_offset <- contig_offset + contig_width + gene_gap

# length bars panel, to the RIGHT of gene map
bar_gap   <- 0.1 * x_span
bar_width <- 0.6  * x_span
len_offset <- gene_offset + gene_width + bar_gap
len_width  <- bar_width

# numeric labels to the RIGHT of length bars
value_gap   <- 0.02 * x_span
value_width <- 0.06 * x_span
value_x     <- len_offset + len_width + value_gap

# total width for xlim_tree
total_width_all <- (panels_total_width +
                      contig_gap + contig_width +
                      gene_gap   + gene_width +
                      bar_gap    + bar_width +
                      value_gap  + value_width)
extra_right <- 0.10 * x_span

# ======================= 7) base tree =======================
# --- 7) base tree ---
p <- ggtree(tr_pruned, layout="rectangular") %<+% tips_tbl_pruned +
  theme_tree() +                                # <— was theme_tree2(); removes the long x-axis
  geom_tippoint(aes(color = family), size = 1.2, alpha = 0.9, show.legend = FALSE) +
  ggtitle("TerL — Crassvirales clade (pruned), non-prophage clades collapsed")

for (nd in collapse_nodes) p <- collapse(p, node = nd)

pd <- p$data

y_min_tip <- min(pd$y[pd$isTip], na.rm = TRUE)
row_gap   <- 0.9                      # vertical spacing between scale rows
y_tree_sc <- y_min_tip - 0.6          # tree scale
y_contig_sc <- y_tree_sc - row_gap    # under "Prophage coordinates"
y_gene_sc   <- y_contig_sc - row_gap  # under "Prophage genomic map"
y_len_sc    <- y_gene_sc   - row_gap  # under "Prophage length"
y_bottom_all <- y_len_sc - 0.6        # bottom limit to include all scales

# --- add a short tree scale bar (only over the tree area) ---
## --- place treescale BELOW the tree --- 
# (right after: pd <- p$data)
tree_xmin <- min(pd$x, na.rm = TRUE)
tree_xmax <- max(pd$x, na.rm = TRUE)
tree_span <- tree_xmax - tree_xmin
y_min_tip <- min(pd$y[pd$isTip], na.rm = TRUE)

BOTTOM_PAD_TIPS <- 4          # <- how far below the lowest tip to place the scale (in "rows")
y_bottom        <- y_min_tip - BOTTOM_PAD_TIPS

p <- p + ggtree::geom_treescale(
  x = tree_xmin + 0.02 * tree_span,
  y = y_bottom + 0.8,         # slightly above the bottom limit so the text isn’t cut
  width = 0.12 * tree_span,   # scale bar length; adjust to taste
  linesize = 0.5,
  fontsize = 3
)

lab_pos <- pd %>%
  dplyr::select(node, x, y) %>%
  dplyr::inner_join(collapsed_info, by="node")

p <- p +
  geom_label(data = lab_pos,
             aes(x = x, y = y, label = collapsed_n),
             size = 3, label.size = 0.2,
             label.padding = unit(0.12, "lines"),
             fill = "white") +
  geom_tippoint(
    data = subset(pd, isTip & !is.na(is_prophage) & is_prophage),
    shape = 21, stroke = 0.2, size = 2, fill = NA, color = "black",
    show.legend = FALSE
  ) +
  geom_tiplab(size = 2, align = TRUE, offset = lab_offset, linesize = 0) +
  scale_color_manual(values = family_palette, na.value = "#BDBDBD", guide = "none") +
  # guides(color = guide_legend(title = "Crassvirales family (tips)",
  #                             override.aes = list(size = 3, alpha = 1), order = 1)) +
  xlim_tree(max(p$data$x, na.rm=TRUE) + lab_offset + panel_shift + total_width_all + extra_right)

# -------- helper to choose a nice kb step --------
nice_kb_step <- function(max_kb) {
  steps <- c(0.1,0.2,0.5,1,2,5,10,20,50,100,200,500,1000)
  steps[max(which(steps <= max_kb/3))] %||% steps[1]
}

`%||%` <- function(a,b) if (!is.finite(a) || length(a)==0) b else a

# Build a robust pretty axis in data units
pretty_axis <- function(max_units, n = 6) {
  if (!is.finite(max_units) || max_units <= 0) max_units <- 1
  br <- pretty(c(0, max_units), n = n)
  br <- br[br >= 0]                 # keep non-negative
  ax_max <- suppressWarnings(max(br, na.rm = TRUE))
  if (!is.finite(ax_max) || ax_max <= 0) {  # fallback if pretty() misbehaves
    ax_max <- max(max_units, 1)
    br <- seq(0, ax_max, length.out = n)
  }
  list(breaks = unique(br), max = ax_max)
}


# ===== axes (ticks) under panels, in Mb/kb (normal, no left nudge) =====

nice_step <- function(max_units) {
  steps <- c(0.1,0.2,0.5,1,2,5,10,20,50,100,200,500,1000)
  steps[max(which(steps <= max_units/3))] %||% steps[1]
}

panel_ticks_normal <- function(offset, width, max_units, step = NULL) {
  if (!is.finite(max_units) || max_units <= 0) max_units <- 1
  if (is.null(step)) step <- nice_step(max_units)
  
  br <- seq(0, max_units, by = step)
  # ensure the rightmost tick is exactly at max_units
  if (abs(tail(br, 1) - max_units) > 1e-9) br <- c(br, max_units)
  
  fmt <- if (step >= 1) "%.0f" else "%.1f"
  x   <- offset + (br / max_units) * width
  
  tibble(units = br, x = x, lab = sprintf(fmt, br))
}

# ======================= 8) side-table (use family_all for everyone) =======================
side_df <- tips_tbl_filled %>%
  transmute(
    label,
    `Crassvirales family` = as.character(family_all),                         # <-- show for all tips
    `Genome origin`       = ifelse(is_prophage, origin,                         NA_character_),
    `Bacterial phylum`    = ifelse(is_prophage, as.character(bacterial_phylum), NA_character_),
    `Bacterial class`     = ifelse(is_prophage, as.character(bacterial_class),  NA_character_)
  ) %>%
  right_join(tibble(label = tr_pruned$tip.label), by = "label") %>%
  arrange(match(label, tr_pruned$tip.label)) %>%
  tibble::column_to_rownames("label")

family_df <- side_df[, "Crassvirales family", drop = FALSE]
origin_df <- side_df[, "Genome origin",       drop = FALSE]
phylum_df <- side_df[, "Bacterial phylum",    drop = FALSE]
class_df  <- side_df[, "Bacterial class",     drop = FALSE]

family_df[[1]] <- factor(family_df[[1]], levels = names(family_palette))
origin_df[[1]] <- factor(origin_df[[1]], levels = names(origin_palette))
phylum_df[[1]] <- factor(phylum_df[[1]], levels = names(phylum_color_map))
class_df[[1]]  <- factor(class_df[[1]],  levels = names(class_color_map))

cat("Origin values in ann:\n");      print(table(ann$origin, useNA="ifany"))
cat("Origin values in side_df:\n");  print(table(side_df[["Genome origin"]], useNA="ifany"))

# ======================= bar data for prophage length =======================
# Get tip y-positions from the plotted tree (after collapse)
pd_tips <- subset(p$data, isTip)

# ============== contig panel data (one row per tip) ==============
# Join the per-tip y back to contig fields
contig_df <- pd_tips %>%
  dplyr::select(label, y) %>%
  dplyr::left_join(
    tips_tbl_pruned %>%
      dplyr::select(label, is_prophage, bacterial_contig_length,
                    prophage_start, prophage_end),
    by = "label"
  ) %>%
  dplyr::mutate(
    contig_len = suppressWarnings(as.numeric(bacterial_contig_length)),
    start_raw  = suppressWarnings(as.numeric(prophage_start)),
    end_raw    = suppressWarnings(as.numeric(prophage_end))
  )

# We'll scale gray bar length by the maximum contig among prophages so rows are comparable
contig_max <- suppressWarnings(max(contig_df$contig_len[contig_df$is_prophage %in% TRUE], na.rm = TRUE))
if (!is.finite(contig_max) || contig_max <= 0) contig_max <- 1

# --- Approach A: pretty axis for contig panel (Mb) ---
contig_max_mb <- contig_max / 1e6
br_c <- pretty(c(0, contig_max_mb))             # tick positions in Mb (nice numbers)
axis_max_mb <- max(br_c)                        # rounded, nice max (>= data max)
axis_max_bp <- axis_max_mb * 1e6                # same max in bp for bar scaling

row_half <- 0.45  # visual thickness of per-tip bars

# Background strip for every row across full contig panel width (light gray)
contig_row_bg_df <- pd_tips %>%
  dplyr::transmute(
    x0   = contig_offset,
    x1   = contig_offset + contig_width,
    ymin = y - 0.5,
    ymax = y + 0.5
  )

# Gray contig bars (length ∝ contig_len)
contig_bar_df <- contig_df %>%
  dplyr::filter(is_prophage %in% TRUE, is.finite(contig_len), contig_len > 0) %>%
  dplyr::mutate(
    bar_len = (contig_len / axis_max_bp) * contig_width,
    x0      = contig_offset,
    x1      = contig_offset + bar_len,
    ymin    = y - row_half,
    ymax    = y + row_half
  )

# Red prophage segment within each gray bar (clamped to [0, contig_len])
proph_rect_df <- contig_bar_df %>%
  dplyr::mutate(
    start_c = pmax(0, pmin(start_raw, contig_len)),
    end_c   = pmax(start_c, pmin(end_raw, contig_len)),
    # map [0, contig_len] to [0, bar_len] then shift by contig_offset
    xr0     = contig_offset + (start_c / contig_len) * bar_len,
    xr1     = contig_offset + (end_c   / contig_len) * bar_len
  ) %>%
  dplyr::filter(is.finite(xr0), is.finite(xr1), xr1 > xr0)

bars_df <- pd_tips %>%
  dplyr::select(label, y) %>%
  dplyr::left_join(tips_tbl_pruned %>% dplyr::select(label, prophage_length), by = "label") %>%
  dplyr::mutate(len_kb = prophage_length / 1000) %>%
  dplyr::filter(!is.na(len_kb))

# scale lengths into the allocated width of the 'len' panel
# scale lengths into the allocated width of the bar region
len_max <- max(bars_df$len_kb, na.rm = TRUE)
if (!is.finite(len_max) || len_max <= 0) len_max <- 1

# len_offset / len_width were defined in the layout block above
bars_df <- bars_df %>%
  dplyr::mutate(
    x0   = len_offset,
    x1   = len_offset + (len_kb / len_max) * len_width,
    ymin = y - 0.45,
    ymax = y + 0.45
  )

len_bg_df <- pd_tips %>%
  dplyr::transmute(
    x0 = len_offset, x1 = len_offset + len_width,
    ymin = y - 0.5,   ymax = y + 0.5
  )


# a slim bar height that fits each tip row (tweak 0.45–0.48 if you want thicker bars)
row_half <- 0.45

bars_df <- bars_df %>%
  dplyr::mutate(
    x0   = len_offset,
    x1   = len_offset + (len_kb / len_max) * len_width,
    ymin = y - row_half,
    ymax = y + row_half
  )

# light background strips for the whole 'len' panel (so NA rows show as light gray)
len_bg_df <- pd_tips %>%
  dplyr::transmute(
    x0 = len_offset, x1 = len_offset + len_width,
    ymin = y - 0.5,   ymax = y + 0.5
  )

# labels placed just beyond each bar; clamp to reserved value column
text_pad <- 0.01 * x_span  # small extra pad after each bar
labels_df <- bars_df %>%
  mutate(
    # propose to put text right after the bar end:
    x_text_raw = x1 + text_pad,
    # but don't let it spill beyond the reserved value column:
    x_text = pmin(x_text_raw, value_x + value_width * 0.95),
    # format length (kb) — pick 0 or 1 decimal as you prefer
    lab = sprintf("%.1f", len_kb)
  )

# ======================= gene map data (per-prophage ORFs) =======================
# The file must have: genome, start, end, strand, yutin
genes_raw <- suppressMessages(readr::read_tsv(genes_file, show_col_types = FALSE))

# normalize / keep only what we need
# ======================= gene map data (per-prophage ORFs) =======================
genes_raw <- suppressMessages(readr::read_tsv(genes_file, show_col_types = FALSE))

genes_norm <- genes_raw %>%
  transmute(
    genome = as.character(genome),
    start  = suppressWarnings(as.numeric(start)),
    end    = suppressWarnings(as.numeric(end)),
    strand = ifelse(strand %in% c("+","-"), strand, NA_character_),
    yutin  = as.character(yutin)
  )

# Join to tips to get row y and prophage length
gene_df <- genes_norm %>%
  inner_join(
    tips_tbl_pruned %>% select(label, genome_id, is_prophage, prophage_length),
    by = c("genome" = "genome_id")
  ) %>%
  inner_join(pd_tips %>% select(label, y), by = "label") %>%
  mutate(pl = suppressWarnings(as.numeric(prophage_length))) %>%
  filter(is_prophage %in% TRUE, is.finite(pl), pl > 0,
         is.finite(start), is.finite(end))

# Global max prophage length (in bp)
pl_max <- suppressWarnings(
  max(tips_tbl_pruned$prophage_length[tips_tbl_pruned$is_prophage %in% TRUE], na.rm = TRUE)
)
if (!is.finite(pl_max) || pl_max <= 0) pl_max <- 1

# --- Approach A: pretty axis for gene map (kb) ---
pl_max_kb <- pl_max / 1e3
br_g <- pretty(c(0, pl_max_kb))          # tick positions in kb
axis_max_kb <- max(br_g)                 # rounded, nice max
pl_axis_max_bp <- axis_max_kb * 1e3      # same max in bp for mapping x0/x1

# Map gene coords to the gene panel and add colors
row_half_gene <- 0.28
gene_df <- gene_df %>%
  mutate(
    # existing fields...
    s = pmax(1, pmin(start, pl)),
    e = pmax(s, pmin(end,   pl)),
    # AFTER  (scale by pretty axis max)
    x0 = gene_offset + ((s - 1) / pl_axis_max_bp) * gene_width,
    x1 = gene_offset + ( e       / pl_axis_max_bp) * gene_width,
    ymin = y - row_half_gene,
    ymax = y + row_half_gene,
    # map to category labels for the legend
    gene_type = ifelse(!is.na(yutin) & yutin != "", "Crassphage gene", "Hypothetical protein")
  ) %>%
  filter(is.finite(x0), is.finite(x1), x1 > x0)


# Background strips for the full gene panel
gene_row_bg_df <- pd_tips %>%
  transmute(
    x0 = gene_offset, x1 = gene_offset + gene_width,
    ymin = y - 0.5,   ymax = y + 0.5
  )

# ---- build arrow-shaped polygons for genes (MUST come after gene_df has x0/x1/ymin/ymax) ----
head_frac <- 0.33
head_min  <- 0.004 * x_span
head_max  <- 0.10  * gene_width

gene_polys <- gene_df %>%
  mutate(
    id    = dplyr::row_number(),
    ymid  = (ymin + ymax) / 2,
    width = pmax(x1 - x0, .Machine$double.eps),
    head  = pmin(pmax(width * head_frac, head_min), head_max)
  ) %>%
  rowwise() %>%
  do({
    g <- .
    if (!g$strand %in% c("+","-")) {
      xs <- c(g$x0, g$x1, g$x1, g$x0, g$x0)
      ys <- c(g$ymin, g$ymin, g$ymax, g$ymax, g$ymin)
    } else if (g$strand == "+") {
      xs <- c(g$x0, g$x1 - g$head, g$x1,        g$x1 - g$head, g$x0)
      ys <- c(g$ymin, g$ymin,      g$ymid,      g$ymax,        g$ymax)
    } else { # strand == "-"
      xs <- c(g$x0 + g$head, g$x1,  g$x1,       g$x0 + g$head, g$x0)
      ys <- c(g$ymin,        g$ymin, g$ymax,    g$ymax,        g$ymid)
    }
    tibble(id = g$id, x = xs, y = ys, gene_type = g$gene_type)
  }) %>%
  ungroup()

# === Build ticks (now that br_c/br_g exist) and draw axes ===

# ticks (CONTIG, Mb)
step_mb <- if (length(br_c) >= 2) diff(br_c)[1] else axis_max_mb
fmt_mb  <- if (is.finite(step_mb) && step_mb >= 1) "%.0f" else "%.1f"
ticks_c <- tibble(
  units = br_c,
  x     = contig_offset + (br_c / axis_max_mb) * contig_width,
  lab   = sprintf(fmt_mb, br_c)
)

# ticks (GENE MAP, kb)
step_kb <- if (length(br_g) >= 2) diff(br_g)[1] else axis_max_kb
fmt_kb  <- if (is.finite(step_kb) && step_kb >= 1) "%.0f" else "%.1f"
ticks_g <- tibble(
  units = br_g,
  x     = gene_offset + (br_g / axis_max_kb) * gene_width,
  lab   = sprintf(fmt_kb, br_g)
)

# guards used by gridlines
GRID_ALPHA <- if (exists("GRID_ALPHA")) GRID_ALPHA else 0.12
ymin_rows  <- min(pd_tips$y, na.rm = TRUE) - 0.5
ymax_rows  <- max(pd_tips$y, na.rm = TRUE) + 0.5

# styling
tick_len   <- 0.25
label_pad  <- 0.45
axis_lwd   <- 0.4
tick_lwd   <- 0.35
label_size <- 3

# baseline axes
p <- p +
  geom_segment(aes(x = contig_offset, xend = contig_offset + contig_width,
                   y = y_contig_sc,  yend = y_contig_sc),
               inherit.aes = FALSE, linewidth = axis_lwd) +
  geom_segment(aes(x = gene_offset,   xend = gene_offset   + gene_width,
                   y = y_gene_sc,     yend = y_gene_sc),
               inherit.aes = FALSE, linewidth = axis_lwd)

# ticks + labels (CONTIG, Mb)
p <- p +
  geom_segment(data = ticks_c,
               aes(x = x, xend = x, y = y_contig_sc - tick_len, yend = y_contig_sc + tick_len),
               inherit.aes = FALSE, linewidth = tick_lwd) +
  geom_text(data = ticks_c,
            aes(x = x, y = y_contig_sc - (tick_len + label_pad), label = lab),
            inherit.aes = FALSE, vjust = 1, size = label_size)

# ticks + labels (GENE MAP, kb)
p <- p +
  geom_segment(data = ticks_g,
               aes(x = x, xend = x, y = y_gene_sc - tick_len, yend = y_gene_sc + tick_len),
               inherit.aes = FALSE, linewidth = tick_lwd) +
  geom_text(data = ticks_g,
            aes(x = x, y = y_gene_sc - (tick_len + label_pad), label = lab),
            inherit.aes = FALSE, vjust = 1, size = label_size)

# OPTIONAL vertical gridlines
p <- p +
  geom_segment(data = ticks_c,
               aes(x = x, xend = x, y = ymin_rows, yend = ymax_rows),
               inherit.aes = FALSE, alpha = GRID_ALPHA) +
  geom_segment(data = ticks_g,
               aes(x = x, xend = x, y = ymin_rows, yend = ymax_rows),
               inherit.aes = FALSE, alpha = GRID_ALPHA)



# --- 9) draw panels — separate legends for each heatmap & keep tree legend ---

# 9a) FAMILY heatmap + its own legend
p <- gheatmap(
  p, family_df,
  offset = as.numeric(panel_offsets["family"]),
  width  = as.numeric(panel_widths["family"]),
  colnames = FALSE, font.size = 3, hjust = 0, color = NA
)

p <- p + scale_fill_manual(
  values = family_palette,
  breaks = names(family_palette),
  drop = FALSE, na.value = "#F2F2F2",
  name = "Crassvirales families",
  guide = guide_legend(
    order = 2, ncol = 1, byrow = TRUE, title.position = "top",
    override.aes = list(alpha = 1)   # make legend swatches opaque
  )
)

# --- read FAMILY panel extents (must be before new_scale_fill) ---
gb <- ggplot_build(p)
tile_layers <- which(vapply(p$layers, function(L) inherits(L$geom, "GeomTile"), logical(1)))
stopifnot(length(tile_layers) >= 1)
fam_built <- gb$data[[tile_layers[1]]]
fam_left  <- min(fam_built$xmin, na.rm = TRUE)
fam_right <- max(fam_built$xmax, na.rm = TRUE)

# --- collapsed clade overlay using the FAMILY fill scale (no legend) ---
collapsed_rows_df <- lab_pos %>%
  dplyr::select(node, y) %>%
  dplyr::inner_join(clade_fam_tbl, by = "node") %>%
  dplyr::mutate(
    fam  = factor(fam, levels = names(family_palette)),
    x0   = fam_left, x1 = fam_right,
    ymin = y - 0.5,  ymax = y + 0.5
  )

p <- p + geom_rect(
  data = collapsed_rows_df,
  aes(xmin = x0, xmax = x1, ymin = ymin, ymax = ymax, fill = fam),
  inherit.aes = FALSE,
  color = NA,
  show.legend = FALSE
)

# --- seed ALL families into the FAMILY legend (still same fill scale) ---
legend_seed <- data.frame(
  fam = factor(names(family_palette), levels = names(family_palette)),
  x   = fam_left - 1,        # anywhere off the panel is fine
  y   = y_bottom - 1
)

p <- p + geom_point(
  data = legend_seed,
  aes(x = x, y = y, fill = fam),
  inherit.aes = FALSE,
  shape = 22, size = 3, alpha = 0,   # invisible on plot
  show.legend = TRUE                  # but registers every level with the guide
)

# >>> Now you can start the next panel with a fresh fill scale <<<
# p <- p + ggnewscale::new_scale_fill()


# 9b) ORIGIN heatmap + its own legend
p <- p + ggnewscale::new_scale_fill()
p <- gheatmap(
  p, origin_df,
  offset = as.numeric(panel_offsets["origin"]),
  width  = as.numeric(panel_widths["origin"]),
  colnames = FALSE, font.size = 3, hjust = 0, color = NA
)
p <- p + scale_fill_manual(
  values = origin_palette,
  breaks = names(origin_palette),
  drop = FALSE, na.value = "#F2F2F2",
  name = "Genome origin",
  guide = guide_legend(order = 3, ncol = 1, byrow = TRUE, title.position = "top")
)

# 9c) PHYLUM heatmap + its own legend
p <- p + ggnewscale::new_scale_fill()
p <- gheatmap(
  p, phylum_df,
  offset = as.numeric(panel_offsets["phylum"]),
  width  = as.numeric(panel_widths["phylum"]),
  colnames = FALSE, font.size = 3, hjust = 0, color = NA
)
p <- p + scale_fill_manual(
  values = phylum_color_map,
  breaks = names(phylum_color_map),
  drop = FALSE, na.value = "#F2F2F2",
  name = "Bacterial phylum",
  guide = guide_legend(order = 4, ncol = 1, byrow = TRUE, title.position = "top")
)

# 9d) CLASS heatmap + its own legend
p <- p + ggnewscale::new_scale_fill()
p <- gheatmap(
  p, class_df,
  offset = as.numeric(panel_offsets["class"]),
  width  = as.numeric(panel_widths["class"]),
  colnames = FALSE, font.size = 3, hjust = 0, color = NA
)
p <- p + scale_fill_manual(
  values = class_color_map,
  breaks = names(class_color_map),
  drop = FALSE, na.value = "#F2F2F2",
  name = "Bacterial class",
  guide = guide_legend(order = 5, ncol = 1, byrow = TRUE, title.position = "top")
)

# Make all legends appear as one vertical column on the right
p <- p + theme(
  legend.position = "right",
  legend.box = "vertical",
  legend.key.height = unit(4, "mm"),
  legend.key.width  = unit(4, "mm"),
  legend.title = element_text(size = 9, face = "bold"),
  legend.text  = element_text(size = 8)
)


# --- panel titles (rotated where requested) ---
# Build once to read the true x extents of each heatmap panel
# --- PANEL TITLES (rotated where requested) ---

# Build once to read the true x extents of each heatmap panel
gb <- ggplot_build(p)
tile_layers <- which(vapply(p$layers, function(L) inherits(L$geom, "GeomTile"), logical(1)))
stopifnot(length(tile_layers) >= 4)

# Centers of the four heatmap panels, in the order added
built_tiles <- gb$data[tile_layers[1:4]]
x_centers <- vapply(built_tiles, function(d) {
  (min(d$xmin, na.rm=TRUE) + max(d$xmax, na.rm=TRUE)) / 2
}, numeric(1))

# Give headroom and disable clipping so titles are not cropped
# Compute a smaller headroom above the top tip
y_max <- max(subset(p$data, isTip)$y, na.rm = TRUE)
tip_n <- Ntip(tr_pruned)

# TUNE THESE THREE NUMBERS to control top space:
TITLE_PAD_TIPS <- 15     # was 10% of tips; try 3–6
TOP_EXPAND     <- 0.1   # was 0.18; try 0.02–0.05
TOP_MARGIN_MM  <- 10      # was 20 mm; try 6–10

y_top <- y_max + TITLE_PAD_TIPS

p <- p +
  coord_cartesian(ylim = c(y_bottom_all, y_top + 0.5), clip = "off") +
  scale_y_continuous(expand = expansion(mult = c(0, TOP_EXPAND))) +
  theme(plot.margin = margin(t = TOP_MARGIN_MM, r = 6, b = 8, l = 6, unit = "mm"))

# (IMPORTANT) Remove this line from your code (it double-adds space):
# p <- p + expand_limits(y = y_top + 1)

# 4 heatmap panel titles (rotated 90°)
titles_df <- data.frame(
  x = x_centers,
  y = y_top,
  label = c("Crassvirales family", "Genome origin", "Bacterial phylum", "Bacterial class")
)
p <- p + geom_text(
  data = titles_df,
  aes(x = x, y = y, label = label),
  angle = 90, vjust = 0, hjust = 0.5,
  size = 3.2, fontface = "bold"
)

# Other panel titles
# Prophage coordinates (contig panel) — rotated 90°
p <- p + annotate(
  "text",
  x = contig_offset + contig_width/2,
  y = y_top,
  label = "Prophage coordinates, mb",
  angle = 90, vjust = 0, hjust = 0.5,
  size = 3.2, fontface = "bold"
)

# Prophage genomic map — horizontal
p <- p + annotate(
  "text",
  x = gene_offset + gene_width/2,
  y = y_top,
  label = "Prophage genomic map, kb",
  vjust = 0, hjust = 0.5,
  size = 3.2, fontface = "bold"
)

# Prophage length, kb — rotated 90°
p <- p + annotate(
  "text",
  x = len_offset + len_width/2,
  y = y_top,
  label = "Prophage length, kb",
  angle = 90, vjust = 0, hjust = 0.5,
  size = 3.2, fontface = "bold"
)

# Ensure the computed y_top is included in limits
p <- p + expand_limits(y = y_top + 1)

# ======================= draw contig panel =======================
# light background across the entire contig panel
# --- PROPHAGE COORDINATES (with legend) ---
p <- p + ggnewscale::new_scale_fill()   # fresh fill scale for this panel

# light background across the entire contig panel (no legend)
p <- p + geom_rect(
  data = contig_row_bg_df,
  aes(xmin = x0, xmax = x1, ymin = ymin, ymax = ymax),
  inherit.aes = FALSE,
  fill = "#F2F2F2", color = NA
)

# gray contig bars (legend label = "Bacterial contig")
p <- p + geom_rect(
  data = contig_bar_df,
  aes(xmin = x0, xmax = x1, ymin = ymin, ymax = ymax, fill = "Bacterial contig"),
  inherit.aes = FALSE,
  color = NA
)

# red prophage segment (legend label = "Prophage")
p <- p + geom_rect(
  data = proph_rect_df,
  aes(xmin = xr0, xmax = xr1, ymin = ymin, ymax = ymax, fill = "Prophage"),
  inherit.aes = FALSE,
  color = NA
)

# legend for this panel
p <- p + scale_fill_manual(
  values = c("Prophage" = "#d62728", "Bacterial contig" = "grey60"),
  breaks = c("Prophage", "Bacterial contig"),
  name   = "Prophage coordinates",
  guide  = guide_legend(order = 6, ncol = 1, byrow = TRUE, title.position = "top")
)


# ======================= draw length bars panel =======================
# background for all tip rows in the length panel
p <- p + geom_rect(
  data = len_bg_df,
  aes(xmin = x0, xmax = x1, ymin = ymin, ymax = ymax),
  inherit.aes = FALSE,
  fill = "#F2F2F2", color = NA
)

# actual bars (use a fixed fill to keep independent of the heatmap fill scale)
p <- p + geom_rect(
  data = bars_df,
  aes(xmin = x0, xmax = x1, ymin = ymin, ymax = ymax),
  inherit.aes = FALSE,
  fill = "grey35", color = NA
)

p <- p + geom_text(
  data = labels_df,
  aes(x = x_text, y = y, label = lab),
  inherit.aes = FALSE,
  hjust = 0, vjust = 0.5,
  size = 2.6
)

# --- NOW add the gene map panel to the right ---
# --- PROPHAGE GENOMIC MAP (with legend) ---
p <- p + ggnewscale::new_scale_fill()   # fresh fill scale for this panel

# background strip (no legend)
p <- p + geom_rect(
  data = gene_row_bg_df,
  aes(xmin = x0, xmax = x1, ymin = ymin, ymax = ymax),
  inherit.aes = FALSE, fill = "#F2F2F2", color = NA
)

# gene arrows with legend
p <- p + geom_polygon(
  data = gene_polys,
  aes(x = x, y = y, group = id, fill = gene_type),
  inherit.aes = FALSE,
  color = NA
)

p <- p + scale_fill_manual(
  values = c("Crassphage gene" = "#d62728", "Hypothetical protein" = "grey30"),
  breaks = c("Crassphage gene", "Hypothetical protein"),
  name   = "Prophage genomic map",
  guide  = guide_legend(order = 7, ncol = 1, byrow = TRUE, title.position = "top")
)

cat("First 10 rows with non-NA origin:\n")
print(which(!is.na(origin_df[[1]]))[1:10])

# ======================= 10) save =======================
ggsave(out_png, p, width=18, height=16, dpi=300)
ggsave(out_svg, p, width=18, height=16)
ggsave(out_pdf, p, width=18, height=16)
message(sprintf("Saved plots:\n  %s\n  %s", out_png, out_pdf))
