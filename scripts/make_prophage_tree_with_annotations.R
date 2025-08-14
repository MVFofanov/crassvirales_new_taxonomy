#!/usr/bin/env Rscript

# ======================= deps =======================
pkgs_cran <- c("ape","ggplot2","dplyr","readr","tibble","tidyr","stringr")
for (pkg in pkgs_cran) if (!requireNamespace(pkg, quietly=TRUE)) install.packages(pkg, repos="https://cloud.r-project.org")
if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager", repos="https://cloud.r-project.org")
for (pkg in c("ggtree","treeio")) if (!requireNamespace(pkg, quietly=TRUE)) BiocManager::install(pkg, ask=FALSE, update=FALSE)
if (!requireNamespace("ggnewscale", quietly=TRUE)) install.packages("ggnewscale", repos="https://cloud.r-project.org")

library(ape)
library(ggplot2)
library(ggtree)
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
  "Bacilli"="#e6194B","Clostridia"="#f58231","Erysipelotrichia"="#911eb4",
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

out_png <- file.path(proj_dir, "TerL_Crassvirales_pruned_collapsed.png")
out_svg <- file.path(proj_dir, "TerL_Crassvirales_pruned_collapsed.svg")
out_pdf <- file.path(proj_dir, "TerL_Crassvirales_pruned_collapsed.pdf")

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
    genome_id               = as.character(genome_id),
    protein_id              = as.character(protein_id)
  ) %>%
  select(label, prophage_family, origin, bacterial_phylum, bacterial_class,
         bacterial_order_family, bacterial_tax_lineage, bacterial_id,
         prophage_length, bacterial_contig_length, genome_id, protein_id)

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

# ======================= 6) geometry for labels & panels =======================
# preview (collapsed) to measure in data units
p0 <- ggtree(tr_pruned, layout="rectangular") %<+% tips_tbl_pruned
for (nd in collapse_nodes) p0 <- collapse(p0, node=nd)
x_span <- diff(range(p0$data$x, na.rm=TRUE))

lab_offset  <- 0.03 * x_span   # distance between tips and labels
panel_shift <- 0.8 * x_span   # push block of panels to the right of labels, 0.060 * x_span
panel_gap   <- 0.15 * x_span   # uniform gap between panels, 0.012 * x_span
panel_w     <- 0.045 * x_span   # SAME width for all panels (change this one number), 0.045 * x_span 

## ---- panel layout (now with a 'len' panel) ----
# widths: keep existing four panels at panel_w; give length panel a bit more room
len_width <- 0.06 * x_span  # tweak to taste

# ---- panel layout: 4 heatmap panels only ----
panel_names  <- c("family","origin","phylum","class")
panel_widths <- setNames(rep(panel_w, length(panel_names)), panel_names)

cum_left       <- c(0, cumsum(head(panel_widths, -1) + panel_gap))
panel_offsets  <- lab_offset + panel_shift + cum_left
names(panel_offsets) <- panel_names

# total width of the 4 heatmap panels
panels_total_width <- sum(panel_widths) + panel_gap * (length(panel_widths) - 1)

# bar region: sits to the RIGHT of the panels
bar_gap   <- 1.2 * x_span     # space between last heatmap panel and bars, 0.06 * x_span
bar_width <- 0.2 * x_span     # width of the barplot strip (tweak to taste), 0.08 * x_span 

# left edge of the bar region
len_offset <- lab_offset + panel_shift + panels_total_width + bar_gap
len_width  <- bar_width

# total width added to the tree canvas
total_width_all <- panels_total_width + bar_gap + bar_width
extra_right     <- 0.10 * x_span

# ======================= 7) base tree =======================
p <- ggtree(tr_pruned, layout="rectangular") %<+% tips_tbl_pruned +
  theme_tree2() +
  geom_tippoint(aes(color = family), size = 1.2, alpha = 0.9) +
  ggtitle("TerL — Crassvirales clade (pruned), non-prophage clades collapsed",
          subtitle = "Triangles = collapsed clades; numbers = # tips inside")

for (nd in collapse_nodes) p <- collapse(p, node = nd)

pd <- p$data
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
  scale_color_manual(values = family_palette, na.value = "#BDBDBD") +
  guides(color = guide_legend(title = "Crassvirales family (tips)",
                              override.aes = list(size = 3, alpha = 1), order = 1)) +
  xlim_tree(max(p$data$x, na.rm=TRUE) + lab_offset + panel_shift + total_width_all + extra_right)

# ======================= 8) side-table (prophages only; exact tip order) =======================
side_df <- tips_tbl_pruned %>%
  transmute(
    label,
    `Crassvirales family` = ifelse(is_prophage, as.character(family),           NA_character_),
    `Genome origin`       = ifelse(is_prophage, origin,                         NA_character_),
    `Bacterial phylum`    = ifelse(is_prophage, as.character(bacterial_phylum), NA_character_),
    `Bacterial class`     = ifelse(is_prophage, as.character(bacterial_class),  NA_character_)
  ) %>%
  right_join(tibble(label = tr_pruned$tip.label), by = "label") %>%
  arrange(match(label, tr_pruned$tip.label)) %>%
  tibble::column_to_rownames("label")

# diagnostics
cat("Non-NA cells in panels:\n"); print(colSums(!is.na(side_df)))

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

# ======================= 9) draw panels (robust, one combined scale) =======================

# helper: add panel-specific prefixes but keep NAs as NA
add_prefix <- function(x, pref) ifelse(is.na(x), NA_character_, paste0(pref, as.character(x)))

# make prefixed copies of each panel column
family_df_pref <- family_df; family_df_pref[[1]] <- add_prefix(family_df_pref[[1]], "FAM:")
origin_df_pref <- origin_df; origin_df_pref[[1]] <- add_prefix(origin_df_pref[[1]], "ORG:")
phylum_df_pref <- phylum_df; phylum_df_pref[[1]] <- add_prefix(phylum_df_pref[[1]], "PHY:")
class_df_pref  <- class_df;  class_df_pref[[1]]  <- add_prefix(class_df_pref[[1]],  "CLS:")

# one big palette = union of all palettes with matching prefixes
panel_palette <- c(
  setNames(family_palette, paste0("FAM:", names(family_palette))),
  setNames(origin_palette, paste0("ORG:", names(origin_palette))),
  setNames(phylum_color_map, paste0("PHY:", names(phylum_color_map))),
  setNames(class_color_map,  paste0("CLS:", names(class_color_map)))
)

# draw the four panels (no ggnewscale; same offsets/widths you computed)
p <- gheatmap(
  p, family_df_pref,
  offset = as.numeric(panel_offsets["family"]),
  width  = as.numeric(panel_widths["family"]),
  colnames = TRUE, colnames_position = "top",
  font.size = 3, hjust = 0, color = NA
)
p <- gheatmap(
  p, origin_df_pref,
  offset = as.numeric(panel_offsets["origin"]),
  width  = as.numeric(panel_widths["origin"]),
  colnames = TRUE, colnames_position = "top",
  font.size = 3, hjust = 0, color = NA
)
p <- gheatmap(
  p, phylum_df_pref,
  offset = as.numeric(panel_offsets["phylum"]),
  width  = as.numeric(panel_widths["phylum"]),
  colnames = TRUE, colnames_position = "top",
  font.size = 3, hjust = 0, color = NA
)
p <- gheatmap(
  p, class_df_pref,
  offset = as.numeric(panel_offsets["class"]),
  width  = as.numeric(panel_widths["class"]),
  colnames = TRUE, colnames_position = "top",
  font.size = 3, hjust = 0, color = NA
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

# (optional) a small title above the bar panel — comment out if not needed
# p <- p + annotate(
#   "text",
#   x = len_offset, y = max(pd_tips$y, na.rm = TRUE) + 2,
#   label = "Prophage length (kb)",
#   hjust = 0, size = 3
# )

# a single fill scale that knows about all prefixed values
p <- p + scale_fill_manual(
  values = panel_palette,
  na.value = "#F2F2F2",
  drop = FALSE,
  name = "Side panels"
)



cat("First 10 rows with non-NA origin:\n")
print(which(!is.na(origin_df[[1]]))[1:10])

# ======================= 10) save =======================
ggsave(out_png, p, width=18, height=14, dpi=300)
ggsave(out_svg, p, width=18, height=14)
ggsave(out_pdf, p, width=18, height=14)
message(sprintf("Saved plots:\n  %s\n  %s", out_png, out_pdf))
