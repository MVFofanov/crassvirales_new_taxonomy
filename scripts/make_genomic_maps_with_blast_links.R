suppressPackageStartupMessages({
  library(tidyverse)
  library(gggenomes)
  library(cowplot)
})

# ============================================================
#  Prophage vs closest bacterial fragment summary figure (2x8)
#  - genes from PHOLD table (colored by functional category)
#  - links from BLASTn (filtered by min aligned length)
#  - pairs are read from user-provided table (order preserved)
# ============================================================

# ------------------ SETTINGS ------------------

WD <- "C:/crassvirales/crassvirales_new_taxonomy/crassvirales_prophages/integrated_prophages/bacterial_flank_analysis"

PHOLD_TSV <- file.path(WD, "all_phold_per_cds_predictions.tsv")
BLAST_TSV <- file.path(WD, "all_vs_all.blastn.tsv")

# NEW: your pairs table (order defines panel order)
PAIRS_TSV <- file.path(WD, "pairs_prophage_with_flanks_vs_bacterial.tsv")

OUT_PNG <- file.path(WD, "prophage_vs_bacterial_2x8.png")
OUT_PDF <- file.path(WD, "prophage_vs_bacterial_2x8.pdf")

N_COL <- 2

# BLAST filters
MIN_ALIGN_LEN <- 200   # show only BLAST homology blocks >= 200 bp
MIN_PID <- 0           # set e.g. 70 if you want stricter %identity

# Self-links (within same contig). Usually OFF for 16-panel summary.
SHOW_SELF_LINKS <- FALSE

# Link styling
COLOR_LINKS_BY_PID <- TRUE

# Make gene arrows more visible
GENE_HEIGHT <- 0.25    # increase if arrows still too small (e.g. 0.35)

FUNC_COLORS <- c(
  "connector" = "#5A5A5A",
  "DNA, RNA and nucleotide metabolism" = "#f000ff",
  "head and packaging" = "#ff008d",
  "integration and excision" = "#E0B0FF",
  "lysis" = "#001eff",
  "moron, auxiliary metabolic gene and host takeover" = "#8900ff",
  "other" = "#4deeea",
  "tail" = "#74ee15",
  "transcription regulation" = "#ffe700",
  "unknown" = "#AAAAAA"
)

# ------------------ READ INPUTS ------------------

# PHOLD
phold <- read_tsv(PHOLD_TSV, show_col_types = FALSE) %>%
  rename(func = .data[["function"]]) %>%
  mutate(
    start = readr::parse_integer(as.character(start)),
    end   = readr::parse_integer(as.character(end)),
    gene_min = pmin(start, end),
    gene_max = pmax(start, end),
    
    strand_num = case_when(
      strand %in% c("+", "plus", 1, "1") ~  1L,
      strand %in% c("-", "minus", -1, "-1") ~ -1L,
      TRUE ~ NA_integer_
    ),
    
    func = case_when(
      is.na(func) | func == "" ~ "unknown",
      func == "unknown function" ~ "unknown",
      TRUE ~ func
    )
  )



# BLAST
blast <- read_tsv(
  BLAST_TSV,
  col_names = c("qseqid","sseqid","qstart","qend","sstart","send",
                "pident","length","bitscore","evalue","qlen","slen","sstrand"),
  show_col_types = FALSE
) %>%
  mutate(
    qstart = as.integer(qstart),
    qend   = as.integer(qend),
    sstart = as.integer(sstart),
    send   = as.integer(send),
    pident = as.numeric(pident),
    length = as.integer(length),
    bitscore = as.numeric(bitscore),
    evalue = as.numeric(evalue),
    qlen = as.integer(qlen),
    slen = as.integer(slen),
    qmin = pmin(qstart, qend),
    qmax = pmax(qstart, qend),
    smin = pmin(sstart, send),
    smax = pmax(sstart, send)
  )

# ------------------ READ PAIRS (ORDER PRESERVED) ------------------

pairs <- read_tsv(PAIRS_TSV, show_col_types = FALSE) %>%
  rename(
    prophage_id = prophage,
    bacterial_id = bacterial
  ) %>%
  mutate(pair_index = row_number())

stopifnot(nrow(pairs) == 16)

message("Loaded pairs: ", nrow(pairs))
print(pairs)

# ------------------ HELPERS ------------------

get_seq_lengths <- function(ids, phold_df, blast_df) {
  seqs_from_blast <- blast_df %>%
    filter(qseqid %in% ids | sseqid %in% ids) %>%
    transmute(seq_id = qseqid, len = qlen) %>%
    bind_rows(blast_df %>% transmute(seq_id = sseqid, len = slen)) %>%
    filter(seq_id %in% ids) %>%
    group_by(seq_id) %>%
    summarize(len = suppressWarnings(max(len, na.rm = TRUE)), .groups = "drop") %>%
    filter(is.finite(len))
  
  seqs_from_genes <- phold_df %>%
    filter(contig_id %in% ids) %>%
    group_by(contig_id) %>%
    summarize(len = suppressWarnings(max(pmax(start, end), na.rm = TRUE)), .groups = "drop") %>%
    transmute(seq_id = contig_id, len = len) %>%
    filter(is.finite(len))
  
  bind_rows(seqs_from_blast, seqs_from_genes) %>%
    group_by(seq_id) %>%
    summarize(len = suppressWarnings(max(len, na.rm = TRUE)), .groups = "drop") %>%
    filter(is.finite(len)) %>%
    rename(length = len)
}

make_links_schema <- function(blast_df, prophage_id, bacterial_id, min_align_len, min_pid) {
  blast_df %>%
    filter(
      (qseqid == prophage_id & sseqid == bacterial_id) |
        (qseqid == bacterial_id & sseqid == prophage_id)
    ) %>%
    filter(length >= min_align_len, pident >= min_pid) %>%
    mutate(
      # normalize so seq_id is prophage, seq_id2 is bacterial
      seq_id  = if_else(qseqid == prophage_id, qseqid, sseqid),
      seq_id2 = if_else(qseqid == prophage_id, sseqid, qseqid),
      
      start  = if_else(qseqid == prophage_id, qmin, smin),
      end    = if_else(qseqid == prophage_id, qmax, smax),
      
      start2 = if_else(qseqid == prophage_id, smin, qmin),
      end2   = if_else(qseqid == prophage_id, smax, qmax)
    ) %>%
    transmute(seq_id, start, end, seq_id2, start2, end2, pident)
}

make_self_links_schema <- function(blast_df, ids, min_align_len, min_pid) {
  blast_df %>%
    filter(qseqid == sseqid, qseqid %in% ids) %>%
    filter(length >= min_align_len, pident >= min_pid) %>%
    transmute(
      seq_id  = qseqid,
      start   = qmin,
      end     = qmax,
      seq_id2 = sseqid,
      start2  = smin,
      end2    = smax,
      pident  = pident
    )
}

make_pair_plot <- function(prophage_id, bacterial_id,
                           phold_df, blast_df,
                           title_txt = NULL,
                           min_align_len = 200,
                           min_pid = 0,
                           show_self_links = FALSE,
                           color_links_by_pid = TRUE,
                           func_colors = FUNC_COLORS,
                           gene_height = GENE_HEIGHT) {
  
  ids <- c(prophage_id, bacterial_id)
  
  # genes (IMPORTANT: use numeric strand column strand_num!)
  genes <- phold_df %>%
    filter(contig_id %in% ids) %>%
    filter(!is.na(gene_min), !is.na(gene_max), !is.na(strand_num)) %>%
    transmute(
      seq_id  = contig_id,
      gene_id = cds_id,
      start   = gene_min,
      end     = gene_max,
      strand  = strand_num,   # <<<<<< THIS makes arrows appear
      func    = func
    )
  
  # seqs with plotting order (prophage top, bacterial bottom)
  seqs <- get_seq_lengths(ids, phold_df, blast_df) %>%
    mutate(bin_id = factor(seq_id, levels = c(prophage_id, bacterial_id))) %>%
    arrange(bin_id)
  
  # links
  links_pb   <- make_links_schema(blast_df, prophage_id, bacterial_id, min_align_len, min_pid)
  links_self <- make_self_links_schema(blast_df, ids, min_align_len, min_pid)
  links <- if (show_self_links) bind_rows(links_pb, links_self) else links_pb
  
  if (is.null(title_txt)) title_txt <- prophage_id
  
  # ---- base: seq tracks only ----
  p <- gggenomes(
    seqs  = seqs %>% transmute(seq_id, length, bin_id),
    genes = genes,
    links = links
  ) +
    geom_seq(linewidth = 0.25)
  
  # ---- add links FIRST (so genes are drawn on top later) ----
  if (color_links_by_pid) {
    p <- p +
      geom_link(aes(color = pident), alpha = 0.35, linewidth = 0.25) +
      scale_color_viridis_c(option = "C", end = 0.95, guide = ggplot2::guide_none())
  } else {
    p <- p + geom_link(alpha = 0.30, linewidth = 0.25, color = "grey35")
  }
  
  # ---- add genes LAST (on top of links) + styling ----
  p <- p +
    geom_gene(aes(fill = func), color = "black", linewidth = 0.12, height = gene_height) +
    scale_fill_manual(values = func_colors, drop = FALSE, guide = ggplot2::guide_none()) +
    labs(title = title_txt) +
    theme_bw(base_size = 8) +
    theme(
      plot.title = element_text(size = 7.5, face = "bold"),
      legend.position = "none",
      panel.grid = element_blank(),
      axis.title = element_blank()
    ) +
    guides(fill = ggplot2::guide_none(), color = ggplot2::guide_none())
  
  p
}


# ------------------ BUILD PANELS (ORDER = PAIRS FILE) ------------------

# You want first 8 in column 1, next 8 in column 2.
# cowplot fills row-wise by default, so we reorder to achieve column-wise filling:
# desired layout:
#   col1: 1..8
#   col2: 9..16
# row-wise order should be: 1,9,2,10,3,11,...,8,16
plot_order <- c(rbind(1:8, 9:16))

pairs_ordered <- pairs %>% slice(plot_order)

plots <- purrr::pmap(
  list(pairs_ordered$prophage_id, pairs_ordered$bacterial_id, pairs_ordered$pair_index),
  function(proph, bact, i) {
    make_pair_plot(
      prophage_id = proph,
      bacterial_id = bact,
      phold_df = phold,
      blast_df = blast,
      title_txt = paste0("Pair ", i),
      min_align_len = MIN_ALIGN_LEN,
      min_pid = MIN_PID,
      show_self_links = SHOW_SELF_LINKS,
      color_links_by_pid = COLOR_LINKS_BY_PID,
      func_colors = FUNC_COLORS,
      gene_height = GENE_HEIGHT
    )
  }
)

final <- cowplot::plot_grid(plotlist = plots, ncol = N_COL, align = "hv")

# ------------------ EXPORT ------------------

ggsave(OUT_PDF, final, width = 12, height = 16, device = cairo_pdf)
ggsave(OUT_PNG, final, width = 12, height = 16, dpi = 600)

message("Saved: ", OUT_PNG)
message("Saved: ", OUT_PDF)
