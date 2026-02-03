suppressPackageStartupMessages({
  library(tidyverse)
  library(gggenomes)
  library(cowplot)
})

# ------------------ SETTINGS ------------------

WD <- "C:/crassvirales/crassvirales_new_taxonomy/crassvirales_prophages/integrated_prophages/bacterial_flank_analysis"

PHOLD_TSV <- file.path(WD, "all_phold_per_cds_predictions.tsv")
BLAST_TSV <- file.path(WD, "all_vs_all.blastn.tsv")
PAIRS_TSV <- file.path(WD, "pairs_prophage_with_flanks_vs_bacterial.tsv")

OUT_PNG <- file.path(WD, "prophage_vs_bacterial_2x8.png")
OUT_PDF <- file.path(WD, "prophage_vs_bacterial_2x8.pdf")

DEBUG_ONE <- file.path(WD, "debug_one_panel.png")

N_COL <- 2

MIN_ALIGN_LEN <- 200
MIN_PID <- 0

SHOW_SELF_LINKS <- FALSE
COLOR_LINKS_BY_PID <- TRUE

GENE_HEIGHT <- 0.7     # visual height of arrows (y units)
SEQ_LINEWIDTH <- 0.35  # thickness of seq baseline
LINK_LINEWIDTH <- 0.25
GENE_OUTLINE <- 0.20

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

phold <- read_tsv(PHOLD_TSV, show_col_types = FALSE) %>%
  rename(func = `function`) %>%
  mutate(
    start = as.numeric(as.character(start)),
    end   = as.numeric(as.character(end))
  ) %>%
  filter(!is.na(start), !is.na(end)) %>%
  mutate(
    start = as.integer(round(start)),
    end   = as.integer(round(end)),
    gene_min = pmin(start, end),
    gene_max = pmax(start, end),
    strand_num = case_when(
      strand == "+" ~  1L,
      strand == "-" ~ -1L,
      TRUE ~ NA_integer_
    ),
    func = case_when(
      is.na(func) | func == "" ~ "unknown",
      func == "unknown function" ~ "unknown",
      TRUE ~ func
    )
  )

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
    qlen = as.integer(qlen),
    slen = as.integer(slen),
    qmin = pmin(qstart, qend),
    qmax = pmax(qstart, qend),
    smin = pmin(sstart, send),
    smax = pmax(sstart, send)
  )

pairs <- read_tsv(PAIRS_TSV, show_col_types = FALSE) %>%
  rename(prophage_id = prophage, bacterial_id = bacterial) %>%
  mutate(pair_index = row_number())

stopifnot(nrow(pairs) == 16)

# ------------------ HELPERS ------------------

make_links_schema <- function(blast_df, prophage_id, bacterial_id, min_align_len, min_pid) {
  blast_df %>%
    filter(
      (qseqid == prophage_id & sseqid == bacterial_id) |
        (qseqid == bacterial_id & sseqid == prophage_id)
    ) %>%
    filter(length >= min_align_len, pident >= min_pid) %>%
    mutate(
      seq_id  = if_else(qseqid == prophage_id, qseqid, sseqid),
      seq_id2 = if_else(qseqid == prophage_id, sseqid, qseqid),
      start   = if_else(qseqid == prophage_id, qmin, smin),
      end     = if_else(qseqid == prophage_id, qmax, smax),
      start2  = if_else(qseqid == prophage_id, smin, qmin),
      end2    = if_else(qseqid == prophage_id, smax, qmax)
    ) %>%
    transmute(
      seq_id  = as.character(seq_id),
      start   = as.integer(start),
      end     = as.integer(end),
      seq_id2 = as.character(seq_id2),
      start2  = as.integer(start2),
      end2    = as.integer(end2),
      pident
    )
}

make_self_links_schema <- function(blast_df, ids, min_align_len, min_pid) {
  blast_df %>%
    filter(qseqid == sseqid, qseqid %in% ids) %>%
    filter(length >= min_align_len, pident >= min_pid) %>%
    transmute(
      seq_id  = as.character(qseqid),
      start   = as.integer(qmin),
      end     = as.integer(qmax),
      seq_id2 = as.character(sseqid),
      start2  = as.integer(smin),
      end2    = as.integer(smax),
      pident
    )
}

make_pair_plot <- function(prophage_id, bacterial_id,
                           phold_df, blast_df,
                           title_txt = NULL,
                           min_align_len = 200,
                           min_pid = 0,
                           show_self_links = FALSE,
                           color_links_by_pid = TRUE,
                           func_colors = FUNC_COLORS) {
  
  ids <- c(prophage_id, bacterial_id)
  
  # --- genes (PHOLD) ---
  genes <- phold_df %>%
    filter(contig_id %in% ids) %>%
    filter(!is.na(gene_min), !is.na(gene_max)) %>%
    transmute(
      seq_id  = as.character(contig_id),
      start   = as.integer(gene_min),
      end     = as.integer(gene_max),
      strand  = as.character(strand),   # "+" / "-"
      type    = "CDS",                  # IMPORTANT for geom_gene()
      func    = func
    )
  
  if (nrow(genes) == 0) stop("No genes for pair: ", prophage_id, " vs ", bacterial_id)
  
  # --- seqs from genes ---
  seqs <- genes %>%
    group_by(seq_id) %>%
    summarize(length = max(end, na.rm = TRUE), .groups = "drop") %>%
    mutate(
      length = as.integer(pmax(length, 1L)),
      bin_id = factor(seq_id, levels = c(prophage_id, bacterial_id))
    ) %>%
    arrange(bin_id)
  
  # --- links ---
  links_pb <- make_links_schema(blast_df, prophage_id, bacterial_id, min_align_len, min_pid)
  links_self <- make_self_links_schema(blast_df, ids, min_align_len, min_pid)
  links <- if (show_self_links) bind_rows(links_pb, links_self) else links_pb
  
  if (is.null(title_txt)) title_txt <- prophage_id
  
  p <- gggenomes(seqs = seqs, genes = genes, links = links) +
    geom_seq(linewidth = SEQ_LINEWIDTH)
  
  # links (geom_link is a filled polygon; use fill mapping)
  if (nrow(links) > 0) {
    if (color_links_by_pid) {
      p <- p +
        geom_link(aes(color = pident), alpha = 0.25, linewidth = LINK_LINEWIDTH) +
        scale_color_viridis_c(option = "C", end = 0.95, guide = "none")
    } else {
      p <- p + geom_link(color = "grey50", alpha = 0.20, linewidth = LINK_LINEWIDTH)
    }
  }
  
  # genes: make them visible (separate strands + bigger size)
  p <- p +
    geom_gene(aes(fill = func),
              position = "strand",   # THIS is the big visibility fix
              size = 2.2,            # increase if still too thin (try 3–4)
              stroke = GENE_OUTLINE, # outline thickness
              colour = "black") +
    scale_fill_manual(values = func_colors, drop = FALSE, guide = "none") +
    labs(title = title_txt) +
    theme_bw(base_size = 10) +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      legend.position = "none",
      panel.grid = element_blank(),
      axis.title = element_blank()
    ) +
    coord_cartesian(clip = "off")
  
  p
}



# ------------------ DEBUG ONE PANEL FIRST ------------------

p1 <- make_pair_plot(
  prophage_id  = pairs$prophage_id[1],
  bacterial_id = pairs$bacterial_id[1],
  phold_df = phold,
  blast_df = blast,
  title_txt = pairs$prophage_id[1],
  min_align_len = MIN_ALIGN_LEN,
  min_pid = MIN_PID,
  show_self_links = SHOW_SELF_LINKS,
  color_links_by_pid = COLOR_LINKS_BY_PID,
  func_colors = FUNC_COLORS
)

p1 <- make_pair_plot(
  prophage_id  = pairs$prophage_id[1],
  bacterial_id = pairs$bacterial_id[1],
  phold_df = phold,
  blast_df = blast
)

p1 %>% track_info
nrow(p1 %>% pull_genes())   # should be > 0
head(p1 %>% pull_genes())

ggsave(DEBUG_ONE, p1, width = 16, height = 4, dpi = 300)
message("Saved debug panel: ", DEBUG_ONE)

# ------------------ BUILD PANELS ------------------

plot_order <- c(rbind(1:8, 9:16))
pairs_ordered <- pairs %>% slice(plot_order)

plots <- purrr::pmap(
  list(pairs_ordered$prophage_id, pairs_ordered$bacterial_id),
  function(proph, bact) {
    make_pair_plot(
      prophage_id = proph,
      bacterial_id = bact,
      phold_df = phold,
      blast_df = blast,
      title_txt = proph,
      min_align_len = MIN_ALIGN_LEN,
      min_pid = MIN_PID,
      show_self_links = SHOW_SELF_LINKS,
      color_links_by_pid = COLOR_LINKS_BY_PID,
      func_colors = FUNC_COLORS
    )
  }
)

final <- cowplot::plot_grid(plotlist = plots, ncol = N_COL, align = "hv")

ggsave(OUT_PDF, final, width = 14, height = 30, device = cairo_pdf)
ggsave(OUT_PNG, final, width = 14, height = 30, dpi = 600)

message("Saved: ", OUT_PNG)
message("Saved: ", OUT_PDF)
