suppressPackageStartupMessages({
  library(tidyverse)
  library(gggenomes)
  library(cowplot)
  library(svglite)
})

# ------------------ SETTINGS ------------------

WD <- "C:/crassvirales/crassvirales_new_taxonomy/crassvirales_prophages/integrated_prophages/bacterial_flank_analysis"

PHOLD_TSV <- file.path(WD, "all_phold_per_cds_predictions.tsv")
BLAST_TSV <- file.path(WD, "all_vs_all.blastn.tsv")
PAIRS_TSV <- file.path(WD, "pairs_prophage_with_flanks_vs_bacterial.tsv")
PROPHAGE_TSV <- file.path(WD, "Crassphage_prophage_analysis_annotation_ncbi_with_samples_and_countries_crassvirales_clade_flanks.tsv")
VIRUS_SUMMARY_TSV <- file.path(WD, "all_virus_summary.tsv")

OUT_PNG <- file.path(WD, "prophage_vs_bacterial_2x8.png")
OUT_SVG <- file.path(WD, "prophage_vs_bacterial_2x8.svg")
OUT_PDF <- file.path(WD, "prophage_vs_bacterial_2x8.pdf")

DEBUG_ONE <- file.path(WD, "debug_one_panel.png")

N_COL <- 2

MIN_ALIGN_LEN <- 200
MIN_PID <- 0

SHOW_SELF_LINKS <- FALSE
COLOR_LINKS_BY_PID <- TRUE

GENE_HEIGHT <- 0.8     # visual height of arrows (y units)
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

# box colors (edit as you like)
BOX_COLORS <- c(
  "Crassvirales" = "#ff7f00",
  "Other"        = "#FFA500"
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
  rename(
    prophage_fragment_id = prophage,
    bacterial_fragment_id = bacterial
  )


# stopifnot(nrow(pairs) == 16)

virus_sum <- read_tsv(VIRUS_SUMMARY_TSV, show_col_types = FALSE) %>%
  mutate(
    # fragment id is left part before "|"
    fragment_id = str_split_fixed(seq_name, "\\|", 2)[,1],
    # prefer the explicit coordinates column if present
    coord_txt   = as.character(coordinates),
    # fallback: parse from seq_name suffix "provirus_123_456"
    suffix_txt  = str_split_fixed(seq_name, "\\|", 2)[,2]
  ) %>%
  mutate(
    # parse "start-end" from coordinates
    box_start = suppressWarnings(as.integer(str_match(coord_txt, "^(\\d+)-(\\d+)$")[,2])),
    box_end   = suppressWarnings(as.integer(str_match(coord_txt, "^(\\d+)-(\\d+)$")[,3]))
  ) %>%
  mutate(
    # fallback parse from suffix (e.g. "provirus_375_106207") if coordinates missing
    box_start = if_else(is.na(box_start),
                        suppressWarnings(as.integer(str_match(suffix_txt, "provirus_(\\d+)[_-](\\d+)")[,2])),
                        box_start),
    box_end   = if_else(is.na(box_end),
                        suppressWarnings(as.integer(str_match(suffix_txt, "provirus_(\\d+)[_-](\\d+)")[,3])),
                        box_end)
  ) %>%
  filter(!is.na(fragment_id), !is.na(box_start), !is.na(box_end)) %>%
  mutate(
    box_min = pmin(box_start, box_end),
    box_max = pmax(box_start, box_end),
    box_group = if_else(str_detect(taxonomy, "Crassvirales"), "Crassvirales", "Other"),
    box_color = unname(BOX_COLORS[box_group])
  ) %>%
  select(fragment_id, box_min, box_max, box_group, box_color, taxonomy)


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

parse_fragment_id <- function(x) {
  m <- stringr::str_match(x, "^(.+)_([0-9]+)-([0-9]+)$")
  if (any(is.na(m))) stop("Cannot parse fragment id: ", x)
  tibble(
    seq_id = x,
    base_contig = m[,2],
    frag_start = as.integer(m[,3]),
    frag_end   = as.integer(m[,4])
  )
}

# calc_prophage_boxes_for_ids <- function(ids, proph_tbl) {
#   frag_meta <- purrr::map_dfr(ids, parse_fragment_id) %>%
#     filter(!is.na(bacterial_id), !is.na(frag_start), !is.na(frag_end)) %>%
#     mutate(
#       frag_min = pmin(frag_start, frag_end),
#       frag_max = pmax(frag_start, frag_end)
#     )
#   
#   if (nrow(frag_meta) == 0) return(tibble())
#   
#   # join prophages by bacterial contig and keep overlaps with fragment
#   boxes <- frag_meta %>%
#     left_join(proph_tbl, by = "bacterial_id") %>%
#     filter(!is.na(prophage_min), !is.na(prophage_max)) %>%
#     filter(prophage_max >= frag_min, prophage_min <= frag_max) %>%  # overlap
#     mutate(
#       abs_box_start = pmax(prophage_min, frag_min),
#       abs_box_end   = pmin(prophage_max, frag_max),
#       rel_start = abs_box_start - frag_min + 1L,
#       rel_end   = abs_box_end   - frag_min + 1L,
#       rel_start = pmax(rel_start, 1L),
#       rel_end   = pmax(rel_end, 1L)
#     ) %>%
#     transmute(seq_id, rel_start, rel_end, genome_id)
#   
#   boxes
# }

calc_virus_boxes_for_fragments <- function(ids, virus_sum_tbl) {
  virus_sum_tbl %>%
    filter(fragment_id %in% ids) %>%
    transmute(
      seq_id = fragment_id,
      box_start = as.integer(box_min),
      box_end   = as.integer(box_max),
      box_group,
      box_color
    )
}


calc_crop_from_links <- function(links_pb, pad = 2000) {
  # links_pb has columns: start2/end2 on bacterial, start/end on prophage
  if (is.null(links_pb) || nrow(links_pb) == 0) return(NULL)
  
  lo <- min(pmin(links_pb$start2, links_pb$end2), na.rm = TRUE)
  hi <- max(pmax(links_pb$start2, links_pb$end2), na.rm = TRUE)
  
  tibble(
    crop_start = max(1L, as.integer(lo - pad)),
    crop_end   = as.integer(hi + pad)
  )
}


make_pair_plot <- function(prophage_id, bacterial_id,
                           phold_df, blast_df, virus_sum_tbl,
                           title_txt = NULL,
                           min_align_len = 200,
                           min_pid = 0,
                           show_self_links = FALSE,
                           color_links_by_pid = TRUE,
                           func_colors = FUNC_COLORS,
                           crop_pad = 10000) {
  
  ids <- c(prophage_id, bacterial_id)
  
  # --- links FIRST (these define what you actually plot) ---
  links_pb <- make_links_schema(blast_df, prophage_id, bacterial_id, min_align_len, min_pid)
  
  # --- crop window from plotted links (bacterial coords are start2/end2) ---
  crop_win <- calc_crop_from_links(links_pb, pad = crop_pad)
  
  # --- genes (PHOLD) ---
  genes_raw <- phold_df %>%
    filter(contig_id %in% ids) %>%
    filter(!is.na(gene_min), !is.na(gene_max)) %>%
    transmute(
      seq_id = as.character(contig_id),
      gene_min = as.integer(gene_min),
      gene_max = as.integer(gene_max),
      strand = as.character(strand),
      func = func
    )
  
  if (nrow(genes_raw) == 0) stop("No genes for pair: ", prophage_id, " vs ", bacterial_id)
  
  # --- crop/shift bacterial genes to crop window ---
  genes <- genes_raw %>%
    mutate(is_bacterial = seq_id == bacterial_id) %>%
    {
      if (!is.null(crop_win)) {
        filter(., !is_bacterial | (gene_max >= crop_win$crop_start & gene_min <= crop_win$crop_end))
      } else .
    } %>%
    mutate(
      start = if_else(is_bacterial & !is.null(crop_win),
                      gene_min - crop_win$crop_start + 1L,
                      gene_min),
      end   = if_else(is_bacterial & !is.null(crop_win),
                      gene_max - crop_win$crop_start + 1L,
                      gene_max)
    ) %>%
    transmute(
      seq_id,
      start = as.integer(start),
      end   = as.integer(end),
      strand,
      type = "CDS",
      func
    )
  
  if (nrow(genes) == 0) stop("All genes removed by crop for pair: ", prophage_id, " vs ", bacterial_id)
  
  # --- shift bacterial link coordinates to the cropped system ---
  if (!is.null(crop_win) && nrow(links_pb) > 0) {
    links_pb <- links_pb %>%
      mutate(
        start2 = as.integer(start2 - crop_win$crop_start + 1L),
        end2   = as.integer(end2   - crop_win$crop_start + 1L)
      )
  }
  
  # self links unchanged (unless you also want to crop them; usually not needed)
  links_self <- make_self_links_schema(blast_df, ids, min_align_len, min_pid)
  links <- if (show_self_links) bind_rows(links_pb, links_self) else links_pb
  
  if (is.null(title_txt)) title_txt <- prophage_id
  
  # --- seqs: enforce bacterial length to crop window length ---
  seqs <- genes %>%
    group_by(seq_id) %>%
    summarize(length = max(end, na.rm = TRUE), .groups = "drop") %>%
    mutate(length = as.integer(pmax(length, 1L)))
  
  if (!is.null(crop_win)) {
    crop_len <- as.integer(crop_win$crop_end - crop_win$crop_start + 1L)
    seqs <- seqs %>%
      mutate(length = if_else(seq_id == bacterial_id, crop_len, length))
  }
  
  seqs <- seqs %>%
    mutate(bin_id = factor(seq_id, levels = c(prophage_id, bacterial_id))) %>%
    arrange(bin_id)
  
  # --- plot ---
  p <- gggenomes(seqs = seqs, genes = genes, links = links) +
    geom_seq(linewidth = SEQ_LINEWIDTH)
  
  # p <- p +
  #   geom_seq_label(
  #     aes(label = seq_id),
  #     size = 2.8,
  #     nudge_y = 0.55
  #   )
  
  # prophage box on prophage fragment (unchanged)
  # --- provirus boxes from all_virus_summary.tsv (relative coords) ---
  box_df <- calc_virus_boxes_for_fragments(ids, virus_sum_tbl)
  
  # shift bacterial boxes into cropped coordinate system (same as genes/links)
  if (!is.null(crop_win) && nrow(box_df) > 0) {
    box_df <- box_df %>%
      mutate(
        is_bacterial = (seq_id == bacterial_id),
        box_start = if_else(is_bacterial, as.integer(box_start - crop_win$crop_start + 1L), box_start),
        box_end   = if_else(is_bacterial, as.integer(box_end   - crop_win$crop_start + 1L), box_end)
      ) %>%
      # keep only boxes that still overlap the cropped view
      filter(!(seq_id == bacterial_id) | (pmax(box_start, box_end) >= 1L)) %>%
      mutate(
        box_start = pmax(1L, box_start),
        box_end   = pmax(1L, box_end)
      ) %>%
      select(-is_bacterial)
  }
  
  if (nrow(box_df) > 0) {
    y_df <- p %>% pull_seqs() %>% distinct(seq_id, y)
    box_df <- box_df %>%
      left_join(y_df, by = "seq_id") %>%
      filter(!is.na(y)) %>%
      mutate(
        ymin = y - 0.45,
        ymax = y + 0.45
      )
    
    p <- p +
      geom_rect(
        data = box_df,
        aes(xmin = box_start, xmax = box_end, ymin = ymin, ymax = ymax),
        inherit.aes = FALSE,
        fill = box_df$box_color,   # vector -> per-row colors
        alpha = 0.18,
        color = box_df$box_color,
        linewidth = 0.7
      )
  }
  
  
  # links
  if (nrow(links) > 0) {
    if (color_links_by_pid) {
      p <- p +
        geom_link(
          aes(color = pident),
          alpha = 0.25,
          linewidth = LINK_LINEWIDTH
        ) +
        scale_color_viridis_c(
          option = "C",
          end = 0.95,
          name = "BLAST pident (%)"
        )
      
    } else {
      p <- p + geom_link(color = "grey50", alpha = 0.20, linewidth = LINK_LINEWIDTH)
    }
  }
  
  p <- p +
    geom_seq_label(
      aes(label = seq_id),
      size = 2.8,
      nudge_y = 0.55
    )
  
  # genes
  p <- p +
    geom_gene(aes(fill = func),
              position = "strand",
              size = 2.2,
              stroke = GENE_OUTLINE,
              colour = "black") +
    scale_fill_manual(
      values = func_colors,
      drop = FALSE,
      name = "PHOLD function"
    ) +
    labs(title = title_txt) +
    theme_bw(base_size = 10) +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      panel.grid = element_blank(),
      axis.title = element_blank()
    ) +
    coord_cartesian(clip = "off")
  
  p
}




# ------------------ DEBUG ONE PANEL FIRST ------------------

p1 <- make_pair_plot(
  prophage_id  = pairs$prophage_fragment_id[1],
  bacterial_id = pairs$bacterial_fragment_id[1],
  phold_df = phold,
  blast_df = blast,
  virus_sum_tbl = virus_sum,
  crop_pad = 10000
)


p1 %>% pull_seqs() %>% select(seq_id, bin_id, y)

parse_fragment_id(pairs$prophage_fragment_id[1])
parse_fragment_id(pairs$bacterial_fragment_id[1])

proph_tbl %>%
  filter(bacterial_id %in% c(
    parse_fragment_id(pairs$prophage_fragment_id[1])$base_contig,
    parse_fragment_id(pairs$bacterial_fragment_id[1])$base_contig
  )) %>%
  select(bacterial_id, prophage_start, prophage_end, genome_id) %>%
  head(5)


#box_df
#unique(box_df$seq_id)

p1 %>% track_info
nrow(p1 %>% pull_genes())   # should be > 0
head(p1 %>% pull_genes())

ggsave(DEBUG_ONE, p1, width = 16, height = 4, dpi = 300)
message("Saved debug panel: ", DEBUG_ONE)

# ------------------ BUILD PANELS ------------------

plot_order <- c(rbind(1:8, 9:16))
# plot_order <- c(rbind(1:4, 5:8, 9:12, 13:16))
pairs_ordered <- pairs %>% slice(plot_order)

plots <- purrr::pmap(
  list(pairs_ordered$prophage_fragment_id, pairs_ordered$bacterial_fragment_id),
  function(proph, bact) {
    make_pair_plot(
      prophage_id  = proph,
      bacterial_id = bact,
      phold_df = phold,
      blast_df = blast,
      virus_sum_tbl = virus_sum,
      title_txt = proph,
      min_align_len = MIN_ALIGN_LEN,
      min_pid = MIN_PID,
      show_self_links = SHOW_SELF_LINKS,
      color_links_by_pid = COLOR_LINKS_BY_PID,
      func_colors = FUNC_COLORS
    )
  }
)

legend_plot <- plots[[1]] + theme(legend.position = "right")
legend_g <- cowplot::get_legend(legend_plot)

plots_nolegend <- purrr::map(plots, ~ .x + theme(legend.position = "none"))

main_grid <- cowplot::plot_grid(plotlist = plots_nolegend, ncol = N_COL, align = "hv")

final <- cowplot::plot_grid(
  main_grid,
  legend_g,
  ncol = 2,
  rel_widths = c(1, 0.18)   # подстрой ширину легенды
)

#final <- cowplot::plot_grid(plotlist = plots, ncol = N_COL, align = "hv")

WIDTH <- 18 #18
HEIGTH <- 12 #14

ggsave(OUT_PNG, final, width = WIDTH, height = HEIGTH, dpi = 600)

ggsave(
  filename = OUT_SVG,
  plot = final,
  width = WIDTH,
  height = HEIGTH,
  device = svglite::svglite,
  bg = "white"
)

#ggsave(OUT_PDF, final, width = 25, height = 20, device = pdf)

message("Saved: ", OUT_PNG)
#message("Saved: ", OUT_PDF)
