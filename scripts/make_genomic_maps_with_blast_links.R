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
TRNA_GFF <- file.path(WD, "all_trnascan_out.gff")
TAX_TSV <- file.path(WD, "all_contigs_virus_summary_ge_10kb_sorted_by_seq_name_ids_uniq_taxonomy.txt")

OUT_PNG <- file.path(WD, "prophage_vs_bacterial_2x8.png")
OUT_SVG <- file.path(WD, "prophage_vs_bacterial_2x8.svg")
OUT_PDF <- file.path(WD, "prophage_vs_bacterial_2x8.pdf")

DEBUG_ONE <- file.path(WD, "debug_one_panel.png")

N_COL <- 2 # 2

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

GENE_LABEL_LEGEND <- tibble::tribble(
  ~meaning,                      ~sym,
  "Terminase large subunit",   "t",
  "Integrase",                   "i",
  "transposase",                 "+",
  "DNA polymerase",              "p",
  "DNA helicase",                "h",
  "tRNA",                        "*"
) %>%
  mutate(
    shape = case_when(
      sym == "*" ~ 8,   # star
      TRUE ~ 3          # plus; we'll print letters via override (below)
    )
  )

gene_legend_grob <- cowplot::ggdraw() +
  cowplot::draw_label(
    "t = Terminase large subunit    i = Integrase   + = Transposase    p = DNA polymerase    h = DNA helicase    * = tRNA",
    x = 0, hjust = 0, size = 16
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

trna <- read_tsv(
  TRNA_GFF,
  comment = "#",
  col_names = c("seq_id","source","type","start","end","score","strand","phase","attributes"),
  show_col_types = FALSE
) %>%
  filter(type == "tRNA") %>%
  transmute(
    seq_id = as.character(seq_id),
    start = as.integer(start),
    end   = as.integer(end),
    strand = as.character(strand),
    x = (start + end) / 2
  )

tax_raw <- read_tsv(TAX_TSV, show_col_types = FALSE) %>%
  transmute(
    accession = str_trim(as.character(accession)),
    taxonomy  = str_trim(as.character(taxonomy))
  )

tax_tbl <- read_tsv(TAX_TSV, show_col_types = FALSE) %>%
  transmute(
    accession = as.character(accession),
    taxonomy  = as.character(taxonomy)
  ) %>%
  mutate(
    # split by ";" and trim whitespace
    tax_list = str_split(taxonomy, "\\s*;\\s*"),
    class = map_chr(tax_list, ~ if (length(.x) >= 4)  .x[[4]] else NA_character_),
    family= map_chr(tax_list, ~ if (length(.x) >= 7)  .x[[7]] else NA_character_),
    genus = map_chr(tax_list, ~ if (length(.x) >= 8)  .x[[8]] else NA_character_),
    tax_short = str_c(class, family, genus, sep = "; ")
  ) %>%
  select(-tax_list)


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

extract_gene_labels <- function(genes_df) {
  genes_df %>%
    mutate(
      product_lc = tolower(product),
      label = case_when(
        # existing
        str_detect(product_lc, "terminase large subunit") ~ "t",
        str_detect(product_lc, "\\bintegrase\\b") ~ "i",
        str_detect(product_lc, "transposase") ~ "+",
        
        # NEW
        str_detect(product_lc, "dna polymerase") ~ "p",
        str_detect(product_lc, "dna helicase")  ~ "h",
        str_detect(product_lc, "major head protein")  ~ "m",
        
        TRUE ~ NA_character_
      ),
      x = (start + end) / 2
    ) %>%
    filter(!is.na(label)) %>%
    select(seq_id, start, end, strand, product, label, x)
}


crop_shift_trna <- function(trna_df, prophage_id, bacterial_id, crop_win) {
  if (is.null(trna_df) || nrow(trna_df) == 0) return(trna_df)
  
  out <- trna_df %>%
    filter(seq_id %in% c(prophage_id, bacterial_id)) %>%
    mutate(is_bacterial = (seq_id == bacterial_id))
  
  if (!is.null(crop_win)) {
    out <- out %>%
      # keep only tRNAs that overlap crop on bacterial
      filter(!is_bacterial | (end >= crop_win$crop_start & start <= crop_win$crop_end)) %>%
      mutate(
        start = if_else(is_bacterial, start - crop_win$crop_start + 1L, start),
        end   = if_else(is_bacterial, end   - crop_win$crop_start + 1L, end),
        x     = if_else(is_bacterial, (start + end) / 2, x)
      )
  }
  
  out
}

seqid_to_accession <- function(x) {
  # "NZ_OZ245719.1_3456868-3729269" -> "NZ_OZ245719.1"
  str_replace(x, "_\\d+-\\d+$", "")
}

make_tax_label_df <- function(p, seqs, tax_tbl, off = 2, x_pad = 2000) {
  # seq positions (y) from gggenomes
  seq_y <- p %>% pull_seqs() %>% distinct(seq_id, y)
  top_y <- max(seq_y$y)  # gggenomes often has reversed y
  
  # right edge per seq from seqs table
  # (seqs has length already; for bacterial it’s crop_len if you set it)
  seq_right <- seqs %>%
    select(seq_id, length) %>%
    mutate(x_right = as.numeric(length))
  
  # join taxonomy
  out <- seq_y %>%
    left_join(seq_right, by = "seq_id") %>%
    mutate(
      accession = seqid_to_accession(seq_id)
    ) %>%
    left_join(tax_tbl, by = "accession") %>%
    mutate(
      is_top = (y == top_y),
      y_lab  = if_else(is_top, y - off, y + off),
      vjust  = if_else(is_top, 1, 0),
      x_lab  = x_right + x_pad,
      label  = if_else(is.na(tax_short) | tax_short == "NA; NA; NA",
                       accession,   # fallback if taxonomy missing
                       tax_short)
    )
  
  out
}

debug_taxonomy_matches <- function(seq_ids, tax_raw) {
  tibble(seq_id = unique(as.character(seq_ids))) %>%
    mutate(
      accession = seqid_to_accession(seq_id),
      acc_plus_NZ  = if_else(str_starts(accession, "NZ_"), accession, paste0("NZ_", accession)),
      acc_minus_NZ = str_replace(accession, "^NZ_", "")
    ) %>%
    left_join(tax_raw %>% rename(tax_exact = taxonomy), by = c("accession" = "accession")) %>%
    left_join(tax_raw %>% rename(tax_plus_NZ = taxonomy), by = c("acc_plus_NZ" = "accession")) %>%
    left_join(tax_raw %>% rename(tax_minus_NZ = taxonomy), by = c("acc_minus_NZ" = "accession")) %>%
    mutate(
      match_type = case_when(
        !is.na(tax_exact) ~ "exact",
        !is.na(tax_plus_NZ) ~ "matched_by_adding_NZ_",
        !is.na(tax_minus_NZ) ~ "matched_by_removing_NZ_",
        TRUE ~ "NO_MATCH"
      ),
      taxonomy_found = coalesce(tax_exact, tax_plus_NZ, tax_minus_NZ)
    ) %>%
    select(seq_id, accession, match_type, taxonomy_found, tax_exact, tax_plus_NZ, tax_minus_NZ)
}

tax_last4 <- function(tax) {
  if (is.na(tax) || tax == "") return(NA_character_)
  parts <- str_split(tax, "\\s*;\\s*")[[1]]
  parts <- parts[parts != ""]
  if (length(parts) == 0) return(NA_character_)
  paste(tail(parts, 4), collapse = "; ")
}


make_pair_plot <- function(prophage_id, bacterial_id,
                           phold_df, blast_df, virus_sum_tbl, trna_df, tax_tbl,
                           prophage_taxonomy = NA_character_,
                           bacterial_taxonomy = NA_character_,
                           prophage_is_mag = NA_character_,
                           bacterial_is_mag = NA_character_,
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
  
  trna_pair <- crop_shift_trna(trna_df, prophage_id, bacterial_id, crop_win)
  
  # --- genes (PHOLD) ---
  genes_raw <- phold_df %>%
    filter(contig_id %in% ids) %>%
    filter(!is.na(gene_min), !is.na(gene_max)) %>%
    transmute(
      seq_id = as.character(contig_id),
      gene_min = as.integer(gene_min),
      gene_max = as.integer(gene_max),
      strand = as.character(strand),
      func = func,
      product = as.character(product)   # <-- add this
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
      func,
      product   # <-- add this
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
  
  # --- sequence labels: top above, bottom below (works with reversed y) ---
  seq_y <- p %>% pull_seqs() %>% distinct(seq_id, y)
  
  top_y <- max(seq_y$y)   # in gggenomes, this is often the TOP track visually
  off   <- 2 #0.65
  
  label_df <- seq_y %>%
    mutate(
      label = seq_id,
      is_top = (y == top_y),
      y_lab  = if_else(is_top, y + off, y - off),   # NOTE: swapped signs
      vjust  = if_else(is_top, 1, 0)                # above: anchor bottom; below: anchor top
    )
  
  p <- p +
    geom_text(
      data = label_df,
      aes(x = 1, y = y_lab, label = label, vjust = vjust),
      inherit.aes = FALSE,
      size = 2.8,
      hjust = 0
    )
  
  # --- right-side taxonomy labels INSIDE panel (aligned + right-justified) ---
  seq_y <- p %>% pull_seqs() %>% distinct(seq_id, y)
  top_y <- max(seq_y$y)
  off   <- 2
  
  x_max <- max(
    seqs$length,
    genes$end,
    if (nrow(links) > 0) max(c(links$end, links$end2), na.rm = TRUE) else 0,
    na.rm = TRUE
  )
  
  x_anchor <- x_max * 0.98   # <-- inside the panel (tune 0.95–0.99)
  
  tax_df <- seq_y %>%
    mutate(
      tax = case_when(
        seq_id == prophage_id  ~ tax_last4(prophage_taxonomy),
        seq_id == bacterial_id ~ tax_last4(bacterial_taxonomy),
        TRUE ~ NA_character_
      ),
      is_top = (y == top_y),
      y_lab  = if_else(is_top, y + off, y - off),
      vjust  = if_else(is_top, 1, 0),
      x_lab  = x_anchor
    ) %>%
    filter(!is.na(tax), tax != "")
  
  p <- p +
    geom_text(
      data = tax_df,
      aes(x = x_lab, y = y_lab, label = tax, vjust = vjust),
      inherit.aes = FALSE,
      size = 2.8, #2.8
      hjust = 1   # <-- right-align so it extends LEFT, not outside
    )
  
  # --- mid-panel MAG/isolate labels (centered) ---
  seq_y <- p %>% pull_seqs() %>% distinct(seq_id, y)
  top_y <- max(seq_y$y)
  off   <- 2
  
  # use same x_max you already compute for taxonomy
  x_max <- max(
    seqs$length,
    genes$end,
    if (nrow(links) > 0) max(c(links$end, links$end2), na.rm = TRUE) else 0,
    na.rm = TRUE
  )
  
  x_mid <- x_max * 0.50  # center of the panel (tune if you want)
  
  mag_df <- seq_y %>%
    mutate(
      mag = case_when(
        seq_id == prophage_id  ~ as.character(prophage_is_mag),
        seq_id == bacterial_id ~ as.character(bacterial_is_mag),
        TRUE ~ NA_character_
      ),
      mag = if_else(is.na(mag) | mag == "", NA_character_, mag),
      is_top = (y == top_y),
      y_lab  = if_else(is_top, y + off, y - off),
      vjust  = if_else(is_top, 1, 0),
      x_lab  = x_mid
    ) %>%
    filter(!is.na(mag))
  
  if (nrow(mag_df) > 0) {
    p <- p +
      geom_text(
        data = mag_df,
        aes(x = x_lab, y = y_lab, label = mag, vjust = vjust),
        inherit.aes = FALSE,
        size = 2.6,
        hjust = 0.5,       # centered text
        fontface = "bold"  # optional, makes it easier to see
      )
  }
  
  
  
  
  
  
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
    labs(title = title_txt,
         x = "Position (bp)") +
    theme_bw(base_size = 10) +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      panel.grid = element_blank(),
      axis.title = element_blank(),
      axis.text.x  = element_text(size = 14),   # tick labels
      axis.title.x = element_text(size = 14),    # if you ever show the title
      axis.ticks.x = element_line(linewidth = 0.6)
    ) +
    coord_cartesian(clip = "off")
  
  # # --- Gene label legend (no colors; just a key with letters) ---
  # GENE_LABEL_LEGEND2 <- GENE_LABEL_LEGEND %>%
  #   mutate(label = paste0(sym, "  =  ", meaning))
  # 
  # p <- p +
  #   geom_point(
  #     data = GENE_LABEL_LEGEND2,
  #     aes(x = Inf, y = Inf, shape = label),
  #     inherit.aes = FALSE,
  #     alpha = 0,
  #     show.legend = TRUE
  #   ) +
  #   scale_shape_discrete(name = "Gene / feature labels") +
  #   guides(shape = guide_legend(override.aes = list(alpha = 1, size = 4)))
  # 
  
  gene_labels <- extract_gene_labels(genes)
  
  if (nrow(gene_labels) > 0) {
    y_df <- p %>% pull_seqs() %>% distinct(seq_id, y)
    
    gene_labels <- gene_labels %>%
      left_join(y_df, by = "seq_id") %>%
      filter(!is.na(y))
    
    p <- p +
      geom_text(
        data = gene_labels,
        aes(x = x, y = y, label = label),
        inherit.aes = FALSE,
        size = 3.2,
        fontface = "bold",
        vjust = -0.9
      )
  }
  
  # --- tRNA labels as "*" ---
  if (!is.null(trna_pair) && nrow(trna_pair) > 0) {
    y_df <- p %>% pull_seqs() %>% distinct(seq_id, y)
    
    trna_lab <- trna_pair %>%
      left_join(y_df, by = "seq_id") %>%
      filter(!is.na(y)) %>%
      mutate(y_star = y + 0.35)  # tweak vertical placement
    
    p <- p +
      geom_text(
        data = trna_lab,
        aes(x = x, y = y_star),
        label = "*",
        inherit.aes = FALSE,
        size = 4.2,
        fontface = "bold"
      )
  }
  
  
  p
}




# ------------------ DEBUG ONE PANEL FIRST ------------------

p1 <- make_pair_plot(
  prophage_id  = pairs$prophage_fragment_id[10],
  bacterial_id = pairs$bacterial_fragment_id[10],
  phold_df = phold,
  blast_df = blast,
  virus_sum_tbl = virus_sum,
  trna_df = trna,
  tax_tbl = tax_tbl,
  crop_pad = 10000
)

seq_ids_p1 <- p1 %>% pull_seqs() %>% pull(seq_id)
debug_taxonomy_matches(seq_ids_p1, tax_raw) %>% print(n = 50)

seq_ids_all <- unique(c(pairs$prophage_fragment_id, pairs$bacterial_fragment_id))
debug_taxonomy_matches(seq_ids_all, tax_raw) %>% print(n = 200)

p1 %>% pull_seqs() %>% select(seq_id, bin_id, y)

parse_fragment_id(pairs$prophage_fragment_id[10])
parse_fragment_id(pairs$bacterial_fragment_id[10])

proph_tbl %>%
  filter(bacterial_id %in% c(
    parse_fragment_id(pairs$prophage_fragment_id[10])$base_contig,
    parse_fragment_id(pairs$bacterial_fragment_id[10])$base_contig
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
  list(
    pairs_ordered$prophage_fragment_id,
    pairs_ordered$bacterial_fragment_id,
    pairs_ordered$prophage_taxonomy,
    pairs_ordered$bacterial_taxonomy,
    pairs_ordered$prophage_is_mag,
    pairs_ordered$bacterial_is_mag
  ),
  function(proph, bact, proph_tax, bact_tax, proph_mag, bact_mag) {
    make_pair_plot(
      prophage_id  = proph,
      bacterial_id = bact,
      phold_df = phold,
      blast_df = blast,
      virus_sum_tbl = virus_sum,
      trna_df = trna,
      tax_tbl = tax_tbl,
      prophage_taxonomy = proph_tax,
      bacterial_taxonomy = bact_tax,
      prophage_is_mag = proph_mag,
      bacterial_is_mag = bact_mag,
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
# --- get ONE legend from the first plot ---
legend_g <- cowplot::get_legend(
  plots[[1]] +
    theme(
      legend.position = "bottom",
      legend.box = "vertical",
      
      # --- text sizes ---
      legend.title = element_text(size = 18, face = "bold"),
      legend.text  = element_text(size = 16),
      
      # --- make the colored boxes bigger ---
      legend.key.size = unit(6, "mm"),
      
      # optional: more spacing
      legend.spacing.x = unit(4, "mm"),
      legend.spacing.y = unit(2, "mm")
    ) +
    guides(
      # optional: bigger colorbar
      color = guide_colorbar(
        barheight = unit(5, "mm"),
        barwidth  = unit(60, "mm"),
        title.position = "top"
      ),
      # optional: arrange PHOLD categories into multiple rows
      fill = guide_legend(nrow = 2, byrow = TRUE)
    )
)

# --- remove legends from panels ---
plots_nolegend <- purrr::map(plots, ~ .x + theme(legend.position = "none"))

# --- main grid of panels ---
main_grid <- cowplot::plot_grid(
  plotlist = plots_nolegend,
  ncol = N_COL,
  align = "hv"
)

# --- combine: grid + legend at bottom ---
final <- cowplot::plot_grid(
  main_grid,
  legend_g,
  gene_legend_grob,
  ncol = 1,
  rel_heights = c(1, 0.12, 0.05)   # tune legend height (e.g. 0.10–0.20)
)


WIDTH <- 18 #18
HEIGTH <- 20 #14

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
