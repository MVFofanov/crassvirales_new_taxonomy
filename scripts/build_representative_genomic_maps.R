#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
  library(gggenomes)
  library(cowplot)
  library(svglite)
})

# ============================================================
# SETTINGS
# ============================================================

WD <- "C:/crassvirales/crassvirales_new_taxonomy/crassvirales_prophages/integrated_prophages/bacterial_flank_analysis/genomic_maps_for_manuscript"

INPUT_DIR  <- file.path(WD, "data")
OUTPUT_DIR <- file.path(WD, "figures")

dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

PAIRS_TSV <- file.path(INPUT_DIR, "pairs_prophage_with_flanks_vs_bacterial_reoriented.tsv")
BLAST_TSV <- file.path(INPUT_DIR, "all_blastn_all_vs_all.tsv")
PHOLD_TSV <- file.path(INPUT_DIR, "all_phold_per_cds_predictions.tsv")
TRNA_GFF <- file.path(INPUT_DIR, "all_trnascan_out.gff")

OUT_PNG <- file.path(OUTPUT_DIR, "representative_prophage_vs_bacterial_maps.png")
OUT_SVG <- file.path(OUTPUT_DIR, "representative_prophage_vs_bacterial_maps.svg")
OUT_PDF <- file.path(OUTPUT_DIR, "representative_prophage_vs_bacterial_maps.pdf")

OUT_PROPHAGE_ONLY_PNG <- file.path(
  OUTPUT_DIR,
  "representative_prophage_only_maps.png"
)

OUT_PROPHAGE_ONLY_SVG <- file.path(
  OUTPUT_DIR,
  "representative_prophage_only_maps.svg"
)

OUT_PROPHAGE_ONLY_PDF <- file.path(
  OUTPUT_DIR,
  "representative_prophage_only_maps.pdf"
)

# Save additional prophage-only figures for these groups.
PROPHAGE_ONLY_GROUPS <- c(
  "Crassvirales",
  "Sister_branch",
  "Metacrassvirales_zeta",
  "Paracrassvirales_epsilon"
)

# Track cropping
TOP_FLANK_BP <- 50000L # 100000L
PROPHAGE_ONLY_FLANK_BP <- 10000L

# BLAST filtering
MIN_ALIGN_LEN <- 200L
MIN_PID <- 60
MAX_CLUSTER_GAP_BP <- 50000L
MAX_BACTERIAL_SPAN_BP <- 500000L
BACTERIAL_CONTEXT_BP <- 10000L

# Layout
N_COL <- 1
PANEL_WIDTH <- 30
PANEL_HEIGHT <- 2.8
DPI <- 600

# Plot appearance
SHOW_COORDINATES <- TRUE
PANEL_TEXT_SIZE <- 4.5
TITLE_TEXT_SIZE <- 16
SEQ_LINEWIDTH <- 0.35
LINK_LINEWIDTH <- 0.22
GENE_OUTLINE <- 0.15
TRNA_MARKER_SIZE <- 5
TRNA_MARKER_COLOR <- "#000000"
TRNA_MARKER_Y_OFFSET <- 0.34
PROPHAGE_COORD_Y_OFFSET <- 0.40
GRID_INTERVAL_BP <- 50000L
GRID_LINE_COLOR <- "#000000" #"#BDBDBD"
GRID_LINEWIDTH <- 0.30
GRID_LINE_ALPHA <- 0.9

# Highlight color for the prophage region on the top track
PROPHAGE_BOX_FILL <- "#ff7f00"
  PROPHAGE_BOX_ALPHA <- 0.18
  PROPHAGE_BOX_LINEWIDTH <- 0.7
  
  PRODUCT_COLORS <- c(
    "terminase large subunit" = "#ff0000", # "#8ccfb3", "#ff0000" # "#d73027",
    "major head protein"      = "#0000ff", # "#0000ff", "#e7298a",
    "portal protein"          = "#22b2cc",
    "tail sheath stabilizer"  = "#8d89b9",
    "DNA polymerase"          = "#1a9850",
    "DNA helicase"            = "#984ea3",
    "integrase"               = "#5E00D6", #"#ff99ff", # "#4575b4",
    "transposase"             = "#ff7f00",
    "excisionase and transcriptional regulator" = "#8B4513",
    "DNA transposition protein" = "#C49A00",
    "hypothetical protein"    = "#FFFFFF"
  )

  FUNCTION_COLORS <- c(
    "connector" = "#5A5A5A",
    "DNA, RNA and nucleotide metabolism" = "#f000ff",
    "head and packaging" = "#ff008d",
    "integration and excision" = "#E0B0FF",
    "lysis" = "#4deeea", #"#003F5C",
    "moron, auxiliary metabolic gene and host takeover" = "#D3D3D3", ##8900ff",
    "other" = "#D3D3D3",
    "tail" = "#74ee15",
    "transcription regulation" = "#ffe700",
    "unknown" = "#FFFFFF",
    "unknown function" = "#FFFFFF"
  )

  GENE_COLORS <- c(PRODUCT_COLORS, FUNCTION_COLORS)
  
  # ============================================================
  # HELPERS
  # ============================================================
  
  fmt_bp <- function(x) {
    format(
      as.integer(round(x)),
      big.mark = ",",
      scientific = FALSE,
      trim = TRUE
    )
  }
  
  is_yes <- function(x) {
    str_to_lower(str_trim(as.character(x))) %in% c("yes", "y", "true", "1")
  }

  make_grid_breaks <- function(x_max) {
    if (is.na(x_max) || x_max < GRID_INTERVAL_BP) {
      return(numeric(0))
    }

    seq(
      from = GRID_INTERVAL_BP,
      to = floor(x_max / GRID_INTERVAL_BP) * GRID_INTERVAL_BP,
      by = GRID_INTERVAL_BP
    )
  }
  
  assign_gene_color_group <- function(product, func) {
    product <- replace_na(as.character(product), "")
    func <- replace_na(as.character(func), "")
    
    product_lc <- str_to_lower(str_squish(product))
    func_lc <- str_to_lower(str_squish(func))

    function_labels <- set_names(
      names(FUNCTION_COLORS),
      str_to_lower(names(FUNCTION_COLORS))
    )
    function_group <- unname(function_labels[func_lc])
    
    case_when(
      str_detect(product_lc, "terminase large subunit") ~ "terminase large subunit",
      str_detect(product_lc, "\\bintegrase\\b") ~ "integrase",
      str_detect(product_lc, "dna polymerase") ~ "DNA polymerase",
      str_detect(product_lc, "dna helicase") ~ "DNA helicase",
      str_detect(product_lc, "transposase") ~ "transposase",
      str_detect(product_lc, "major head protein") ~ "major head protein",
      str_detect(product_lc, "portal protein") ~ "portal protein",
      str_detect(product_lc, "tail sheath stabilizer") ~ "tail sheath stabilizer",
      str_detect(
        product_lc,
        "excisionase and transcriptional regulator"
      ) ~ "excisionase and transcriptional regulator",
      str_detect(
        product_lc,
        "DNA transposition protein"
      ) ~ "DNA transposition protein",
      product_lc == "hypothetical protein" ~ "hypothetical protein",
      !is.na(function_group) ~ function_group,
      TRUE ~ "other"
    )
  }
  
  validate_columns <- function(df, required, object_name) {
    missing <- setdiff(required, names(df))
    if (length(missing) > 0) {
      stop(
        object_name, " is missing required columns: ",
        paste(missing, collapse = ", ")
      )
    }
  }
  
  crop_and_shift_genes <- function(phold_df, seq_id, crop_start, crop_end) {
    phold_df %>%
      filter(contig_id == seq_id) %>%
      mutate(
        gene_min = pmin(start, end),
        gene_max = pmax(start, end)
      ) %>%
      filter(gene_max >= crop_start, gene_min <= crop_end) %>%
      mutate(
        clipped_start = pmax(gene_min, crop_start),
        clipped_end = pmin(gene_max, crop_end),
        plot_start = clipped_start - crop_start + 1L,
        plot_end = clipped_end - crop_start + 1L,
        gene_fill_group = assign_gene_color_group(product, func)
      ) %>%
      transmute(
        seq_id = contig_id,
        start = as.integer(plot_start),
        end = as.integer(plot_end),
        strand = as.character(strand),
        type = "CDS",
        product = as.character(product),
        func = as.character(func),
        gene_fill_group
      )
  }

  crop_and_shift_trnas <- function(trna_df, seq_id, crop_start, crop_end) {
    trna_df %>%
      filter(
        contig_id == seq_id,
        position >= crop_start,
        position <= crop_end
      ) %>%
      transmute(
        seq_id = contig_id,
        position = as.integer(position - crop_start + 1L)
      )
  }
  
  prepare_pair_links <- function(blast_df, row_info, top_start, top_end) {
    top_id <- row_info$prophage_nucleotide_id_to_use[[1]]
    bottom_id <- row_info$bacterial_to_use[[1]]
    target_prophage_id <- row_info$prophage_id[[1]]
    
    blast_df %>%
      filter(
        prophage_id == !!target_prophage_id,
        prophage_annotation_id == !!top_id,
        bacterial_annotation_id == !!bottom_id,
        qseqid == !!top_id,
        sseqid == !!bottom_id,
        length >= MIN_ALIGN_LEN,
        pident >= MIN_PID
      ) %>%
      mutate(
        qmin = pmin(qstart, qend),
        qmax = pmax(qstart, qend),
        smin = pmin(sstart, send),
        smax = pmax(sstart, send)
      ) %>%
      filter(qmax >= top_start, qmin <= top_end)
  }

  select_bacterial_locus <- function(hits, bacterial_length) {
    if (nrow(hits) == 0) {
      return(list(
        hits = hits,
        crop_start = NA_integer_,
        crop_end = NA_integer_
      ))
    }

    hits <- hits %>%
      arrange(smin, smax) %>%
      mutate(hit_row = row_number())

    cluster_id <- integer(nrow(hits))
    current_cluster <- 0L
    current_end <- -Inf

    for (i in seq_len(nrow(hits))) {
      if (i == 1L || hits$smin[[i]] - current_end > MAX_CLUSTER_GAP_BP) {
        current_cluster <- current_cluster + 1L
        current_end <- hits$smax[[i]]
      } else {
        current_end <- max(current_end, hits$smax[[i]])
      }
      cluster_id[[i]] <- current_cluster
    }

    hits$cluster_id <- cluster_id

    cluster_summary <- hits %>%
      group_by(cluster_id) %>%
      summarise(
        aligned_bp = sum(length, na.rm = TRUE),
        total_bitscore = sum(bitscore, na.rm = TRUE),
        cluster_start = min(smin),
        cluster_end = max(smax),
        .groups = "drop"
      ) %>%
      arrange(
        desc(aligned_bp),
        desc(total_bitscore),
        cluster_start
      )

    selected_cluster <- cluster_summary$cluster_id[[1]]
    selected_hits <- hits %>%
      filter(cluster_id == selected_cluster)

    cluster_start <- min(selected_hits$smin)
    cluster_end <- max(selected_hits$smax)
    cluster_span <- cluster_end - cluster_start + 1L

    if (
      cluster_span + 2L * BACTERIAL_CONTEXT_BP <=
        MAX_BACTERIAL_SPAN_BP
    ) {
      crop_start <- max(1L, cluster_start - BACTERIAL_CONTEXT_BP)
      crop_end <- min(
        bacterial_length,
        cluster_end + BACTERIAL_CONTEXT_BP
      )
    } else {
      # Find the fixed-width window containing the greatest total number of
      # aligned bacterial bases from the selected cluster.
      candidate_starts <- unique(c(
        selected_hits$smin,
        pmax(
          1L,
          selected_hits$smax - MAX_BACTERIAL_SPAN_BP + 1L
        )
      ))

      candidate_starts <- pmin(
        candidate_starts,
        max(1L, bacterial_length - MAX_BACTERIAL_SPAN_BP + 1L)
      )

      window_scores <- map_dbl(
        candidate_starts,
        function(window_start) {
          window_end <- min(
            bacterial_length,
            window_start + MAX_BACTERIAL_SPAN_BP - 1L
          )
          sum(
            pmax(
              0,
              pmin(selected_hits$smax, window_end) -
                pmax(selected_hits$smin, window_start) + 1
            )
          )
        }
      )

      crop_start <- as.integer(
        candidate_starts[[which.max(window_scores)]]
      )
      crop_end <- min(
        bacterial_length,
        crop_start + MAX_BACTERIAL_SPAN_BP - 1L
      )
    }

    selected_hits <- selected_hits %>%
      filter(smax >= crop_start, smin <= crop_end)

    list(
      hits = selected_hits,
      crop_start = as.integer(crop_start),
      crop_end = as.integer(crop_end)
    )
  }
  
  make_coordinate_labels <- function(seq_y, crop_tbl) {
    seq_y %>%
      left_join(crop_tbl, by = "seq_id") %>%
      mutate(
        top_y = max(y),
        is_top = y == top_y,
        y_lab = if_else(is_top, y + 0.38, y - 0.38),
        vjust = if_else(is_top, 1, 0),
        left_x = 1,
        right_x = plot_length,
        left_lab = fmt_bp(crop_start),
        right_lab = fmt_bp(crop_end)
      )
  }
  
  make_one_panel <- function(
    row_info,
    phold_df,
    trna_df,
    blast_df,
    common_x_max = NULL
  ) {
    top_id <- row_info$prophage_nucleotide_id_to_use[[1]]
    bottom_id <- row_info$bacterial_to_use[[1]]
    
    top_length <- as.integer(row_info$prophage_nucleotide_length[[1]])
    bottom_length <- as.integer(row_info$bacterial_length[[1]])
    prophage_start <- as.integer(row_info$prophage_new_start[[1]])
    prophage_end <- as.integer(row_info$prophage_new_end[[1]])
    
    prophage_min <- min(prophage_start, prophage_end)
    prophage_max <- max(prophage_start, prophage_end)
    
    top_crop_start <- max(1L, prophage_min - TOP_FLANK_BP)
    top_crop_end <- min(top_length, prophage_max + TOP_FLANK_BP)
    
    all_pair_links_raw <- prepare_pair_links(
      blast_df = blast_df,
      row_info = row_info,
      top_start = top_crop_start,
      top_end = top_crop_end
    )

    selected_locus <- select_bacterial_locus(
      hits = all_pair_links_raw,
      bacterial_length = bottom_length
    )

    pair_links_raw <- selected_locus$hits
    
    if (nrow(pair_links_raw) == 0) {
      show_bottom_track <- FALSE
      warning(
        "No prophage-to-bacterial BLAST links found for representative: ",
        row_info$prophage_id[[1]],
        ". The bacterial track will be omitted."
      )
      bottom_crop_start <- NA_integer_
      bottom_crop_end <- NA_integer_
    } else {
      show_bottom_track <- TRUE
      bottom_crop_start <- selected_locus$crop_start
      bottom_crop_end <- selected_locus$crop_end

      message(
        "  Selected bacterial locus: ",
        fmt_bp(bottom_crop_start),
        "-",
        fmt_bp(bottom_crop_end),
        " (",
        nrow(pair_links_raw),
        "/",
        nrow(all_pair_links_raw),
        " qualifying BLAST hits retained)"
      )
    }
    
    top_genes <- crop_and_shift_genes(
      phold_df, top_id, top_crop_start, top_crop_end
    )

    bottom_genes <- if (show_bottom_track) {
      crop_and_shift_genes(
        phold_df, bottom_id, bottom_crop_start, bottom_crop_end
      )
    } else {
      top_genes[0, ]
    }

    genes <- bind_rows(top_genes, bottom_genes)

    top_trnas <- crop_and_shift_trnas(
      trna_df, top_id, top_crop_start, top_crop_end
    )

    bottom_trnas <- if (show_bottom_track) {
      crop_and_shift_trnas(
        trna_df, bottom_id, bottom_crop_start, bottom_crop_end
      )
    } else {
      top_trnas[0, ]
    }

    trna_markers <- bind_rows(top_trnas, bottom_trnas)
    
    if (nrow(top_genes) == 0) {
      warning("No PHOLD genes found for top track: ", top_id)
    }
    if (show_bottom_track && nrow(bottom_genes) == 0) {
      warning("No PHOLD genes found for bottom track: ", bottom_id)
    }
    
    crop_tbl <- tibble(
      seq_id = top_id,
      crop_start = top_crop_start,
      crop_end = top_crop_end
    )

    if (show_bottom_track) {
      crop_tbl <- bind_rows(
        crop_tbl,
        tibble(
          seq_id = bottom_id,
          crop_start = bottom_crop_start,
          crop_end = bottom_crop_end
        )
      )
    }

    track_levels <- if (show_bottom_track) {
      c(top_id, bottom_id)
    } else {
      top_id
    }

    crop_tbl <- crop_tbl %>%
      mutate(plot_length = crop_end - crop_start + 1L)
    
    seqs <- crop_tbl %>%
      transmute(
        seq_id,
        length = as.integer(plot_length),
        bin_id = factor(seq_id, levels = track_levels)
      ) %>%
      arrange(bin_id)
    
    links <- pair_links_raw %>%
      filter(
        smax >= bottom_crop_start,
        smin <= bottom_crop_end
      ) %>%
      transmute(
        seq_id = top_id,
        start = as.integer(pmax(qmin, top_crop_start) - top_crop_start + 1L),
        end = as.integer(pmin(qmax, top_crop_end) - top_crop_start + 1L),
        seq_id2 = bottom_id,
        start2 = as.integer(pmax(smin, bottom_crop_start) - bottom_crop_start + 1L),
        end2 = as.integer(pmin(smax, bottom_crop_end) - bottom_crop_start + 1L),
        pident
      ) %>%
      filter(start <= end, start2 <= end2)
    
    top_box <- tibble(
      seq_id = top_id,
      box_start = prophage_min - top_crop_start + 1L,
      box_end = prophage_max - top_crop_start + 1L
    )

    panel_x_max <- if (is.null(common_x_max)) {
      max(seqs$length)
    } else {
      common_x_max
    }
    
    p <- gggenomes(seqs = seqs, genes = genes, links = links) +
      geom_vline(
        xintercept = make_grid_breaks(panel_x_max),
        color = GRID_LINE_COLOR,
        linetype = "dashed",
        linewidth = GRID_LINEWIDTH,
        alpha = GRID_LINE_ALPHA
      ) +
      geom_seq(linewidth = SEQ_LINEWIDTH)
    
    seq_y <- p %>% pull_seqs() %>% distinct(seq_id, y)
    seq_y <- seq_y %>% ungroup()

    trna_markers <- trna_markers %>%
      left_join(seq_y, by = "seq_id") %>%
      mutate(
        marker_y = if_else(
          seq_id == top_id,
          y + TRNA_MARKER_Y_OFFSET,
          y - TRNA_MARKER_Y_OFFSET
        )
      )
    
    top_box <- top_box %>%
      left_join(seq_y, by = "seq_id") %>%
      mutate(
        ymin = y - 0.32,
        ymax = y + 0.32
      )

    prophage_box_coord_labels <- bind_rows(
      top_box %>%
        transmute(
          x = box_start,
          y_lab = y + PROPHAGE_COORD_Y_OFFSET,
          label = fmt_bp(prophage_min),
          hjust = 0
        ),
      top_box %>%
        transmute(
          x = box_end,
          y_lab = y + PROPHAGE_COORD_Y_OFFSET,
          label = fmt_bp(prophage_max),
          hjust = 1
        )
    )
    
    p <- p +
      geom_rect(
        data = top_box,
        aes(
          xmin = box_start,
          xmax = box_end,
          ymin = ymin,
          ymax = ymax
        ),
        inherit.aes = FALSE,
        fill = PROPHAGE_BOX_FILL,
        color = PROPHAGE_BOX_FILL,
        alpha = PROPHAGE_BOX_ALPHA,
        linewidth = PROPHAGE_BOX_LINEWIDTH
      ) +
      geom_text(
        data = prophage_box_coord_labels,
        aes(x = x, y = y_lab, label = label, hjust = hjust),
        inherit.aes = FALSE,
        vjust = 1,
        size = PANEL_TEXT_SIZE
      )
    
    if (nrow(links) > 0) {
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
    }
    
    p <- p +
      geom_gene(
        aes(fill = gene_fill_group),
        position = "strand",
        size = 10,
        stroke = GENE_OUTLINE,
        colour = "black"
      ) +
      scale_fill_manual(
        values = GENE_COLORS,
        drop = FALSE,
        name = "Gene class"
      )

    if (nrow(trna_markers) > 0) {
      p <- p +
        geom_text(
          data = trna_markers,
          aes(x = position, y = marker_y),
          label = "*",
          inherit.aes = FALSE,
          size = TRNA_MARKER_SIZE,
          color = TRNA_MARKER_COLOR,
          fontface = "bold",
          vjust = 0.5
        )
    }
    
    # Track labels
    seq_labels <- seq_y %>%
      mutate(
        top_y = max(y),
        is_top = y == top_y,
        label = case_when(
          seq_id == top_id ~ paste0(
            "Prophage-containing contig: ",
            top_id
          ),
          seq_id == bottom_id ~ paste0(
            "Bacterial contig without a prophage: ",
            bottom_id
          ),
          TRUE ~ seq_id
        ),
        # Anchor the bacterial label just inside the paired panel's lower
        # border (the visible y-range starts at 0.48).
        y_lab = if_else(is_top, y + 0.50, 0.53),
        vjust = if_else(is_top, 1, 0)
      )
    
    p <- p +
      geom_text(
        data = seq_labels,
        aes(x = 1, y = y_lab, label = label, vjust = vjust),
        inherit.aes = FALSE,
        hjust = 0,
        size = PANEL_TEXT_SIZE,
        fontface = "bold"
      )
    
    if (SHOW_COORDINATES) {
      coord_df <- make_coordinate_labels(seq_y, crop_tbl)
      
      p <- p +
        geom_text(
          data = coord_df,
          aes(x = left_x, y = y_lab, label = left_lab, vjust = vjust),
          inherit.aes = FALSE,
          hjust = 0,
          size = PANEL_TEXT_SIZE
        ) +
        geom_text(
          data = coord_df,
          aes(x = right_x, y = y_lab, label = right_lab, vjust = vjust),
          inherit.aes = FALSE,
          hjust = 1,
          size = PANEL_TEXT_SIZE
        )
    }
    
    title_txt <- paste0(
      row_info$prophage_id[[1]],
      "   |   ",
      row_info$group[[1]]
    )
    
    p +
      labs(title = title_txt, x = NULL, y = NULL) +
      theme_bw(base_size = 10) +
      theme(
        plot.title = element_text(
          size = TITLE_TEXT_SIZE,
          face = "bold",
          margin = margin(b = 5)
        ),
        panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "none",
        plot.margin = margin(4, 4, 4, 4)
      ) +
      # Keep the paired-panel height, but remove excess vertical whitespace.
      # The limits include both tracks and their names/coordinate labels.
      coord_cartesian(
        xlim = if (is.null(common_x_max)) {
          NULL
        } else {
          c(1, common_x_max)
        },
        ylim = c(0.48, 2.52),
        expand = FALSE,
        clip = "off"
      )
  }
  
  make_prophage_only_panel <- function(
    row_info,
    phold_df,
    trna_df,
    common_x_max = NULL
  ) {
    top_id <- row_info$prophage_nucleotide_id_to_use[[1]]
    
    top_length <- as.integer(
      row_info$prophage_nucleotide_length[[1]]
    )
    
    prophage_start <- as.integer(
      row_info$prophage_new_start[[1]]
    )
    
    prophage_end <- as.integer(
      row_info$prophage_new_end[[1]]
    )
    
    prophage_min <- min(prophage_start, prophage_end)
    prophage_max <- max(prophage_start, prophage_end)
    
    top_crop_start <- max(
      1L,
      prophage_min - PROPHAGE_ONLY_FLANK_BP
    )
    
    top_crop_end <- min(
      top_length,
      prophage_max + PROPHAGE_ONLY_FLANK_BP
    )
    
    plot_length <- top_crop_end - top_crop_start + 1L
    
    top_genes <- crop_and_shift_genes(
      phold_df = phold_df,
      seq_id = top_id,
      crop_start = top_crop_start,
      crop_end = top_crop_end
    )

    top_trnas <- crop_and_shift_trnas(
      trna_df = trna_df,
      seq_id = top_id,
      crop_start = top_crop_start,
      crop_end = top_crop_end
    )
    
    if (nrow(top_genes) == 0) {
      warning(
        "No PHOLD genes found for prophage-only track: ",
        top_id
      )
    }
    
    seqs <- tibble(
      seq_id = top_id,
      length = as.integer(plot_length),
      bin_id = factor(top_id, levels = top_id)
    )
    
    crop_tbl <- tibble(
      seq_id = top_id,
      crop_start = top_crop_start,
      crop_end = top_crop_end,
      plot_length = plot_length
    )
    
    top_box <- tibble(
      seq_id = top_id,
      box_start = prophage_min - top_crop_start + 1L,
      box_end = prophage_max - top_crop_start + 1L
    )

    panel_x_max <- if (is.null(common_x_max)) {
      plot_length
    } else {
      common_x_max
    }
    
    p <- gggenomes(
      seqs = seqs,
      genes = top_genes
    ) +
      geom_vline(
        xintercept = make_grid_breaks(panel_x_max),
        color = GRID_LINE_COLOR,
        linetype = "dashed",
        linewidth = GRID_LINEWIDTH,
        alpha = GRID_LINE_ALPHA
      ) +
      geom_seq(linewidth = SEQ_LINEWIDTH)
    
    seq_y <- p %>%
      pull_seqs() %>%
      distinct(seq_id, y) %>%
      ungroup()

    trna_markers <- top_trnas %>%
      left_join(seq_y, by = "seq_id") %>%
      mutate(marker_y = y + TRNA_MARKER_Y_OFFSET)
    
    top_box <- top_box %>%
      left_join(seq_y, by = "seq_id") %>%
      mutate(
        ymin = y - 0.32,
        ymax = y + 0.32
      )

    prophage_box_coord_labels <- bind_rows(
      top_box %>%
        transmute(
          x = box_start,
          y_lab = y + PROPHAGE_COORD_Y_OFFSET,
          label = fmt_bp(prophage_min),
          hjust = 0
        ),
      top_box %>%
        transmute(
          x = box_end,
          y_lab = y + PROPHAGE_COORD_Y_OFFSET,
          label = fmt_bp(prophage_max),
          hjust = 1
        )
    )
    
    p <- p +
      geom_rect(
        data = top_box,
        aes(
          xmin = box_start,
          xmax = box_end,
          ymin = ymin,
          ymax = ymax
        ),
        inherit.aes = FALSE,
        fill = PROPHAGE_BOX_FILL,
        color = PROPHAGE_BOX_FILL,
        alpha = PROPHAGE_BOX_ALPHA,
        linewidth = PROPHAGE_BOX_LINEWIDTH
      ) +
      geom_text(
        data = prophage_box_coord_labels,
        aes(x = x, y = y_lab, label = label, hjust = hjust),
        inherit.aes = FALSE,
        vjust = 1,
        size = PANEL_TEXT_SIZE
      ) +
      geom_gene(
        aes(fill = gene_fill_group),
        position = "strand",
        size = 10,
        stroke = GENE_OUTLINE,
        colour = "black"
      ) +
      scale_fill_manual(
        values = GENE_COLORS,
        drop = FALSE,
        name = "Gene class"
      )

    if (nrow(trna_markers) > 0) {
      p <- p +
        geom_text(
          data = trna_markers,
          aes(x = position, y = marker_y),
          label = "*",
          inherit.aes = FALSE,
          size = TRNA_MARKER_SIZE,
          color = TRNA_MARKER_COLOR,
          fontface = "bold",
          vjust = 0.5
        )
    }
    
    track_label <- seq_y %>%
      mutate(
        label = paste0(
          "Prophage-containing contig: ",
          top_id
        ),
        y_lab = y + 0.50
      )
    
    p <- p +
      geom_text(
        data = track_label,
        aes(
          x = 1,
          y = y_lab,
          label = label
        ),
        inherit.aes = FALSE,
        hjust = 0,
        vjust = 1,
        size = PANEL_TEXT_SIZE,
        fontface = "bold"
      )
    
    if (SHOW_COORDINATES) {
      coordinate_labels <- seq_y %>%
        left_join(crop_tbl, by = "seq_id") %>%
        mutate(
          y_lab = y - 0.38,
          left_x = 1,
          right_x = plot_length,
          left_lab = fmt_bp(crop_start),
          right_lab = fmt_bp(crop_end)
        )
      
      p <- p +
        geom_text(
          data = coordinate_labels,
          aes(
            x = left_x,
            y = y_lab,
            label = left_lab
          ),
          inherit.aes = FALSE,
          hjust = 0,
          vjust = 0,
          size = PANEL_TEXT_SIZE
        ) +
        geom_text(
          data = coordinate_labels,
          aes(
            x = right_x,
            y = y_lab,
            label = right_lab
          ),
          inherit.aes = FALSE,
          hjust = 1,
          vjust = 0,
          size = PANEL_TEXT_SIZE
        )
    }
    
    title_txt <- paste0(
      row_info$prophage_id[[1]],
      "   |   ",
      row_info$group[[1]]
    )
    
    p +
      labs(
        title = title_txt,
        x = NULL,
        y = NULL
      ) +
      theme_bw(base_size = 10) +
      theme(
        plot.title = element_text(
          size = TITLE_TEXT_SIZE,
          face = "bold",
          margin = margin(b = 5)
        ),
        panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "none",
        plot.margin = margin(4, 4, 4, 4)
      ) +
      # Keep the panel height, but use more of it for the single genomic track.
      # The limits include the upper track label and lower coordinate labels.
      coord_cartesian(
        xlim = if (is.null(common_x_max)) {
          NULL
        } else {
          c(1, common_x_max)
        },
        ylim = c(0.55, 1.52),
        expand = FALSE,
        clip = "off"
      )
  }
  
  # ============================================================
  # READ INPUTS
  # ============================================================
  
  pairs <- read_tsv(PAIRS_TSV, show_col_types = FALSE) %>%
    mutate(.input_order = row_number())
  
  required_pairs_columns <- c(
    "prophage_id",
    "prophage_nucleotide_id_to_use",
    "bacterial_to_use",
    "representative",
    "group",
    "prophage_nucleotide_length",
    "bacterial_length",
    "prophage_new_start",
    "prophage_new_end"
  )
  validate_columns(pairs, required_pairs_columns, "PAIRS_TSV")
  
  representatives <- pairs %>%
    filter(is_yes(representative)) %>%
    arrange(.input_order)
  
  if (nrow(representatives) == 0) {
    stop("No rows with representative == Yes were found in PAIRS_TSV")
  }
  
  message("Representative rows selected: ", nrow(representatives))
  
  phold <- read_tsv(PHOLD_TSV, show_col_types = FALSE) %>%
    mutate(
      contig_id = as.character(contig_id),
      start = as.integer(start),
      end = as.integer(end),
      strand = as.character(strand),
      product = as.character(product),
      func = as.character(`function`)
    )
  
  required_phold_columns <- c(
    "contig_id", "start", "end", "strand", "product", "function"
  )
  validate_columns(phold, required_phold_columns, "PHOLD_TSV")

  trnas <- read_tsv(
    TRNA_GFF,
    comment = "#",
    col_names = c(
      "contig_id",
      "source",
      "feature_type",
      "start",
      "end",
      "score",
      "strand",
      "phase",
      "attributes"
    ),
    col_types = cols(.default = col_character()),
    progress = FALSE
  ) %>%
    filter(feature_type == "tRNA") %>%
    transmute(
      contig_id = as.character(contig_id),
      position = as.integer(start)
    ) %>%
    distinct(contig_id, position)

  message("tRNA features loaded: ", nrow(trnas))
  
  required_blast_columns <- c(
    "pair_id",
    "prophage_id",
    "prophage_annotation_id",
    "bacterial_annotation_id",
    "qseqid",
    "sseqid",
    "pident",
    "length",
    "mismatch",
    "gapopen",
    "qstart",
    "qend",
    "sstart",
    "send",
    "evalue",
    "bitscore"
  )
  
  blast <- read_tsv(
    BLAST_TSV,
    show_col_types = FALSE
  )
  
  validate_columns(
    blast,
    required_blast_columns,
    "BLAST_TSV"
  )
  
  blast <- blast %>%
    mutate(
      across(
        c(
          pident,
          length,
          mismatch,
          gapopen,
          qstart,
          qend,
          sstart,
          send,
          evalue,
          bitscore
        ),
        as.numeric
      ),
      length = as.integer(length),
      qstart = as.integer(qstart),
      qend = as.integer(qend),
      sstart = as.integer(sstart),
      send = as.integer(send)
    )
  
  
  # ============================================================
  # BUILD PANELS IN ORIGINAL TABLE ORDER
  # ============================================================

  representative_prophage_common_x_max <- representatives %>%
    transmute(
      prophage_min = pmin(
        as.integer(prophage_new_start),
        as.integer(prophage_new_end)
      ),
      prophage_max = pmax(
        as.integer(prophage_new_start),
        as.integer(prophage_new_end)
      ),
      crop_start = pmax(
        1L,
        prophage_min - PROPHAGE_ONLY_FLANK_BP
      ),
      crop_end = pmin(
        as.integer(prophage_nucleotide_length),
        prophage_max + PROPHAGE_ONLY_FLANK_BP
      ),
      plot_length = crop_end - crop_start + 1L
    ) %>%
    pull(plot_length) %>%
    max(na.rm = TRUE)

  paired_common_x_max <- map_dbl(
    seq_len(nrow(representatives)),
    function(i) {
      row_info <- representatives[i, ]
      top_length <- as.integer(
        row_info$prophage_nucleotide_length[[1]]
      )
      bottom_length <- as.integer(row_info$bacterial_length[[1]])
      prophage_min <- min(
        as.integer(row_info$prophage_new_start[[1]]),
        as.integer(row_info$prophage_new_end[[1]])
      )
      prophage_max <- max(
        as.integer(row_info$prophage_new_start[[1]]),
        as.integer(row_info$prophage_new_end[[1]])
      )
      top_crop_start <- max(1L, prophage_min - TOP_FLANK_BP)
      top_crop_end <- min(
        top_length,
        prophage_max + TOP_FLANK_BP
      )
      top_plot_length <- top_crop_end - top_crop_start + 1L

      qualifying_hits <- prepare_pair_links(
        blast_df = blast,
        row_info = row_info,
        top_start = top_crop_start,
        top_end = top_crop_end
      )
      selected_locus <- select_bacterial_locus(
        hits = qualifying_hits,
        bacterial_length = bottom_length
      )

      if (nrow(selected_locus$hits) == 0) {
        return(as.numeric(top_plot_length))
      }

      bottom_plot_length <-
        selected_locus$crop_end -
        selected_locus$crop_start +
        1L

      max(top_plot_length, bottom_plot_length)
    }
  ) %>%
    max(na.rm = TRUE)

  message(
    "Representative prophage-only common horizontal scale: ",
    fmt_bp(representative_prophage_common_x_max),
    " bp"
  )
  message(
    "Paired BLAST-map common horizontal scale: ",
    fmt_bp(paired_common_x_max),
    " bp"
  )
  
  plots <- map(
    seq_len(nrow(representatives)),
    function(i) {
      message(
        "Building paired panel ", i, "/", nrow(representatives), ": ",
        representatives$prophage_id[[i]]
      )
      
      make_one_panel(
        row_info = representatives[i, ],
        phold_df = phold,
        trna_df = trnas,
        blast_df = blast,
        common_x_max = paired_common_x_max
      )
    }
  )
  
  prophage_only_plots <- map(
    seq_len(nrow(representatives)),
    function(i) {
      message(
        "Building prophage-only panel ", i, "/", nrow(representatives), ": ",
        representatives$prophage_id[[i]]
      )
      
      make_prophage_only_panel(
        row_info = representatives[i, ],
        phold_df = phold,
        trna_df = trnas,
        common_x_max = representative_prophage_common_x_max
      )
    }
  )
  
  # Build one shared legend from the first panel
  legend_source <- make_one_panel(
    row_info = representatives[1, ],
    phold_df = phold,
    trna_df = trnas,
    blast_df = blast,
    common_x_max = paired_common_x_max
  ) +
    theme(
      legend.position = "bottom",
      legend.box = "horizontal",
      legend.title = element_text(size = 14, face = "bold"),
      legend.text = element_text(size = 12),
      legend.key.size = unit(5, "mm")
    ) +
    guides(
      fill = guide_legend(
        title = "Gene class",
        nrow = 2,
        byrow = TRUE,
        order = 1
      ),
      color = guide_colorbar(
        title = "BLAST pident (%)",
        barheight = unit(5, "mm"),
        barwidth = unit(60, "mm"),
        title.position = "top",
        order = 2
      )
    )
  
  legend_g <- cowplot::get_legend(legend_source)
  
  
  prophage_only_legend_source <- make_prophage_only_panel(
    row_info = representatives[1, ],
    phold_df = phold,
    trna_df = trnas,
    common_x_max = representative_prophage_common_x_max
  ) +
    theme(
      legend.position = "bottom",
      legend.box = "horizontal",
      legend.title = element_text(size = 14, face = "bold"),
      legend.text = element_text(size = 12),
      legend.key.size = unit(5, "mm")
    ) +
    guides(
      fill = guide_legend(
        title = "Gene class",
        nrow = 2,
        byrow = TRUE
      )
    )
  
  prophage_only_legend_g <- cowplot::get_legend(
    prophage_only_legend_source
  )
  
  main_grid <- cowplot::plot_grid(
    plotlist = plots,
    ncol = N_COL,
    align = "v"
  )
  
  final_plot <- cowplot::plot_grid(
    main_grid,
    legend_g,
    ncol = 1,
    rel_heights = c(1, 0.05)
  )
  
  
  prophage_only_grid <- cowplot::plot_grid(
    plotlist = prophage_only_plots,
    ncol = N_COL,
    align = "v"
  )
  
  final_prophage_only_plot <- cowplot::plot_grid(
    prophage_only_grid,
    prophage_only_legend_g,
    ncol = 1,
    rel_heights = c(1, 0.05)
  )
  
  output_height <- PANEL_HEIGHT * ceiling(length(plots) / N_COL) + 1.8
  output_width <- PANEL_WIDTH
  
  
  PROPHAGE_ONLY_PANEL_HEIGHT <- 1.8
  prophage_only_output_height <-
    PROPHAGE_ONLY_PANEL_HEIGHT *
    ceiling(length(prophage_only_plots) / N_COL) +
    1.5
  prophage_only_output_width <- PANEL_WIDTH
  
  # ============================================================
  # SAVE OUTPUTS
  # ============================================================
  
  ggsave(
    OUT_PNG,
    final_plot,
    width = output_width,
    height = output_height,
    dpi = DPI,
    bg = "white",
    limitsize = FALSE
  )
  
  ggsave(
    OUT_SVG,
    final_plot,
    width = output_width,
    height = output_height,
    device = svglite::svglite,
    bg = "white",
    limitsize = FALSE
  )
  
  ggsave(
    OUT_PDF,
    final_plot,
    width = output_width,
    height = output_height,
    device = cairo_pdf,
    bg = "white",
    limitsize = FALSE
  )
  
  
  ggsave(
    OUT_PROPHAGE_ONLY_PNG,
    final_prophage_only_plot,
    width = prophage_only_output_width,
    height = prophage_only_output_height,
    dpi = DPI,
    bg = "white",
    limitsize = FALSE
  )
  
  ggsave(
    OUT_PROPHAGE_ONLY_SVG,
    final_prophage_only_plot,
    width = prophage_only_output_width,
    height = prophage_only_output_height,
    device = svglite::svglite,
    bg = "white",
    limitsize = FALSE
  )
  
  ggsave(
    OUT_PROPHAGE_ONLY_PDF,
    final_prophage_only_plot,
    width = prophage_only_output_width,
    height = prophage_only_output_height,
    device = cairo_pdf,
    bg = "white",
    limitsize = FALSE
  )

  # Save one prophage-only figure per requested group. Unlike the combined
  # figure above, each group figure contains all group members, irrespective of
  # the value in the representative column.
  group_output_paths <- list()

  for (group_name in PROPHAGE_ONLY_GROUPS) {
    group_rows <- pairs %>%
      filter(group == group_name) %>%
      arrange(.input_order)

    if (nrow(group_rows) == 0) {
      warning(
        "No prophages found for group: ",
        group_name,
        ". No group-specific figure was saved."
      )
      next
    }

    message(
      "Building group-specific prophage-only figure for ",
      group_name,
      " (", nrow(group_rows), " prophages)"
    )

    group_common_x_max <- group_rows %>%
      transmute(
        prophage_min = pmin(
          as.integer(prophage_new_start),
          as.integer(prophage_new_end)
        ),
        prophage_max = pmax(
          as.integer(prophage_new_start),
          as.integer(prophage_new_end)
        ),
        crop_start = pmax(
          1L,
          prophage_min - PROPHAGE_ONLY_FLANK_BP
        ),
        crop_end = pmin(
          as.integer(prophage_nucleotide_length),
          prophage_max + PROPHAGE_ONLY_FLANK_BP
        ),
        plot_length = crop_end - crop_start + 1L
      ) %>%
      pull(plot_length) %>%
      max(na.rm = TRUE)

    message(
      "  Common horizontal scale: ",
      fmt_bp(group_common_x_max),
      " bp"
    )

    group_plots <- map(
      seq_len(nrow(group_rows)),
      function(i) {
        make_prophage_only_panel(
          row_info = group_rows[i, ],
          phold_df = phold,
          trna_df = trnas,
          common_x_max = group_common_x_max
        )
      }
    )

    group_grid <- cowplot::plot_grid(
      plotlist = group_plots,
      ncol = N_COL,
      align = "v"
    )

    final_group_plot <- cowplot::plot_grid(
      group_grid,
      prophage_only_legend_g,
      ncol = 1,
      rel_heights = c(1, 0.05)
    )

    group_output_height <-
      PROPHAGE_ONLY_PANEL_HEIGHT *
      ceiling(length(group_plots) / N_COL) +
      1.5

    group_file_stem <- file.path(
      OUTPUT_DIR,
      paste0("representative_prophage_only_", group_name, "_maps")
    )

    group_paths <- c(
      png = paste0(group_file_stem, ".png"),
      svg = paste0(group_file_stem, ".svg"),
      pdf = paste0(group_file_stem, ".pdf")
    )

    ggsave(
      group_paths[["png"]],
      final_group_plot,
      width = prophage_only_output_width,
      height = group_output_height,
      dpi = DPI,
      bg = "white",
      limitsize = FALSE
    )

    ggsave(
      group_paths[["svg"]],
      final_group_plot,
      width = prophage_only_output_width,
      height = group_output_height,
      device = svglite::svglite,
      bg = "white",
      limitsize = FALSE
    )

    ggsave(
      group_paths[["pdf"]],
      final_group_plot,
      width = prophage_only_output_width,
      height = group_output_height,
      device = cairo_pdf,
      bg = "white",
      limitsize = FALSE
    )

    group_output_paths[[group_name]] <- group_paths
  }
  
  message("Saved: ", OUT_PNG)
  message("Saved: ", OUT_SVG)
  message("Saved: ", OUT_PDF)
  message("Saved: ", OUT_PROPHAGE_ONLY_PNG)
  message("Saved: ", OUT_PROPHAGE_ONLY_SVG)
  message("Saved: ", OUT_PROPHAGE_ONLY_PDF)
  walk(
    unlist(group_output_paths, use.names = FALSE),
    ~ message("Saved: ", .x)
  )
