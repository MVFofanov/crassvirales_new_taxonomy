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
TOP_FLANK_BP <- 100000L
BOTTOM_FLANK_BP <- 0L
PROPHAGE_ONLY_FLANK_BP <- 10000L

# BLAST filtering
MIN_ALIGN_LEN <- 200L
MIN_PID <- 0

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

# Highlight color for the prophage region on the top track
PROPHAGE_BOX_FILL <- "#ff7f00"
  PROPHAGE_BOX_ALPHA <- 0.18
  PROPHAGE_BOX_LINEWIDTH <- 0.7
  
  FUNC_COLORS <- c(
    "Terminase large subunit" = "#d73027",
    "Integrase"               = "#4575b4",
    "DNA polymerase"          = "#1a9850",
    "DNA helicase"            = "#984ea3",
    "Transposase"             = "#ff7f00",
    "Major head protein"      = "#e7298a",
    "Hypothetical protein"    = "#D3D3D3",
    "Other annotated"         = "#4D4D4D"
  )
  
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
  
  assign_gene_color_group <- function(product, func) {
    product <- replace_na(as.character(product), "")
    func <- replace_na(as.character(func), "")
    
    product_lc <- str_to_lower(product)
    func_lc <- str_to_lower(func)
    
    case_when(
      str_detect(product_lc, "terminase large subunit") ~ "Terminase large subunit",
      str_detect(product_lc, "\\bintegrase\\b") ~ "Integrase",
      str_detect(product_lc, "dna polymerase") ~ "DNA polymerase",
      str_detect(product_lc, "dna helicase") ~ "DNA helicase",
      str_detect(product_lc, "transposase") ~ "Transposase",
      str_detect(product_lc, "major head protein") ~ "Major head protein",
      product_lc == "hypothetical protein" |
        func_lc == "unknown" |
        func_lc == "unknown function" ~ "Hypothetical protein",
      TRUE ~ "Other annotated"
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
  
  make_one_panel <- function(row_info, phold_df, blast_df) {
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
    
    pair_links_raw <- prepare_pair_links(
      blast_df = blast_df,
      row_info = row_info,
      top_start = top_crop_start,
      top_end = top_crop_end
    )
    
    if (nrow(pair_links_raw) == 0) {
      warning(
        "No prophage-to-bacterial BLAST links found for representative: ",
        row_info$prophage_id[[1]],
        ". The panel will use the full bacterial track."
      )
      bottom_crop_start <- 1L
      bottom_crop_end <- bottom_length
    } else {
      bottom_crop_start <- max(1L, min(pair_links_raw$smin) - BOTTOM_FLANK_BP)
      bottom_crop_end <- min(bottom_length, max(pair_links_raw$smax) + BOTTOM_FLANK_BP)
    }
    
    top_genes <- crop_and_shift_genes(
      phold_df, top_id, top_crop_start, top_crop_end
    )
    bottom_genes <- crop_and_shift_genes(
      phold_df, bottom_id, bottom_crop_start, bottom_crop_end
    )
    genes <- bind_rows(top_genes, bottom_genes)
    
    if (nrow(top_genes) == 0) {
      warning("No PHOLD genes found for top track: ", top_id)
    }
    if (nrow(bottom_genes) == 0) {
      warning("No PHOLD genes found for bottom track: ", bottom_id)
    }
    
    crop_tbl <- tibble(
      seq_id = c(top_id, bottom_id),
      crop_start = c(top_crop_start, bottom_crop_start),
      crop_end = c(top_crop_end, bottom_crop_end)
    ) %>%
      mutate(plot_length = crop_end - crop_start + 1L)
    
    seqs <- crop_tbl %>%
      transmute(
        seq_id,
        length = as.integer(plot_length),
        bin_id = factor(seq_id, levels = c(top_id, bottom_id))
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
    
    p <- gggenomes(seqs = seqs, genes = genes, links = links) +
      geom_seq(linewidth = SEQ_LINEWIDTH)
    
    seq_y <- p %>% pull_seqs() %>% distinct(seq_id, y)
    
    top_box <- top_box %>%
      left_join(seq_y, by = "seq_id") %>%
      mutate(
        ymin = y - 0.32,
        ymax = y + 0.32
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
        values = FUNC_COLORS,
        drop = FALSE,
        name = "Gene class"
      )
    
    # Track labels
    seq_labels <- seq_y %>%
      mutate(
        top_y = max(y),
        is_top = y == top_y,
        label = if_else(
          is_top,
          paste0("Prophage-containing contig: ", top_id),
          paste0("Related bacterial contig: ", bottom_id)
        ),
        y_lab = if_else(is_top, y + 0.50, y - 0.50),
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
      coord_cartesian(clip = "off")
  }
  
  make_prophage_only_panel <- function(row_info, phold_df) {
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
    
    p <- gggenomes(
      seqs = seqs,
      genes = top_genes
    ) +
      geom_seq(linewidth = SEQ_LINEWIDTH)
    
    seq_y <- p %>%
      pull_seqs() %>%
      distinct(seq_id, y)
    
    top_box <- top_box %>%
      left_join(seq_y, by = "seq_id") %>%
      mutate(
        ymin = y - 0.32,
        ymax = y + 0.32
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
      geom_gene(
        aes(fill = gene_fill_group),
        position = "strand",
        size = 10,
        stroke = GENE_OUTLINE,
        colour = "black"
      ) +
      scale_fill_manual(
        values = FUNC_COLORS,
        drop = FALSE,
        name = "Gene class"
      )
    
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
        blast_df = blast
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
        phold_df = phold
      )
    }
  )
  
  # Build one shared legend from the first panel
  legend_source <- make_one_panel(
    row_info = representatives[1, ],
    phold_df = phold,
    blast_df = blast
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
    phold_df = phold
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

    group_plots <- map(
      seq_len(nrow(group_rows)),
      function(i) {
        make_prophage_only_panel(
          row_info = group_rows[i, ],
          phold_df = phold
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
