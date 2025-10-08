#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required. Please install it with install.packages('ggplot2').", call. = FALSE)
  }
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("Package 'dplyr' is required. Please install it with install.packages('dplyr').", call. = FALSE)
  }
})

parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  inputs <- character()
  input_dir <- NULL
  output <- "pi_trend.png"
  plot_title <- NULL
  dpi <- 150
  highlights <- character()
  highlight_bed <- NULL

  i <- 1
  while (i <= length(args)) {
    arg <- args[[i]]

    if (grepl("^--input=", arg)) {
      inputs <- c(inputs, sub("^--input=", "", arg))
      i <- i + 1
      next
    }

    if (arg %in% c("--input", "-i")) {
      if (i == length(args)) stop("--input requires a value", call. = FALSE)
      inputs <- c(inputs, args[[i + 1]])
      i <- i + 2
      next
    }

    if (grepl("^--input-dir=", arg)) {
      input_dir <- sub("^--input-dir=", "", arg)
      i <- i + 1
      next
    }

    if (arg == "--input-dir") {
      if (i == length(args)) stop("--input-dir requires a value", call. = FALSE)
      input_dir <- args[[i + 1]]
      i <- i + 2
      next
    }

    if (grepl("^--output=", arg)) {
      output <- sub("^--output=", "", arg)
      i <- i + 1
      next
    }

    if (arg %in% c("--output", "-o")) {
      if (i == length(args)) stop("--output requires a value", call. = FALSE)
      output <- args[[i + 1]]
      i <- i + 2
      next
    }

    if (grepl("^--title=", arg)) {
      plot_title <- sub("^--title=", "", arg)
      i <- i + 1
      next
    }

    if (arg == "--title") {
      if (i == length(args)) stop("--title requires a value", call. = FALSE)
      plot_title <- args[[i + 1]]
      i <- i + 2
      next
    }

    if (grepl("^--dpi=", arg)) {
      dpi <- as.integer(sub("^--dpi=", "", arg))
      i <- i + 1
      next
    }

    if (arg == "--dpi") {
      if (i == length(args)) stop("--dpi requires a value", call. = FALSE)
      dpi <- as.integer(args[[i + 1]])
      i <- i + 2
      next
    }

    if (grepl("^--highlight=", arg)) {
      highlights <- c(highlights, sub("^--highlight=", "", arg))
      i <- i + 1
      next
    }

    if (arg %in% c("--highlight", "-H")) {
      if (i == length(args)) stop("--highlight requires a value", call. = FALSE)
      highlights <- c(highlights, args[[i + 1]])
      i <- i + 2
      next
    }

    if (grepl("^--highlight-bed=", arg)) {
      highlight_bed <- sub("^--highlight-bed=", "", arg)
      i <- i + 1
      next
    }

    if (arg == "--highlight-bed") {
      if (i == length(args)) stop("--highlight-bed requires a value", call. = FALSE)
      highlight_bed <- args[[i + 1]]
      i <- i + 2
      next
    }

    stop(sprintf("Unknown argument: %s", arg), call. = FALSE)
  }

  if (length(inputs) == 0 && is.null(input_dir)) {
    stop("At least one --input or --input-dir specification is required", call. = FALSE)
  }

  list(
    inputs = inputs,
    input_dir = input_dir,
    output = output,
    title = plot_title,
    dpi = dpi,
    highlights = highlights,
    highlight_bed = highlight_bed
  )
}

parse_input_spec <- function(spec) {
  parts <- strsplit(spec, "=", fixed = TRUE)[[1]]
  if (length(parts) >= 2) {
    label <- parts[1]
    path <- paste(parts[-1], collapse = "=")
    label <- trimws(label)
    path <- trimws(path)
    if (nchar(label) == 0) {
      stop(sprintf("Invalid label in input spec: %s", spec), call. = FALSE)
    }
    return(list(label = label, path = path))
  }
  list(label = NULL, path = spec)
}

parse_highlight_spec <- function(spec) {
  pattern <- "^([^:]+):(\\d+)-(\\d+)$"
  m <- regexec(pattern, spec)
  res <- regmatches(spec, m)[[1]]
  if (length(res) != 4) {
    stop(sprintf("Invalid highlight specification '%s'. Use chrom:start-end", spec), call. = FALSE)
  }
  data.frame(
    chrom = res[2],
    start = as.numeric(res[3]),
    end = as.numeric(res[4]),
    stringsAsFactors = FALSE
  )
}

read_highlight_bed <- function(path) {
  if (!file.exists(path)) {
    stop(sprintf("Highlight BED file not found: %s", path), call. = FALSE)
  }

  bed <- utils::read.table(path, header = FALSE, sep = "\t", stringsAsFactors = FALSE, comment.char = "")
  if (ncol(bed) < 3) {
    stop("Highlight BED file must have at least three columns (chrom, start, end)", call. = FALSE)
  }

  bed <- bed[, 1:3]
  colnames(bed) <- c("chrom", "start", "end")
  bed$start <- as.numeric(bed$start)
  bed$end <- as.numeric(bed$end)
  bed
}

read_pi_table <- function(path, label_override = NULL) {
  if (!file.exists(path)) {
    stop(sprintf("Input file not found: %s", path), call. = FALSE)
  }

  df <- utils::read.delim(path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

  required_cols <- c("REGION", "PICA_OUTPUT")
  missing <- setdiff(required_cols, names(df))
  if (length(missing) > 0) {
    stop(sprintf("File %s is missing required columns: %s", path, paste(missing, collapse = ", ")), call. = FALSE)
  }

  region_pattern <- "(?:[^#]+#\\d+#)?([^:]+):(\\d+)-(\\d+)"
  coords <- regmatches(df$REGION, regexec(region_pattern, df$REGION))
  if (any(vapply(coords, length, integer(1)) != 4)) {
    stop(sprintf("Failed to parse REGION values in %s", path), call. = FALSE)
  }

  chrom <- vapply(coords, function(x) x[2], character(1))
  start <- as.numeric(vapply(coords, function(x) x[3], character(1)))
  end <- as.numeric(vapply(coords, function(x) x[4], character(1)))

  df$chrom <- chrom
  df$start <- start
  df$end <- end
  df$midpoint <- (start + end) / 2

  df$pi <- as.numeric(sub(" .*$", "", df$PICA_OUTPUT))

  if ("SUBSET" %in% names(df)) {
    df$label <- as.character(df$SUBSET)
  } else if (!is.null(label_override)) {
    df$label <- label_override
  } else {
    df$label <- tools::file_path_sans_ext(basename(path))
  }
  df$source <- path
  df
}

chrom_sort_key <- function(chrom) {
  chrom <- gsub("^chr", "", chrom, ignore.case = TRUE)
  suppressWarnings(num <- as.numeric(chrom))
  ifelse(is.na(num), paste0("Z_", chrom), sprintf("%02.f", num))
}

compute_offsets <- function(df, gap = 5e5) {
  chrom_order <- df |>
    dplyr::distinct(chrom) |>
    dplyr::mutate(order_key = chrom_sort_key(chrom)) |>
    dplyr::arrange(order_key)

  chrom_max <- df |> dplyr::group_by(chrom) |> dplyr::summarise(chrom_end = max(end), .groups = "drop")

  chrom_order <- dplyr::left_join(chrom_order, chrom_max, by = "chrom")

  offsets <- numeric(nrow(chrom_order))
  cumulative <- 0
  for (i in seq_len(nrow(chrom_order))) {
    offsets[i] <- cumulative
    cumulative <- cumulative + chrom_order$chrom_end[i] + gap
  }

  chrom_order$offset <- offsets
  chrom_order
}

format_bp_value <- function(bp) {
  bp <- as.numeric(bp)
  out <- rep(NA_character_, length(bp))
  valid <- !is.na(bp)
  if (!any(valid)) {
    return(out)
  }
  out[valid] <- sprintf("%d bp", round(bp[valid]))
  mega <- valid & bp >= 1e6
  out[mega] <- sprintf("%.2f Mb", bp[mega] / 1e6)
  kilo <- valid & !mega & bp >= 1e3
  out[kilo] <- sprintf("%.1f kb", bp[kilo] / 1e3)
  out
}

format_plain_coord <- function(bp) {
  bp <- as.numeric(bp)
  fmt <- format(round(bp), scientific = FALSE, trim = TRUE)
  fmt[is.na(bp)] <- NA_character_
  fmt
}

format_mb_coord <- function(bp, digits = 2) {
  bp <- as.numeric(bp)
  fmt <- sprintf(paste0("%.", digits, "f"), bp / 1e6)
  fmt[is.na(bp)] <- NA_character_
  fmt
}

plot_pi_trend <- function(df, output, title = NULL, dpi = 150, highlights = NULL) {
  df <- df |> dplyr::arrange(factor(chrom, levels = unique(chrom)))

  chrom_offsets <- compute_offsets(df)
  df <- dplyr::left_join(df, chrom_offsets[, c("chrom", "offset")], by = "chrom")
  df$genome_pos <- df$midpoint + df$offset
  df$genome_pos_mb <- df$genome_pos / 1e6

  pop_levels <- unique(df$label)
  base_cols <- c("#d62728", "#1f77b4", "#2ca02c", "#9467bd", "#ff7f0e")
  if (length(pop_levels) <= length(base_cols)) {
    colour_values <- base_cols[seq_along(pop_levels)]
  } else {
    extra_needed <- length(pop_levels) - length(base_cols)
    extra_cols <- grDevices::hcl.colors(extra_needed, palette = "Dark3")
    colour_values <- c(base_cols, extra_cols)
  }
  names(colour_values) <- pop_levels
  df$label <- factor(df$label, levels = pop_levels)

  region_summary <- df |>
    dplyr::group_by(chrom) |>
    dplyr::summarise(
      offset = dplyr::first(offset),
      start_bp = min(start),
      end_bp = max(end),
      .groups = "drop"
    ) |>
    dplyr::arrange(offset)

  total_length_bp <- sum(region_summary$end_bp - region_summary$start_bp)
  total_length_text <- format_bp_value(total_length_bp)[1]

  window_lengths <- df$end - df$start
  window_lengths <- window_lengths[!is.na(window_lengths)]
  if (length(window_lengths) == 0) {
    window_size_text <- "unknown"
  } else {
    unique_sizes <- sort(unique(window_lengths))
    if (length(unique_sizes) == 1) {
      window_size_text <- format_bp_value(unique_sizes)[1]
    } else {
      window_size_text <- sprintf(
        "%s–%s",
        format_bp_value(min(unique_sizes))[1],
        format_bp_value(max(unique_sizes))[1]
      )
    }
  }

  subtitle_text <- sprintf(
    "Total region length: %s; Window size: %s",
    total_length_text,
    window_size_text
  )

  tick_list <- lapply(seq_len(nrow(region_summary)), function(i) {
    row <- region_summary[i, ]
    local_breaks <- pretty(c(row$start_bp, row$end_bp), n = 4)
    local_breaks <- sort(unique(c(row$start_bp, row$end_bp, local_breaks)))
    local_breaks <- local_breaks[local_breaks >= row$start_bp & local_breaks <= row$end_bp]
    data.frame(
      pos = (row$offset + local_breaks) / 1e6,
      label = format_mb_coord(local_breaks),
      stringsAsFactors = FALSE
    )
  })

  if (is.null(title)) {
    chrom_names <- unique(region_summary$chrom)
    region_start <- min(df$start, na.rm = TRUE)
    region_end <- max(df$end, na.rm = TRUE)
    if (length(chrom_names) == 1) {
      title <- sprintf(
        "PICA Output: Nucleotide Diversity Across %s:%s-%s",
        chrom_names,
        format_plain_coord(region_start),
        format_plain_coord(region_end)
      )
    } else {
      title <- "PICA Output: Nucleotide Diversity Across Multiple Regions"
    }
  }

  tick_df <- if (length(tick_list) > 0) {
    dplyr::bind_rows(tick_list) |>
      dplyr::arrange(pos) |>
      dplyr::distinct(pos, .keep_all = TRUE)
  } else {
    data.frame(pos = numeric(0), label = character(0))
  }

  chrom_axis_title <- {
    chrom_names <- unique(region_summary$chrom)
    stripped <- sub("^chr", "", chrom_names, ignore.case = TRUE)
    stripped <- sub("^chromosome ", "", stripped, ignore.case = TRUE)
    if (length(stripped) == 1) {
      sprintf("Genomic Position (Mb, chr%s)", stripped)
    } else {
      "Genomic Position (Mb)"
    }
  }

  vlines <- chrom_offsets$offset
  vlines <- vlines[vlines != 0] / 1e6

  highlight_df <- NULL
  if (!is.null(highlights) && nrow(highlights) > 0) {
    highlight_df <- highlights |>
      dplyr::left_join(chrom_offsets[, c("chrom", "offset")], by = "chrom") |>
      dplyr::filter(!is.na(offset)) |>
      dplyr::mutate(
        xmin = (start + offset) / 1e6,
        xmax = (end + offset) / 1e6,
        ymin = -Inf,
        ymax = Inf
      )
  }

  plt <- ggplot2::ggplot(df, ggplot2::aes(x = genome_pos_mb, y = pi, colour = label))

  if (!is.null(highlight_df) && nrow(highlight_df) > 0) {
    plt <- plt + ggplot2::geom_rect(
      data = highlight_df,
      inherit.aes = FALSE,
      ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
      fill = "gold",
      alpha = 0.2
    )
  }

  label_df <- NULL
  if (!is.null(highlight_df) && nrow(highlight_df) > 0) {
    ymax <- max(df$pi, na.rm = TRUE)
    label_df <- highlight_df |>
      dplyr::mutate(
        label = sprintf("%s:%.2f-%.2f Mb", chrom, start / 1e6, end / 1e6),
        x = (xmin + xmax) / 2,
        y = ymax * 1.05
      )
    plt <- plt + ggplot2::expand_limits(y = max(label_df$y, ymax, na.rm = TRUE))
  }

  plt <- plt +
    ggplot2::geom_line(linewidth = 0.9) +
    ggplot2::geom_point(size = 2.4, stroke = 0) +
    ggplot2::scale_y_log10() +
    ggplot2::annotation_logticks(sides = "l", short = grid::unit(0.15, "cm"), mid = grid::unit(0.3, "cm"), long = grid::unit(0.45, "cm")) +
    ggplot2::theme_minimal(base_size = 14) +
    ggplot2::theme(
      plot.background = ggplot2::element_rect(fill = "white", colour = NA),
      panel.background = ggplot2::element_rect(fill = "white", colour = NA),
      legend.background = ggplot2::element_rect(fill = "white", colour = "grey80"),
      legend.key = ggplot2::element_rect(fill = "white", colour = NA),
      panel.grid.major.x = ggplot2::element_line(colour = "grey90", linewidth = 0.3),
      panel.grid.minor.x = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_line(colour = "grey90", linewidth = 0.3),
      panel.grid.minor.y = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(face = "bold", size = 16, hjust = 0.5),
      plot.subtitle = ggplot2::element_text(size = 12, hjust = 0.5),
      legend.position = c(0.98, 0.98),
      legend.justification = c(1, 1)
    ) +
    ggplot2::labs(
      title = title,
      subtitle = subtitle_text,
      x = chrom_axis_title,
      y = expression(paste("Nucleotide Diversity (" * pi * ") - Log Scale")),
      colour = "Population"
    ) +
    ggplot2::scale_x_continuous(
      breaks = tick_df$pos,
      labels = tick_df$label,
      expand = ggplot2::expansion(mult = c(0.01, 0.03))
    ) +
    ggplot2::scale_colour_manual(values = colour_values) +
    ggplot2::guides(colour = ggplot2::guide_legend(title = "Population")) +
    ggplot2::coord_cartesian(clip = "off")

  if (!is.null(label_df) && nrow(label_df) > 0) {
    plt <- plt +
      ggplot2::geom_text(
        data = label_df,
        inherit.aes = FALSE,
        ggplot2::aes(x = x, y = y, label = label),
        colour = "grey20",
        size = 3,
        angle = 90,
        vjust = -0.2
      )
  }

  if (length(vlines) > 0) {
    plt <- plt + ggplot2::geom_vline(xintercept = vlines, colour = "grey80", linetype = "dashed", linewidth = 0.4)
  }
  message(sprintf("Saving plot to %s", output))
  ggplot2::ggsave(filename = output, plot = plt, width = 12, height = 6.5, dpi = dpi, units = "in")
}

main <- function() {
  opts <- parse_args()

  specs <- lapply(opts$inputs, parse_input_spec)

  if (!is.null(opts$input_dir)) {
    if (!dir.exists(opts$input_dir)) {
      stop(sprintf("Input directory not found: %s", opts$input_dir), call. = FALSE)
    }
    dir_files <- list.files(opts$input_dir, full.names = TRUE)
    if (length(dir_files) == 0) {
      stop(sprintf("No files found in input directory: %s", opts$input_dir), call. = FALSE)
    }
    file_info <- file.info(dir_files)
    dir_files <- dir_files[!is.na(file_info$isdir) & !file_info$isdir]
    if (length(dir_files) == 0) {
      stop(sprintf("No regular files found in input directory: %s", opts$input_dir), call. = FALSE)
    }
    dir_specs <- lapply(sort(dir_files), function(path) list(label = NULL, path = path))
    specs <- c(specs, dir_specs)
  }

  if (length(specs) == 0) {
    stop("No input files provided", call. = FALSE)
  }

  tables <- lapply(specs, function(spec) read_pi_table(spec$path, spec$label))
  combined <- dplyr::bind_rows(tables)

  highlight_df <- NULL
  if (length(opts$highlights) > 0) {
    highlight_df <- do.call(rbind, lapply(opts$highlights, parse_highlight_spec))
  }

  if (!is.null(opts$highlight_bed)) {
    bed_regions <- read_highlight_bed(opts$highlight_bed)
    if (is.null(highlight_df)) {
      highlight_df <- bed_regions
    } else {
      highlight_df <- rbind(highlight_df, bed_regions)
    }
  }

  plot_pi_trend(combined, opts$output, opts$title, opts$dpi, highlight_df)
}

if (identical(environment(), globalenv())) {
  suppressPackageStartupMessages({
    if (!requireNamespace("dplyr", quietly = TRUE)) {
      stop("Package 'dplyr' is required. Please install it with install.packages('dplyr').", call. = FALSE)
    }
  })

  main()
}
