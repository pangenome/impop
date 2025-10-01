#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required. Install it with install.packages('ggplot2').", call. = FALSE)
  }
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("Package 'dplyr' is required. Install it with install.packages('dplyr').", call. = FALSE)
  }
})

parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  inputs <- character()
  output <- "fst_trend.png"
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

  if (length(inputs) == 0) {
    stop("At least one --input specification is required", call. = FALSE)
  }

  list(
    inputs = inputs,
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
    label <- trimws(parts[1])
    path <- trimws(paste(parts[-1], collapse = "="))
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

read_fst_table <- function(path, label_override = NULL) {
  if (!file.exists(path)) {
    stop(sprintf("Input file not found: %s", path), call. = FALSE)
  }

  df <- utils::read.delim(path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

  required_cols <- c("REGION", "FST")
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

  fst_vals <- df$FST
  fst_vals[fst_vals %in% c("NA", "")] <- NA
  df$fst <- suppressWarnings(as.numeric(fst_vals))

  if (!is.null(label_override)) {
    label <- label_override
  } else {
    label <- tools::file_path_sans_ext(basename(path))
  }

  df$label <- label
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

  chrom_max <- df |>
    dplyr::group_by(chrom) |>
    dplyr::summarise(chrom_end = max(end), .groups = "drop")

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

plot_fst_trend <- function(df, output, title = NULL, dpi = 150, highlights = NULL) {
  df <- dplyr::arrange(df, factor(chrom, levels = unique(chrom)))

  chrom_offsets <- compute_offsets(df)
  df <- dplyr::left_join(df, chrom_offsets[, c("chrom", "offset")], by = "chrom")
  df$genome_pos <- df$midpoint + df$offset

  comparison_levels <- unique(df$label)
  base_cols <- c("#ff595e", "#ffca3a", "#8ac926", "#1982c4", "#6a4c93")
  if (length(comparison_levels) <= length(base_cols)) {
    colour_values <- base_cols[seq_along(comparison_levels)]
  } else {
    extra_needed <- length(comparison_levels) - length(base_cols)
    extra_cols <- grDevices::hcl.colors(extra_needed, palette = "Dark3")
    colour_values <- c(base_cols, extra_cols)
  }
  names(colour_values) <- comparison_levels
  df$label <- factor(df$label, levels = comparison_levels)

  centers <- df |>
    dplyr::group_by(chrom) |>
    dplyr::summarise(center = mean(genome_pos), .groups = "drop") |>
    dplyr::mutate(order_key = chrom_sort_key(chrom)) |>
    dplyr::arrange(order_key) |>
    dplyr::left_join(chrom_offsets[, c("chrom", "offset")], by = "chrom") |>
    dplyr::mutate(label = sprintf("%s\n%.1f Mb", chrom, (center - offset) / 1e6))

  vlines <- chrom_offsets$offset
  vlines <- vlines[vlines != 0]

  highlight_df <- NULL
  if (!is.null(highlights) && nrow(highlights) > 0) {
    highlight_df <- highlights |>
      dplyr::left_join(chrom_offsets[, c("chrom", "offset")], by = "chrom") |>
      dplyr::filter(!is.na(offset)) |>
      dplyr::mutate(
        xmin = start + offset,
        xmax = end + offset,
        ymin = -Inf,
        ymax = Inf
      )
  }

  plt <- ggplot2::ggplot(df, ggplot2::aes(x = genome_pos, y = fst, colour = label))

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
    ymax <- max(df$fst, na.rm = TRUE)
    label_df <- highlight_df |>
      dplyr::mutate(
        label = sprintf("%s:%.2f-%.2f Mb", chrom, start / 1e6, end / 1e6),
        x = (xmin + xmax) / 2,
        y = ymax * 1.05
      )
    plt <- plt + ggplot2::expand_limits(y = max(label_df$y, ymax, na.rm = TRUE))
  }

  plt <- plt +
    ggplot2::geom_line(linewidth = 0.6) +
    ggplot2::geom_point(size = 1.2) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      plot.background = ggplot2::element_rect(fill = "white", colour = NA),
      panel.background = ggplot2::element_rect(fill = "white", colour = NA),
      legend.background = ggplot2::element_rect(fill = "white", colour = NA),
      legend.key = ggplot2::element_rect(fill = "white", colour = NA)
    ) +
    ggplot2::labs(x = "Genomic position (window midpoint)", y = expression(F[ST]), colour = "Comparison") +
    ggplot2::scale_x_continuous(breaks = centers$center, labels = centers$label) +
    ggplot2::scale_colour_manual(values = colour_values) +
    ggplot2::guides(colour = ggplot2::guide_legend(title = "Comparison")) +
    ggplot2::theme(panel.grid.major.x = ggplot2::element_blank(), panel.grid.minor.x = ggplot2::element_blank())

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
    plt <- plt + ggplot2::geom_vline(xintercept = vlines, colour = "grey70", linetype = "dashed", linewidth = 0.3)
  }

  if (!is.null(title)) {
    plt <- plt + ggplot2::ggtitle(title)
  }

  ggplot2::ggsave(output, plot = plt, dpi = dpi, width = 10, height = 5, units = "in")
}

main <- function() {
  args <- parse_args()

  highlight_regions <- NULL
  if (length(args$highlights) > 0) {
    highlight_regions <- do.call(rbind, lapply(args$highlights, parse_highlight_spec))
  }

  if (!is.null(args$highlight_bed)) {
    bed <- read_highlight_bed(args$highlight_bed)
    highlight_regions <- if (is.null(highlight_regions)) bed else rbind(highlight_regions, bed)
  }

  tables <- lapply(args$inputs, function(spec) {
    parsed <- parse_input_spec(spec)
    read_fst_table(parsed$path, label_override = parsed$label)
  })

  combo <- dplyr::bind_rows(tables)
  combo <- combo[!is.na(combo$fst), , drop = FALSE]

  if (nrow(combo) == 0) {
    stop("No valid Fst values found after filtering.", call. = FALSE)
  }

  plot_fst_trend(combo, args$output, title = args$title, dpi = args$dpi, highlights = highlight_regions)
}

main()
