#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required. Please install it with install.packages('ggplot2').", call. = FALSE)
  }
})

parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  inputs <- character()
  output <- "pi_trend.png"
  plot_title <- NULL
  dpi <- 150

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

    stop(sprintf("Unknown argument: %s", arg), call. = FALSE)
  }

  if (length(inputs) == 0) {
    stop("At least one --input specification is required", call. = FALSE)
  }

  list(inputs = inputs, output = output, title = plot_title, dpi = dpi)
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

  if (!is.null(label_override)) {
    label <- label_override
  } else if ("SUBSET" %in% names(df) && length(unique(df$SUBSET)) == 1) {
    label <- unique(df$SUBSET)
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
  chrom_order <- df |> dplyr::distinct(chrom) |> dplyr::mutate(order_key = chrom_sort_key(chrom)) |
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

plot_pi_trend <- function(df, output, title = NULL, dpi = 150) {
  df <- df |> dplyr::arrange(factor(chrom, levels = unique(chrom)))

  chrom_offsets <- compute_offsets(df)
  df <- dplyr::left_join(df, chrom_offsets[, c("chrom", "offset")], by = "chrom")
  df$genome_pos <- df$midpoint + df$offset

  centers <- df |> dplyr::group_by(chrom) |> dplyr::summarise(center = mean(genome_pos), .groups = "drop") |
    dplyr::arrange(chrom_sort_key(chrom))

  vlines <- chrom_offsets$offset
  vlines <- vlines[vlines != 0]

  plt <- ggplot2::ggplot(df, ggplot2::aes(x = genome_pos, y = pi, colour = label)) +
    ggplot2::geom_line(linewidth = 0.6) +
    ggplot2::geom_point(size = 1.2) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::labs(x = "Genomic position (window midpoint)", y = expression(pi), colour = "Population") +
    ggplot2::scale_x_continuous(breaks = centers$center, labels = centers$chrom) +
    ggplot2::guides(colour = ggplot2::guide_legend(title = "Population")) +
    ggplot2::theme(panel.grid.major.x = ggplot2::element_blank(), panel.grid.minor.x = ggplot2::element_blank())

  if (length(vlines) > 0) {
    plt <- plt + ggplot2::geom_vline(xintercept = vlines, colour = "grey70", linetype = "dashed", linewidth = 0.3)
  }

  if (!is.null(title)) {
    plt <- plt + ggplot2::ggtitle(title)
  }

  message(sprintf("Saving plot to %s", output))
  ggplot2::ggsave(filename = output, plot = plt, width = 11, height = 4, dpi = dpi, units = "in")
}

main <- function() {
  opts <- parse_args()

  specs <- lapply(opts$inputs, parse_input_spec)
  tables <- lapply(specs, function(spec) read_pi_table(spec$path, spec$label))
  combined <- dplyr::bind_rows(tables)

  plot_pi_trend(combined, opts$output, opts$title, opts$dpi)
}

if (identical(environment(), globalenv())) {
  suppressPackageStartupMessages({
    if (!requireNamespace("dplyr", quietly = TRUE)) {
      stop("Package 'dplyr' is required. Please install it with install.packages('dplyr').", call. = FALSE)
    }
  })

  main()
}
