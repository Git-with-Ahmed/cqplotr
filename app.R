suppressPackageStartupMessages({
  library(shiny)
  library(ggplot2)
  library(tibble)
  library(munsell)
  library(DT)
  library(rhandsontable)
})

example_cq_data <- function() {
  data.frame(
    Well = c(
      "A01", "A02", "A03",
      "B01", "B02", "B03",
      "C01", "C02", "C03",
      "D01", "D02", "D03"
    ),
    Cq = c(
      18.2, 18.5, 18.1,
      18.7, 18.9, 18.6,
      22.4, 22.1, 22.7,
      20.8, 20.5, 20.9
    ),
    stringsAsFactors = FALSE
  )
}

example_plate_data <- function() {
  plate <- create_plate_template()
  plate[plate$Row == "A" & plate$Type == "Sample", c("1", "2", "3")] <- "Cnt"
  plate[plate$Row == "A" & plate$Type == "Gene", c("1", "2", "3")] <- "GAPDH"
  plate[plate$Row == "B" & plate$Type == "Sample", c("1", "2", "3")] <- "LPS"
  plate[plate$Row == "B" & plate$Type == "Gene", c("1", "2", "3")] <- "GAPDH"
  plate[plate$Row == "C" & plate$Type == "Sample", c("1", "2", "3")] <- "Cnt"
  plate[plate$Row == "C" & plate$Type == "Gene", c("1", "2", "3")] <- "IL6"
  plate[plate$Row == "D" & plate$Type == "Sample", c("1", "2", "3")] <- "LPS"
  plate[plate$Row == "D" & plate$Type == "Gene", c("1", "2", "3")] <- "IL6"
  plate
}

# ---- Helpers ----
format_well_id <- function(x) {
  cleaned <- toupper(trimws(as.character(x)))
  matches <- regexec("^([A-H])0*([1-9]|1[0-2])$", cleaned)
  parts <- regmatches(cleaned, matches)
  out <- rep(NA_character_, length(cleaned))
  for (i in seq_along(parts)) {
    if (length(parts[[i]]) >= 3) {
      out[i] <- paste0(parts[[i]][2], sprintf("%02d", as.integer(parts[[i]][3])))
    }
  }
  out
}

clean_plate_value <- function(x) {
  cleaned <- trimws(as.character(x))
  cleaned[cleaned == ""] <- NA_character_
  cleaned
}

create_cq_template <- function() {
  grid <- expand.grid(
    Row = LETTERS[1:8],
    Column = sprintf("%02d", 1:12),
    stringsAsFactors = FALSE
  )
  data.frame(
    Well = paste0(grid$Row, grid$Column),
    Cq = NA_real_,
    stringsAsFactors = FALSE
  )
}

fill_cq_template <- function(values_df) {
  template <- create_cq_template()
  if (!is.null(values_df) && nrow(values_df)) {
    matches <- match(values_df$Well, template$Well)
    valid <- !is.na(matches)
    template$Cq[matches[valid]] <- values_df$Cq[valid]
  }
  template
}

create_plate_template <- function() {
  rows <- LETTERS[1:8]
  template <- data.frame(
    Row = rep(rows, each = 2),
    Type = rep(c("Sample", "Gene"), times = length(rows)),
    stringsAsFactors = FALSE
  )
  for (idx in 1:12) {
    template[[as.character(idx)]] <- NA_character_
  }
  template
}

empty_processed_table <- function() {
  data.frame(
    Well = character(0),
    Sample = character(0),
    Target = character(0),
    Cq = numeric(0),
    stringsAsFactors = FALSE
  )
}

normalize_plate_for_manual <- function(plate_df) {
  if (!"Plate1" %in% names(plate_df)) {
    return(plate_df)
  }
  cleaned <- as.data.frame(plate_df, stringsAsFactors = FALSE)
  plate_val <- trimws(as.character(cleaned$Plate1))
  cleaned$Row <- substr(toupper(plate_val), 1, 1)
  cleaned$Type <- ifelse(grepl("x$", tolower(plate_val)), "Sample", "Gene")
  extra_cols <- setdiff(names(cleaned), c("Plate1", "Row", "Type"))
  cleaned[, c("Row", "Type", extra_cols), drop = FALSE]
}

hot_input_to_tbl <- function(hot_value, fallback) {
  if (is.null(hot_value)) {
    return(fallback)
  }
  as.data.frame(rhandsontable::hot_to_r(hot_value), stringsAsFactors = FALSE)
}

parse_manual_colors <- function(value) {
  if (is.null(value) || !nzchar(value)) {
    return(character(0))
  }
  pieces <- unlist(strsplit(value, ","), use.names = FALSE)
  pieces <- trimws(pieces)
  pieces[pieces != ""]
}

get_palette <- function(choice, n, manual_colors) {
  if (n <= 0) {
    return(character(0))
  }
  normalize_palette <- function(palette, min_n = 20) {
    palette <- as.character(palette)
    palette <- palette[!is.na(palette) & nzchar(palette)]
    if (!length(palette)) {
      return(c("#ffffff", rep("#d0d0d0", max(min_n - 2, 0)), "#b3b3b3"))
    }
    if (length(palette) < min_n) {
      base_palette <- palette
      base_palette[1] <- "#ffffff"
      ramp_n <- min_n - 2
      ramp <- grDevices::colorRampPalette(base_palette)(ramp_n + 1)[-1]
      palette <- c("#ffffff", ramp, "#b3b3b3")
    } else {
      palette[1] <- "#ffffff"
      palette[length(palette)] <- "#b3b3b3"
    }
    palette
  }
  palettes <- list(
    grayscale = c(
      "#0d0d0d", "#262626", "#595959", "#7f7f7f",
      "#a1a1a1", "#bababa", "#d4d4d4", "#ededed"
    ),
    bright = c(
      "#003a7d", "#008dff", "#ff73b6", "#c701ff",
      "#4ecb8d", "#ff9d3a", "#f9e858", "#d83034"
    ),
    muted = c(
      "#c8c8c8", "#f0c571", "#59a89c", "#0b81a2",
      "#e25759", "#9d2c00", "#7E4794", "#36b700"
    ),
    pairs = c(
      "#8fd7d7", "#00b0be", "#ff8ca1", "#f45f74",
      "#bdd373", "#98c127", "#ffcd8e", "#ffb255"
    ),
    aesthetic1 = c(
      "#ffffff", "#a8c3e6", "#F08784", "#A6CDA7", "#D7CE99", "#B1D5CB"
    ),
    aesthetic2 = c(
      "#ffffff", "#F68CB4", "#85BFC0", "#A094BC", "#C4B5D6", "#B6E3F8"
    )
  )

  choice <- if (is.null(choice) || !nzchar(choice)) "bright" else tolower(choice)
  if (choice == "manual") {
    palette <- manual_colors
    if (!length(palette)) {
      palette <- palettes$bright
    }
  } else if (choice %in% names(palettes)) {
    palette <- palettes[[choice]]
  } else {
    palette <- palettes$bright
  }

  palette <- normalize_palette(palette, min_n = 20)
  rep(palette, length.out = n)
}

safe_sd <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) {
    return(NA_real_)
  }
  stats::sd(x)
}

darken_hex <- function(col, factor = 0.7) {
  if (is.null(col) || !length(col)) {
    return(character(0))
  }
  rgb <- grDevices::col2rgb(col)
  rgb <- rgb * factor
  rgb[rgb < 0] <- 0
  rgb[rgb > 255] <- 255
  grDevices::rgb(rgb[1, ], rgb[2, ], rgb[3, ], maxColorValue = 255)
}

read_cq_data <- function(source) {
  if (is.data.frame(source)) {
    raw <- as.data.frame(source, stringsAsFactors = FALSE)
  } else {
    if (!file.exists(source)) {
      stop("Cq file could not be found.")
    }
    extension <- tolower(tools::file_ext(source))
    if (extension %in% c("csv", "txt")) {
      raw <- read.csv(source, stringsAsFactors = FALSE, check.names = FALSE)
    } else {
      stop("Only CSV/TXT Cq files are supported in the Shinylive build.")
    }
  }
  if (!nrow(raw)) {
    stop("The Cq results file is empty.")
  }

  nm <- names(raw)
  well_idx <- which(grepl("well", nm, ignore.case = TRUE))[1]
  if (is.na(well_idx) || !length(well_idx)) {
    well_idx <- 1
  }
  cq_idx <- which(grepl("^(cq|ct)", nm, ignore.case = TRUE))[1]
  if (is.na(cq_idx) || !length(cq_idx)) {
    stop("Could not find a column that looks like Cq/Ct values.")
  }

  well_vals <- raw[[well_idx]]
  cq_vals <- raw[[cq_idx]]
  out <- data.frame(
    Well = format_well_id(well_vals),
    Cq = suppressWarnings(as.numeric(cq_vals)),
    stringsAsFactors = FALSE
  )
  out <- out[!is.na(out$Well) & !is.na(out$Cq), , drop = FALSE]
  out
}

read_plate_map <- function(source) {
  if (is.data.frame(source)) {
    plate <- as.data.frame(source, stringsAsFactors = FALSE)
  } else {
    if (!file.exists(source)) {
      stop("Plate file could not be found.")
    }
    extension <- tolower(tools::file_ext(source))
    if (!extension %in% c("csv", "txt")) {
      stop("Only CSV/TXT plate files are supported in the Shinylive build.")
    }
    plate <- read.csv(source, stringsAsFactors = FALSE, check.names = FALSE)
  }
  if (!ncol(plate)) {
    stop("Plate layout file has no columns.")
  }
  has_plate1 <- "Plate1" %in% names(plate)
  has_row_type <- all(c("Row", "Type") %in% names(plate))
  if (!has_plate1 && !has_row_type) {
    names(plate)[1] <- "Plate1"
  }
  plate
}

build_full_table <- function(cq_df, plate_df) {
  plate_df <- as.data.frame(plate_df, stringsAsFactors = FALSE)
  if ("Plate1" %in% names(plate_df)) {
    plate_val <- trimws(as.character(plate_df$Plate1))
    plate_df$RowLetter <- substr(plate_val, 1, 1)
    plate_df$Type <- ifelse(grepl("x$", tolower(plate_val)), "Sample", "Target")
    id_cols <- c("Plate1", "RowLetter", "Type")
  } else {
    row_vals <- trimws(as.character(plate_df$Row))
    plate_df$RowLetter <- substr(toupper(row_vals), 1, 1)
    type_raw <- tolower(trimws(as.character(plate_df$Type)))
    plate_df$Type <- ifelse(
      grepl("sample", type_raw),
      "Sample",
      ifelse(grepl("gene|target", type_raw), "Target", NA_character_)
    )
    id_cols <- c("Row", "Type", "RowLetter")
  }

  well_cols <- setdiff(names(plate_df), id_cols)
  if (!length(well_cols)) {
    stop("Plate layout is missing well columns.")
  }

  long_list <- lapply(seq_len(nrow(plate_df)), function(i) {
    values <- as.character(unlist(plate_df[i, well_cols], use.names = FALSE))
    data.frame(
      RowLetter = plate_df$RowLetter[i],
      Type = plate_df$Type[i],
      WellNumber = well_cols,
      Value = values,
      stringsAsFactors = FALSE
    )
  })
  tidy_plate <- do.call(rbind, long_list)
  tidy_plate <- tidy_plate[!is.na(tidy_plate$RowLetter) & !is.na(tidy_plate$Type), , drop = FALSE]

  sample_df <- tidy_plate[tidy_plate$Type == "Sample", c("RowLetter", "WellNumber", "Value"), drop = FALSE]
  target_df <- tidy_plate[tidy_plate$Type == "Target", c("RowLetter", "WellNumber", "Value"), drop = FALSE]
  names(sample_df)[names(sample_df) == "Value"] <- "Sample"
  names(target_df)[names(target_df) == "Value"] <- "Target"
  plate_wide <- merge(sample_df, target_df, by = c("RowLetter", "WellNumber"), all = TRUE)

  plate_wide$Sample <- clean_plate_value(plate_wide$Sample)
  plate_wide$Target <- clean_plate_value(plate_wide$Target)
  well_num <- suppressWarnings(as.integer(plate_wide$WellNumber))
  plate_wide$Well <- paste0(plate_wide$RowLetter, sprintf("%02d", well_num))
  final_df <- plate_wide[, c("Well", "Sample", "Target"), drop = FALSE]
  final_df <- final_df[order(final_df$Well), , drop = FALSE]

  cq_df <- as.data.frame(cq_df, stringsAsFactors = FALSE)
  final_df <- merge(final_df, cq_df, by = "Well", all.x = TRUE)
  final_df <- final_df[!is.na(final_df$Cq) & !is.na(final_df$Sample) & !is.na(final_df$Target), , drop = FALSE]
  final_df
}

calculate_fold_changes <- function(full_table, housekeeping_gene, control_sample) {
  full_table <- as.data.frame(full_table, stringsAsFactors = FALSE)
  housekeeping <- full_table[full_table$Target == housekeeping_gene, , drop = FALSE]

  if (!nrow(housekeeping)) {
    stop("No wells were annotated with the selected housekeeping gene.")
  }

  hk_mean <- aggregate(Cq ~ Sample, data = housekeeping, FUN = function(x) mean(x, na.rm = TRUE))
  names(hk_mean)[names(hk_mean) == "Cq"] <- "hk_mean"

  delta_df <- merge(full_table, hk_mean, by = "Sample", all.x = TRUE)
  delta_df$delta_Cq <- delta_df$Cq - delta_df$hk_mean
  delta_df <- delta_df[!is.na(delta_df$hk_mean) & !is.na(delta_df$delta_Cq), , drop = FALSE]

  control_ref <- delta_df[delta_df$Sample == control_sample, , drop = FALSE]
  control_delta <- aggregate(delta_Cq ~ Target, data = control_ref, FUN = function(x) mean(x, na.rm = TRUE))
  names(control_delta)[names(control_delta) == "delta_Cq"] <- "control_delta"
  control_delta <- control_delta[!is.na(control_delta$control_delta), , drop = FALSE]

  if (!nrow(control_delta)) {
    stop("The chosen control sample is missing housekeeping gene values.")
  }

  ordered_samples <- c("Cnt", "LPS", "AC", "AC+LPS")
  observed_samples <- unique(full_table$Sample)
  observed_samples <- observed_samples[!is.na(observed_samples)]
  sample_levels <- c(
    intersect(ordered_samples, observed_samples),
    setdiff(observed_samples, ordered_samples)
  )

  delta_delta_df <- merge(delta_df, control_delta, by = "Target", all.x = TRUE)
  delta_delta_df$delta_delta_Cq <- delta_delta_df$delta_Cq - delta_delta_df$control_delta
  delta_delta_df$fold_change <- 2^(-delta_delta_df$delta_delta_Cq)
  delta_delta_df$Sample <- factor(delta_delta_df$Sample, levels = sample_levels)
  delta_delta_df <- delta_delta_df[
    !is.na(delta_delta_df$control_delta) &
      !is.na(delta_delta_df$fold_change) &
      is.finite(delta_delta_df$fold_change),
    ,
    drop = FALSE
  ]

  if (!nrow(delta_delta_df)) {
    stop("No fold changes could be computed. Check housekeeping gene and control sample selections.")
  }

  summary_mean <- aggregate(fold_change ~ Sample + Target, data = delta_delta_df, FUN = function(x) mean(x, na.rm = TRUE))
  summary_sd <- aggregate(fold_change ~ Sample + Target, data = delta_delta_df, FUN = safe_sd)
  summary_n <- aggregate(fold_change ~ Sample + Target, data = delta_delta_df, FUN = length)
  names(summary_mean)[names(summary_mean) == "fold_change"] <- "mean_fold"
  names(summary_sd)[names(summary_sd) == "fold_change"] <- "sd_fold"
  names(summary_n)[names(summary_n) == "fold_change"] <- "replicates"

  summary_data <- merge(summary_mean, summary_sd, by = c("Sample", "Target"), all = TRUE)
  summary_data <- merge(summary_data, summary_n, by = c("Sample", "Target"), all = TRUE)
  summary_data$se_fold <- summary_data$sd_fold / sqrt(summary_data$replicates)
  summary_data <- summary_data[, c("Sample", "Target", "mean_fold", "se_fold", "replicates"), drop = FALSE]

  list(
    delta_delta = delta_delta_df,
    summary = summary_data
  )
}

# ---- UI ----
ui <- navbarPage(
  title = tags$div(
    class = "navbar-title",
    tags$span("CqPlotR", class = "app-title"),
    tags$span("qPCR analysis and visualization in one Shiny app", class = "app-tagline")
  ),
  header = tags$head(
    tags$link(
      rel = "stylesheet",
      href = "https://fonts.googleapis.com/css2?family=DM+Serif+Display&family=Source+Sans+3:wght@400;600;700&display=swap"
    ),
    tags$style(
      HTML(
        paste(
          ":root {",
          "  --nav-bg: #2D3142;",
          "  --nav-text: #ffffff;",
          "  --sidebar-bg: #ADACB5;",
          "  --table-shade: #ADACB5;",
          "  --page-ink: #2D3142;",
          "}",
          "body {",
          "  background: radial-gradient(1200px 600px at 12% 0%, #fbf8f3 0%, #f1ece3 40%, #e6e0d6 100%);",
          "  color: var(--page-ink);",
          "  font-family: \"Source Sans 3\", \"Segoe UI\", sans-serif;",
          "}",
          ".navbar-default {",
          "  background-color: var(--nav-bg);",
          "  border: none;",
          "}",
          ".navbar-default .navbar-brand,",
          ".navbar-default .navbar-brand:hover,",
          ".navbar-default .navbar-nav > li > a {",
          "  color: var(--nav-text);",
          "}",
          ".navbar-brand {",
          "  padding-top: 12px;",
          "  padding-bottom: 12px;",
          "  height: auto;",
          "  display: flex;",
          "  align-items: center;",
          "}",
          ".navbar-header {",
          "  display: flex;",
          "  align-items: center;",
          "}",
          ".navbar-collapse {",
          "  display: flex;",
          "  align-items: center;",
          "}",
          ".navbar-nav {",
          "  margin-left: auto;",
          "}",
          ".navbar-title {",
          "  display: flex;",
          "  flex-direction: column;",
          "  justify-content: center;",
          "  line-height: 1.05;",
          "}",
          ".app-title {",
          "  font-family: \"DM Serif Display\", serif;",
          "  font-size: 28px;",
          "  letter-spacing: 0.04em;",
          "}",
          ".app-tagline {",
          "  font-size: 13px;",
          "  opacity: 0.85;",
          "  margin-top: 2px;",
          "}",
          ".navbar-default .navbar-nav > li > a {",
          "  font-size: 20px;",
          "  font-weight: 700;",
          "  letter-spacing: 0.08em;",
          "  text-transform: uppercase;",
          "  padding: 20px 18px;",
          "}",
          ".navbar-default .navbar-nav > li > a:hover {",
          "  color: #f9e858;",
          "  background-color: transparent;",
          "}",
          ".navbar-default .navbar-nav > .active > a,",
          ".navbar-default .navbar-nav > .active > a:hover {",
          "  background-color: rgba(255, 255, 255, 0.12);",
          "  border-radius: 10px;",
          "  color: var(--nav-text);",
          "}",
          ".btn,",
          ".btn-default,",
          ".btn-secondary,",
          ".btn-primary {",
          "  background-color: var(--nav-bg);",
          "  border-color: var(--nav-bg);",
          "  color: #ffffff;",
          "  font-weight: 700;",
          "  letter-spacing: 0.02em;",
          "}",
          ".btn:hover,",
          ".btn:focus,",
          ".btn-default:hover,",
          ".btn-default:focus,",
          ".btn-secondary:hover,",
          ".btn-secondary:focus,",
          ".btn-primary:hover,",
          ".btn-primary:focus {",
          "  background-color: #1f2230;",
          "  border-color: #1f2230;",
          "  color: #ffffff;",
          "}",
          ".well {",
          "  background-color: var(--sidebar-bg);",
          "  border: none;",
          "  border-radius: 16px;",
          "  box-shadow: 0 12px 24px rgba(45, 49, 66, 0.15);",
          "  color: var(--page-ink);",
          "}",
          ".well .form-control,",
          ".well .selectize-input {",
          "  background-color: rgba(255, 255, 255, 0.85);",
          "  border: 1px solid rgba(45, 49, 66, 0.2);",
          "}",
          "details {",
          "  background-color: rgba(255, 255, 255, 0.75);",
          "  border-radius: 14px;",
          "  padding: 14px 18px;",
          "  margin-bottom: 14px;",
          "  box-shadow: 0 10px 24px rgba(45, 49, 66, 0.12);",
          "}",
          "details > summary {",
          "  font-weight: 700;",
          "  letter-spacing: 0.03em;",
          "}",
          "table.dataTable tbody tr {",
          "  background-color: var(--table-shade);",
          "}",
          "table.dataTable thead th {",
          "  background-color: #2D3142;",
          "  color: #ffffff;",
          "}",
          ".table.dataTable {",
          "  border: 1px solid rgba(45, 49, 66, 0.35);",
          "}",
          ".handsontable .htCore td,",
          ".handsontable .htCore th {",
          "  background-color: #ededed;",
          "  border: 1px solid rgba(45, 49, 66, 0.35);",
          "}",
          ".handsontable thead th,",
          ".handsontable .ht_clone_top th,",
          ".handsontable .ht_clone_left th {",
          "  background-color: #2D3142;",
          "  color: #ffffff;",
          "  border: 1px solid rgba(45, 49, 66, 0.35);",
          "}",
          ".handsontable .htCore td {",
          "  color: var(--page-ink);",
          "}",
          ".nav-tabs > li > a,",
          ".nav-tabs .nav-link {",
          "  color: #1f2230;",
          "  border: 1px solid #2D3142;",
          "  font-weight: 700;",
          "  letter-spacing: 0.02em;",
          "}",
          ".nav-tabs > li.active > a,",
          ".nav-tabs > li.active > a:focus,",
          ".nav-tabs > li.active > a:hover,",
          ".nav-tabs .nav-link.active {",
          "  background-color: #2D3142;",
          "  color: #ffffff;",
          "  border-color: #2D3142;",
          "}",
          ".nav-tabs > li > a:hover,",
          ".nav-tabs .nav-link:hover {",
          "  background-color: rgba(45, 49, 66, 0.08);",
          "  border-color: #2D3142;",
          "  color: #1f2230;",
          "}",
          ".preset-toggle .radio-inline {",
          "  background-color: #b3b3b3;",
          "  border: 1px solid #8c8c8c;",
          "  border-radius: 999px;",
          "  color: #1f2230;",
          "  font-weight: 700;",
          "  letter-spacing: 0.02em;",
          "  margin-right: 6px;",
          "  padding: 6px 12px;",
          "  cursor: pointer;",
          "  user-select: none;",
          "}",
          ".preset-toggle .radio-inline input {",
          "  margin-right: 6px;",
          "}",
          ".preset-toggle .radio-inline.is-selected {",
          "  background-color: #2D3142;",
          "  border-color: #2D3142;",
          "  color: #ffffff;",
          "}",
          "#fold_plot img,",
          "#fold_plot canvas {",
          "  width: 100% !important;",
          "  height: 100% !important;",
          "}",
          ".code-box {",
          "  margin-top: 14px;",
          "}",
          ".code-box pre {",
          "  max-height: 280px;",
          "  overflow: auto;",
          "  white-space: pre;",
          "  background-color: #f6f6f6;",
          "  border-radius: 8px;",
          "  padding: 10px 12px;",
          "}",
          "hr {",
          "  border-top: 1px solid rgba(45, 49, 66, 0.2);",
          "}",
          sep = "\n"
        )
      )
    ),
    tags$script(
      HTML(
        "$(document).on('change', 'input[name=paper_preset]', function() {",
        "  var $labels = $('input[name=paper_preset]').closest('label');",
        "  $labels.removeClass('is-selected');",
        "  $(this).closest('label').addClass('is-selected');",
        "});",
        "$(function(){",
        "  var $checked = $('input[name=paper_preset]:checked');",
        "  if ($checked.length) {",
        "    $checked.closest('label').addClass('is-selected');",
        "  }",
        "});"
      )
    ),
    tags$script(
      HTML(
        "(function(){",
        "  function downloadBlob(blob, filename){",
        "    var link = document.createElement('a');",
        "    link.href = URL.createObjectURL(blob);",
        "    link.download = filename;",
        "    document.body.appendChild(link);",
        "    link.click();",
        "    document.body.removeChild(link);",
        "    setTimeout(function(){ URL.revokeObjectURL(link.href); }, 1000);",
        "  }",
        "  function serializeSvg(svg){",
        "    var serializer = new XMLSerializer();",
        "    var svgText = serializer.serializeToString(svg);",
        "    if (!svgText.match(/^<svg[^>]+xmlns=/)) {",
        "      svgText = svgText.replace('<svg', '<svg xmlns=\"http://www.w3.org/2000/svg\"');",
        "    }",
        "    if (!svgText.match(/^<svg[^>]+xmlns:xlink=/)) {",
        "      svgText = svgText.replace('<svg', '<svg xmlns:xlink=\"http://www.w3.org/1999/xlink\"');",
        "    }",
        "    return svgText;",
        "  }",
        "  function getPlotNode(){",
        "    return document.getElementById('fold_plot');",
        "  }",
        "  function getSvg(node){",
        "    return node ? node.querySelector('svg') : null;",
        "  }",
        "  function getImg(node){",
        "    return node ? node.querySelector('img') : null;",
        "  }",
        "  function elementSize(el){",
        "    var rect = el.getBoundingClientRect();",
        "    var w = el.naturalWidth || el.width || 0;",
        "    var h = el.naturalHeight || el.height || 0;",
        "    if ((!w || !h) && el.viewBox && el.viewBox.baseVal) {",
        "      w = w || el.viewBox.baseVal.width;",
        "      h = h || el.viewBox.baseVal.height;",
        "    }",
        "    if (!w || !h) {",
        "      w = rect.width;",
        "      h = rect.height;",
        "    }",
        "    return {width: w || 800, height: h || 600};",
        "  }",
        "  function svgToImage(svgText, cb){",
        "    var img = new Image();",
        "    img.onload = function(){ cb(null, img); };",
        "    img.onerror = function(){ cb(new Error('Could not load SVG.')); };",
        "    img.src = 'data:image/svg+xml;charset=utf-8,' + encodeURIComponent(svgText);",
        "  }",
        "  function imageToCanvas(img, width, height){",
        "    var canvas = document.createElement('canvas');",
        "    canvas.width = Math.round(width);",
        "    canvas.height = Math.round(height);",
        "    var ctx = canvas.getContext('2d');",
        "    ctx.drawImage(img, 0, 0, width, height);",
        "    return canvas;",
        "  }",
        "  function downloadPng(){",
        "    var node = getPlotNode();",
        "    var svg = getSvg(node);",
        "    var img = getImg(node);",
        "    var source = svg || img;",
        "    if (!source) {",
        "      alert('Plot not found for export.');",
        "      return;",
        "    }",
        "    var size = elementSize(source);",
        "    var useImage = function(image){",
        "      var canvas = imageToCanvas(image, size.width, size.height);",
        "      canvas.toBlob(function(blob){",
        "        if (!blob) {",
        "          alert('Export failed.');",
        "          return;",
        "        }",
        "        downloadBlob(blob, 'qpcr_plot.png');",
        "      }, 'image/png');",
        "    };",
        "    if (svg) {",
        "      var svgText = serializeSvg(svg);",
        "      svgToImage(svgText, function(err, image){",
        "        if (err) {",
        "          alert('SVG export failed.');",
        "          return;",
        "        }",
        "        useImage(image);",
        "      });",
        "    } else if (img) {",
        "      if (img.complete) {",
        "        useImage(img);",
        "      } else {",
        "        img.onload = function(){ useImage(img); };",
        "        img.onerror = function(){ alert('Image export failed.'); };",
        "      }",
        "    }",
        "  }",
        "  document.addEventListener('click', function(evt){",
        "    var btn = evt.target && evt.target.closest ? evt.target.closest('#download_plot_client') : null;",
        "    if (btn) {",
        "      evt.preventDefault();",
        "      downloadPng();",
        "    }",
        "  });",
        "})();"
      )
    ),
    
  ),
  theme = bslib::bs_theme(bootswatch = "sandstone"),
  tabPanel(
    "Data",
    sidebarLayout(
      sidebarPanel(
        p(
          "Pick a data source below, then upload files or work directly inside the manual entry tables. ",
          "This app works with 96 well format BIORAD Real Time System Quantification Cq Results."
        ),
        uiOutput("data_error"),
        radioButtons(
          "data_mode",
          label = "Data source",
          choices = c(
            "Use bundled example" = "example",
            "Upload or enter data" = "upload_enter"
          ),
          selected = "example"
        ),
        conditionalPanel(
          condition = "input.data_mode == 'upload_enter'",
          helpText("Choose how you want to provide Cq values and the plate layout."),
          radioButtons(
            "manual_mode",
            label = "Upload options",
            choices = c(
              "Upload both files" = "upload_both",
              "Upload Cq file, enter plate manually" = "upload_cq",
              "Upload plate file, enter Cq manually" = "upload_plate",
              "Enter both manually" = "enter_both"
            ),
            selected = "upload_both"
          ),
          conditionalPanel(
            condition = "input.manual_mode == 'upload_both' || input.manual_mode == 'upload_cq'",
            fileInput(
              "cq_file",
              label = "Cq results (.csv or .txt)",
              accept = c(".csv", ".txt")
            )
          ),
          conditionalPanel(
            condition = "input.manual_mode == 'upload_both' || input.manual_mode == 'upload_plate'",
            fileInput(
              "plate_file",
              label = "Plate layout (.csv or .txt)",
              accept = c(".csv", ".txt")
            )
          )
        )
      ),
      mainPanel(
        tags$details(
          open = TRUE,
          tags$summary("Manual entry"),
          p("Edit wells directly or paste from a spreadsheet. Use the buttons to reset or load the bundled example."),
          conditionalPanel(
            condition = "input.data_mode == 'upload_enter' && (input.manual_mode == 'upload_plate' || input.manual_mode == 'enter_both')",
            h4("Cq value table"),
            fluidRow(
              column(
                width = 6,
                actionButton("reset_manual_cq", "Reset blank Cq template", class = "btn-secondary btn-sm")
              ),
              column(
                width = 6,
                actionButton("load_example_manual_cq", "Load example Cq data", class = "btn-secondary btn-sm")
              )
            ),
            tags$br(),
            rHandsontableOutput("manual_cq_table"),
            tags$hr()
          ),
          conditionalPanel(
            condition = "input.data_mode == 'upload_enter' && (input.manual_mode == 'upload_cq' || input.manual_mode == 'enter_both')",
            h4("Plate layout table"),
            p("Fill the Sample and Gene rows for each lettered plate row."),
            fluidRow(
              column(
                width = 6,
                actionButton("reset_manual_plate", "Reset plate template", class = "btn-secondary btn-sm")
              ),
              column(
                width = 6,
                actionButton("load_example_manual_plate", "Load example plate", class = "btn-secondary btn-sm")
              )
            ),
            tags$br(),
            rHandsontableOutput("manual_plate_table")
          )
        ),
        tags$hr(),
        tags$details(
          open = TRUE,
          tags$summary("Processed wells"),
          DTOutput("full_table")
        )
      )
    )
  ),
  tabPanel(
    "Summary",
    sidebarLayout(
      sidebarPanel(
        selectInput("control_sample", "Control sample", choices = NULL),
        selectInput("housekeeping_gene", "Housekeeping gene", choices = NULL)
      ),
      mainPanel(
        uiOutput("summary_message"),
        DTOutput("summary_table"),
        tags$br(),
        downloadButton("download_summary", "Download summary (CSV)")
      )
    )
  ),
  tabPanel(
    "Barplots",
    fluidRow(
      column(
        width = 3,
        div(
          class = "well",
          checkboxGroupInput("target_filter", "Targets to display", choices = NULL),
          selectizeInput(
            "sample_order",
            label = "Sample order (drag to reorder, delete to remove)",
            choices = character(0),
            selected = character(0),
            multiple = TRUE,
            options = list(
              plugins = list("drag_drop", "remove_button"),
              persist = TRUE,
              create = FALSE,
              sortField = I("null")
            )
          ),
          bslib::navset_card_tab(
            tabPanel(
              "Size",
              sliderInput("plot_height", "Plot height (px)", min = 400, max = 900, value = 500, step = 50),
              sliderInput("plot_aspect", "Aspect ratio (height/width)", min = 0.6, max = 1.4, value = 0.9, step = 0.05),
              sliderInput("facet_spacing", "Facet spacing (lines)", min = 0.2, max = 2, value = 1, step = 0.1),
              sliderInput("x_expand", "X padding", min = 0, max = 1, value = 0.2, step = 0.05)
            ),
            tabPanel(
              "Aesthetics",
              div(
                class = "preset-row",
                strong("Presets"),
                div(
                  class = "preset-toggle",
                  radioButtons(
                    "paper_preset",
                    label = NULL,
                    choices = c(
                      "Reviewer1" = "paper_ready_1",
                      "Reviewer2" = "paper_ready_2",
                      "Reviewer3" = "paper_ready_3"
                    ),
                    selected = "paper_ready_1",
                    inline = TRUE
                  )
                )
              ),
              tags$br(),
              checkboxInput("show_points", "Show points", value = TRUE),
              sliderInput("point_size", "Point size", min = 1, max = 5, value = 3, step = 0.5),
              sliderInput("jitter_width", "Jitter width", min = 0, max = 0.5, value = 0.25, step = 0.05),
              checkboxInput("show_sample_labels", "Show sample labels", value = TRUE),
              checkboxInput("legend_on", "Legend on", value = FALSE),
              selectInput(
                "palette_choice",
                "Colors",
                choices = c(
                  "Grayscale" = "grayscale",
                  "Bright" = "bright",
                  "Muted" = "muted",
                  "Alt Light/Dark Pairs" = "pairs",
                  "Publish" = "aesthetic1",
                  "Perish" = "aesthetic2",
                  "Custom manual" = "manual"
                ),
                selected = "aesthetic1"
              ),
              conditionalPanel(
                condition = "input.palette_choice == 'manual'",
                textInput("manual_palette", "Manual colors (comma-separated)", "#003a7d,#008dff,#ff73b6")
              ),
              selectInput(
                "theme_choice",
                "Theme",
                choices = c(
                  "Gray" = "gray",
                  "BW" = "bw",
                  "Linedraw" = "linedraw",
                  "Light" = "light",
                  "Dark" = "dark",
                  "Minimal" = "minimal",
                  "Classic" = "classic",
                  "Void" = "void"
                ),
                selected = "classic"
              ),
              sliderInput("base_font_size", "Base font size", min = 12, max = 60, value = 20, step = 2)
            ),
            tabPanel(
              "Axes",
              checkboxInput("rotate_x_labels", "Rotate x labels", value = FALSE),
              conditionalPanel(
                condition = "input.rotate_x_labels",
                sliderInput("x_label_angle", "X label angle", min = 0, max = 90, value = 45, step = 5)
              ),
              radioButtons(
                "y_limit_mode",
                "Y-axis limits",
                choices = c("Auto" = "auto", "Manual" = "manual"),
                inline = TRUE
              ),
              conditionalPanel(
                condition = "input.y_limit_mode == 'manual'",
                numericInput("y_min", "Y min", value = 0),
                numericInput("y_max", "Y max", value = 2)
              ),
              checkboxInput("log_scale", "Show fold change on log2 scale", value = FALSE)
            )
          ),
          actionButton("render_plot", "Refresh plot", class = "btn-primary")
        )
      ),
      column(
        width = 6,
        uiOutput("plot_message"),
        uiOutput("plot_ui")
      ),
      column(
        width = 3,
        div(
          class = "well",
          h4("Download plot"),
          radioButtons(
            "export_quality",
            "Render quality",
            choices = c("Low" = "1", "Medium" = "2", "High" = "3"),
            selected = "2",
            inline = TRUE
          ),
          actionButton("download_plot_client", "Download PNG", class = "btn-primary")
        ),
        div(
          class = "well code-box",
          h4("Reproducible ggplot2 code"),
          verbatimTextOutput("plot_code")
        )
      )
    )
  )
)

# ---- Server ----
server <- function(input, output, session) {
  # Manual entry state
  manual_cq <- reactiveVal(create_cq_template())
  manual_plate <- reactiveVal(create_plate_template())
  preset_style <- reactiveVal("paper_ready_1")
  sample_order_state <- reactiveVal(character(0))
  last_sample_set <- reactiveVal(character(0))
  data_error_msg <- reactiveVal(NULL)

  apply_paper_ready_preset <- function(choice) {
    preset_style(choice)
    updateCheckboxInput(session, "show_points", value = TRUE)
    updateSliderInput(session, "point_size", value = 4)
    updateSliderInput(session, "jitter_width", value = 0.12)
    updateCheckboxInput(session, "show_sample_labels", value = TRUE)
    updateCheckboxInput(session, "legend_on", value = FALSE)
    updateSelectInput(session, "theme_choice", selected = "classic")
    updateSliderInput(session, "base_font_size", value = 16)
    updateCheckboxInput(session, "rotate_x_labels", value = TRUE)
    updateSliderInput(session, "x_label_angle", value = 35)
    updateSliderInput(session, "x_expand", value = 0.35)
    updateSliderInput(session, "plot_aspect", value = 0.85)
    updateSliderInput(session, "facet_spacing", value = 0.8)
  }

  observeEvent(input$paper_preset, {
    choice <- input$paper_preset
    if (is.null(choice) || !nzchar(choice)) {
      return()
    }
    apply_paper_ready_preset(choice)
  }, ignoreInit = TRUE)

  session$onFlushed(function() {
    apply_paper_ready_preset("paper_ready_1")
  }, once = TRUE)

  observeEvent(input$reset_manual_cq, {
    manual_cq(create_cq_template())
  })

  observeEvent(input$load_example_manual_cq, {
    example_vals <- read_cq_data(example_cq_data())
    manual_cq(fill_cq_template(example_vals))
  })

  observeEvent(input$reset_manual_plate, {
    manual_plate(create_plate_template())
  })

  observeEvent(input$load_example_manual_plate, {
    example_plate <- read_plate_map(example_plate_data())
    manual_plate(normalize_plate_for_manual(example_plate))
  })

  observeEvent(input$manual_cq_table, {
    manual_cq(hot_input_to_tbl(input$manual_cq_table, manual_cq()))
  })

  observeEvent(input$manual_plate_table, {
    manual_plate(hot_input_to_tbl(input$manual_plate_table, manual_plate()))
  })

  # Data source selection
  cq_values <- reactive({
    mode <- req(input$data_mode)
    if (mode == "example") {
      read_cq_data(example_cq_data())
    } else {
      manual_mode <- req(input$manual_mode)
      if (manual_mode %in% c("upload_both", "upload_cq")) {
        req(input$cq_file)
        read_cq_data(input$cq_file$datapath)
      } else {
        read_cq_data(manual_cq())
      }
    }
  })

  plate_layout <- reactive({
    mode <- req(input$data_mode)
    if (mode == "example") {
      read_plate_map(example_plate_data())
    } else {
      manual_mode <- req(input$manual_mode)
      if (manual_mode %in% c("upload_both", "upload_plate")) {
        req(input$plate_file)
        read_plate_map(input$plate_file$datapath)
      } else {
        read_plate_map(manual_plate())
      }
    }
  })

  # Derived tables
  processed_wells <- reactive({
    data_error_msg(NULL)
    errors <- character(0)

    cq_df <- tryCatch(
      cq_values(),
      error = function(e) {
        errors <<- c(errors, conditionMessage(e))
        NULL
      }
    )

    plate_df <- tryCatch(
      plate_layout(),
      error = function(e) {
        errors <<- c(errors, conditionMessage(e))
        NULL
      }
    )

    if (length(errors)) {
      data_error_msg(paste(unique(errors), collapse = " "))
      return(empty_processed_table())
    }

    tryCatch(
      build_full_table(cq_df, plate_df),
      error = function(e) {
        data_error_msg(conditionMessage(e))
        empty_processed_table()
      }
    )
  })

  # Input updates based on processed wells
  observeEvent(processed_wells(), {
    data <- processed_wells()
    samples <- sort(unique(as.character(data$Sample)))
    samples <- samples[!is.na(samples) & nzchar(samples)]

    if (length(samples)) {
      selected_sample <- input$control_sample
      if (is.null(selected_sample) || !selected_sample %in% samples) {
        selected_sample <- samples[1]
      }
      updateSelectInput(session, "control_sample", choices = samples, selected = selected_sample)
    } else {
      updateSelectInput(session, "control_sample", choices = character(0), selected = NULL)
    }

    targets <- sort(unique(as.character(data$Target)))
    targets <- targets[!is.na(targets) & nzchar(targets)]

    if (length(targets)) {
      preferred_hk <- if ("GAPDH" %in% targets) "GAPDH" else targets[1]
      selected_hk <- isolate(input$housekeeping_gene)
      if (is.null(selected_hk) || !selected_hk %in% targets) {
        selected_hk <- preferred_hk
      }
      updateSelectInput(session, "housekeeping_gene", choices = targets, selected = selected_hk)
      existing_target_selection <- isolate(input$target_filter)
      new_selection <- intersect(existing_target_selection, targets)
      if (!length(new_selection)) {
        new_selection <- targets
      }
      updateCheckboxGroupInput(session, "target_filter", choices = targets, selected = new_selection)
    } else {
      updateSelectInput(session, "housekeeping_gene", choices = character(0), selected = NULL)
      updateCheckboxGroupInput(session, "target_filter", choices = character(0), selected = NULL)
    }

    last_samples <- last_sample_set()
    if (!identical(samples, last_samples)) {
      sample_order_state(samples)
      updateSelectizeInput(
        session,
        "sample_order",
        choices = samples,
        selected = samples,
        server = FALSE
      )
      last_sample_set(samples)
    }
  }, ignoreNULL = FALSE)

  observeEvent(input$sample_order, {
    if (is.null(input$sample_order)) {
      sample_order_state(character(0))
      return()
    }
    sample_order_state(as.character(input$sample_order))
  }, ignoreInit = TRUE)

  fold_results <- reactive({
    data <- processed_wells()
    if (!nrow(data)) {
      return(list(
        ok = FALSE,
        message = "No wells with Sample/Target labels and Cq values are available. Upload data or enter values inside the manual entry tables.",
        data = data
      ))
    }

    housekeeping_gene <- input$housekeeping_gene
    control_sample <- input$control_sample

    if (is.null(housekeeping_gene) || !nzchar(housekeeping_gene)) {
      return(list(ok = FALSE, message = "Select a housekeeping gene.", data = data))
    }
    if (is.null(control_sample) || !nzchar(control_sample)) {
      return(list(ok = FALSE, message = "Select a control sample.", data = data))
    }
    if (!housekeeping_gene %in% data$Target) {
      return(list(ok = FALSE, message = "Housekeeping gene is not found in the processed wells.", data = data))
    }
    if (!control_sample %in% data$Sample) {
      return(list(ok = FALSE, message = "Control sample is not found in the processed wells.", data = data))
    }

    results <- tryCatch(
      calculate_fold_changes(data, housekeeping_gene, control_sample),
      error = function(e) e
    )

    if (inherits(results, "error")) {
      return(list(ok = FALSE, message = conditionMessage(results), data = data))
    }

    list(ok = TRUE, value = results, data = data)
  })

  summary_all <- reactive({
    res <- fold_results()
    if (!res$ok) {
      return(list(ok = FALSE, message = res$message, data = res$data))
    }
    list(ok = TRUE, data = res$value$summary)
  })

  filtered_summary <- reactive({
    summary_state <- summary_all()
    if (!summary_state$ok) {
      return(list(ok = FALSE, message = summary_state$message, data = summary_state$data))
    }
    summary <- summary_state$data
    targets <- input$target_filter
    if (!is.null(targets) && length(targets)) {
      summary <- summary[summary$Target %in% targets, , drop = FALSE]
    }
    list(ok = TRUE, data = summary)
  })

  filtered_points <- reactive({
    res <- fold_results()
    if (!res$ok) {
      return(list(ok = FALSE, message = res$message, data = res$data))
    }
    points <- res$value$delta_delta
    targets <- input$target_filter
    if (!is.null(targets) && length(targets)) {
      points <- points[points$Target %in% targets, , drop = FALSE]
    }
    list(ok = TRUE, data = points)
  })

  plot_object <- reactive({
    input$render_plot
    summary_state <- filtered_summary()
    validate(need(summary_state$ok, summary_state$message))
    points_state <- filtered_points()
    validate(need(points_state$ok, points_state$message))
    summary <- summary_state$data
    points <- points_state$data
    validate(need(nrow(summary) > 0, "Select at least one target to display."))

    preset_choice <- preset_style()
    paper_ready <- !is.null(preset_choice) && preset_choice %in% c("paper_ready_1", "paper_ready_2", "paper_ready_3")
    paper_ready_2 <- identical(preset_choice, "paper_ready_2")
    paper_ready_3 <- identical(preset_choice, "paper_ready_3")
    size_scale <- suppressWarnings(as.numeric(input$export_quality))
    if (!is.finite(size_scale) || size_scale <= 0) {
      size_scale <- 1
    }
    scaled_base_size <- input$base_font_size * size_scale
    scaled_point_size <- input$point_size * size_scale
    scaled_stroke <- 0.5 * size_scale
    available_samples <- as.character(sort(unique(summary$Sample)))
    sample_levels <- sample_order_state()
    if (!length(sample_levels)) {
      sample_levels <- available_samples
    }
    if (length(sample_levels)) {
      summary <- summary[as.character(summary$Sample) %in% sample_levels, , drop = FALSE]
      points <- points[as.character(points$Sample) %in% sample_levels, , drop = FALSE]
    }
    summary$Sample <- factor(as.character(summary$Sample), levels = sample_levels)
    points$Sample <- factor(as.character(points$Sample), levels = sample_levels)
    manual_colors <- parse_manual_colors(input$manual_palette)
    palette_values <- get_palette(input$palette_choice, length(sample_levels), manual_colors)
    color_values <- setNames(palette_values, sample_levels)
    fill_values <- color_values
    point_fill_values <- NULL
    if (paper_ready) {
      fill_values <- setNames(as.character(color_values), sample_levels)
      if (paper_ready_3) {
        fill_values[] <- "#ffffff"
      }
      point_base_colors <- if (paper_ready_3) color_values else fill_values
      point_fill_values <- setNames(darken_hex(point_base_colors, factor = 0.7), sample_levels)
      if (paper_ready_2) {
        point_fill_values[] <- "#ffffff"
      }
      summary$bar_fill <- fill_values[as.character(summary$Sample)]
      points$point_fill <- point_fill_values[as.character(points$Sample)]
    }

    y_limits <- NULL
    if (identical(input$y_limit_mode, "manual")) {
      y_min <- input$y_min
      y_max <- input$y_max
      validate(need(is.finite(y_min) && is.finite(y_max), "Provide numeric Y-axis limits."))
      validate(need(y_min < y_max, "Y min must be smaller than Y max."))
      if (isTRUE(input$log_scale)) {
        validate(need(y_min > 0, "Log scale requires Y min > 0."))
      }
      y_limits <- c(y_min, y_max)
    }

    if (isTRUE(input$log_scale)) {
      log_values <- c(points$fold_change, summary$mean_fold - summary$se_fold, summary$mean_fold + summary$se_fold)
      log_values <- log_values[is.finite(log_values)]
      validate(need(length(log_values), "Log scale requires finite fold change values."))
      validate(need(all(log_values > 0), "Log scale requires all fold-change values to be > 0."))
    }

    axis_color <- if (paper_ready) "#000000" else "#222222"
    axis_line_width <- (if (paper_ready) 0.7 else 0.4) * size_scale
    angle <- if (isTRUE(input$rotate_x_labels)) input$x_label_angle else 0
    x_text <- if (isTRUE(input$show_sample_labels)) {
      element_text(color = axis_color, angle = angle, hjust = if (angle > 0) 1 else 0.5)
    } else {
      element_blank()
    }
    x_ticks <- if (isTRUE(input$show_sample_labels)) {
      element_line(color = axis_color, linewidth = axis_line_width)
    } else {
      element_blank()
    }

    base_family <- "sans"
    theme_base <- switch(
      input$theme_choice,
      gray = theme_gray(base_size = scaled_base_size, base_family = base_family),
      bw = theme_bw(base_size = scaled_base_size, base_family = base_family),
      linedraw = theme_linedraw(base_size = scaled_base_size, base_family = base_family),
      light = theme_light(base_size = scaled_base_size, base_family = base_family),
      dark = theme_dark(base_size = scaled_base_size, base_family = base_family),
      minimal = theme_minimal(base_size = scaled_base_size, base_family = base_family),
      classic = theme_classic(base_size = scaled_base_size, base_family = base_family),
      void = theme_void(base_size = scaled_base_size, base_family = base_family),
      theme_classic(base_size = scaled_base_size, base_family = base_family)
    )

    bar_width <- if (paper_ready) 0.45 else 0.55
    bar_line_width <- (if (paper_ready) 0.4 else 0.7) * size_scale
    errorbar_line_width <- (if (paper_ready) 0.5 else 0.6) * size_scale
    panel_fill <- "white"

    if (paper_ready) {
      p <- ggplot() +
        geom_col(
          data = summary,
          aes(x = Sample, y = mean_fold, fill = bar_fill),
          width = bar_width,
          color = "black",
          linewidth = bar_line_width
        ) +
        geom_errorbar(
          data = summary,
          aes(
            x = Sample,
            ymin = mean_fold - se_fold,
            ymax = mean_fold + se_fold
          ),
          width = 0.2,
          linewidth = errorbar_line_width,
          color = "black"
        )
    } else {
      p <- ggplot() +
        geom_col(
          data = summary,
          aes(x = Sample, y = mean_fold, color = Sample),
          width = bar_width,
          fill = "white",
          linewidth = bar_line_width
        ) +
        geom_errorbar(
          data = summary,
          aes(
            x = Sample,
            ymin = mean_fold - se_fold,
            ymax = mean_fold + se_fold,
            color = Sample
          ),
          width = 0.2,
          linewidth = errorbar_line_width
        )
    }

    if (isTRUE(input$show_points)) {
      jitter_width <- input$jitter_width
      if (paper_ready) {
        if (paper_ready_2) {
          p <- p +
            geom_point(
              data = points,
              aes(
                x = Sample,
                y = fold_change,
                shape = Sample
              ),
              position = position_jitterdodge(
                jitter.width = jitter_width,
                dodge.width = 0.55
              ),
              size = scaled_point_size,
              stroke = scaled_stroke,
              color = "black",
              fill = "white"
            )
        } else {
          p <- p +
            geom_point(
              data = points,
              aes(
                x = Sample,
                y = fold_change,
                fill = point_fill
              ),
              position = position_jitterdodge(
                jitter.width = jitter_width,
                dodge.width = 0.55
              ),
              shape = 21,
              size = scaled_point_size,
              stroke = scaled_stroke,
              color = "black"
            )
        }
      } else {
        p <- p +
          geom_jitter(
            data = points,
            aes(x = Sample, y = fold_change, color = Sample),
            width = jitter_width,
            height = 0,
            shape = 16,
            size = scaled_point_size
          )
      }
    }

    facet_formals <- names(formals(ggplot2::facet_wrap))
    facet_layer <- if ("axes" %in% facet_formals) {
      ggplot2::facet_wrap(~Target, scales = "free_y", axes = "all")
    } else if ("axis" %in% facet_formals) {
      ggplot2::facet_wrap(~Target, scales = "free_y", axis = "all")
    } else if (requireNamespace("ggh4x", quietly = TRUE)) {
      ggh4x::facet_wrap2(~Target, scales = "free_y", axes = "all", remove_labels = "none")
    } else {
      ggplot2::facet_wrap(~Target, scales = "free_y")
    }

    p <- p +
      facet_layer +
      scale_x_discrete(expand = expansion(mult = c(input$x_expand, input$x_expand))) +
      coord_cartesian(ylim = y_limits, clip = "off") +
      labs(
        x = "Sample",
        y = expression("Fold Change ("*2^{-Delta*Delta*Cq}*")")
      ) +
      theme_base +
      theme(
        legend.position = if (isTRUE(input$legend_on)) "right" else "none",
        panel.background = element_rect(fill = panel_fill, color = NA),
        plot.background = element_rect(fill = "white", color = NA),
        axis.line = element_line(color = axis_color, linewidth = axis_line_width),
        axis.ticks.x = x_ticks,
        axis.ticks.y = element_line(color = axis_color, linewidth = axis_line_width),
        axis.text.x = x_text,
        axis.text.y = element_text(color = axis_color),
        axis.title = element_text(color = axis_color),
        panel.spacing = grid::unit(input$facet_spacing, "lines"),
        aspect.ratio = input$plot_aspect,
        strip.background = element_blank(),
        strip.text = element_text(face = "bold")
      )

    if (paper_ready) {
      p <- p + scale_fill_identity(guide = "none")
    } else {
      p <- p + scale_color_manual(values = color_values, drop = FALSE)
    }

    if (paper_ready_2) {
      shape_values <- rep(c(21, 22, 23, 24, 25), length.out = length(sample_levels))
      p <- p + scale_shape_manual(values = setNames(shape_values, sample_levels))
    }

    if (paper_ready) {
      p <- p + theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
      )
    }

    if (paper_ready) {
      # no significance brackets or asterisks
    }

    if (isTRUE(input$log_scale)) {
      p <- p + scale_y_continuous(trans = "log2")
    } else if (is.null(y_limits)) {
      p <- p + scale_y_continuous(expand = c(0, 0), limits = c(0, NA))
    } else {
      p <- p + scale_y_continuous(expand = c(0, 0))
    }

    p
  })

  output$fold_plot <- renderPlot({
    plot_object()
  },
  width = function() {
    scale_val <- suppressWarnings(as.numeric(input$export_quality))
    if (!is.finite(scale_val) || scale_val <= 0) {
      scale_val <- 1
    }
    base_width <- session$clientData$output_fold_plot_width
    if (is.null(base_width) || !is.finite(base_width) || base_width <= 0) {
      base_width <- 800
    }
    base_width * scale_val
  },
  height = function() {
    scale_val <- suppressWarnings(as.numeric(input$export_quality))
    if (!is.finite(scale_val) || scale_val <= 0) {
      scale_val <- 1
    }
    base_height <- input$plot_height
    if (is.null(base_height) || !is.finite(base_height) || base_height <= 0) {
      base_height <- 500
    }
    base_height * scale_val
  })

  output$plot_ui <- renderUI({
    plotOutput("fold_plot", height = input$plot_height)
  })

  # Client-side downloads are handled in JavaScript for Shinylive compatibility.

  output$plot_code <- renderText({
    summary_state <- filtered_summary()
    points_state <- filtered_points()
    if (!summary_state$ok || !points_state$ok) {
      return("Plot code will appear once valid data is available.")
    }
    summary <- summary_state$data
    points <- points_state$data
    if (!nrow(summary)) {
      return("Plot code will appear once valid data is available.")
    }

    preset_choice <- preset_style()
    paper_ready <- !is.null(preset_choice) && preset_choice %in% c("paper_ready_1", "paper_ready_2", "paper_ready_3")
    paper_ready_2 <- identical(preset_choice, "paper_ready_2")
    paper_ready_3 <- identical(preset_choice, "paper_ready_3")
    size_scale <- suppressWarnings(as.numeric(input$export_quality))
    if (!is.finite(size_scale) || size_scale <= 0) {
      size_scale <- 1
    }

    sample_levels <- sample_order_state()
    if (!length(sample_levels)) {
      sample_levels <- as.character(sort(unique(summary$Sample)))
    }
    if (length(sample_levels)) {
      summary <- summary[as.character(summary$Sample) %in% sample_levels, , drop = FALSE]
      points <- points[as.character(points$Sample) %in% sample_levels, , drop = FALSE]
    }
    summary$Sample <- factor(as.character(summary$Sample), levels = sample_levels)
    points$Sample <- factor(as.character(points$Sample), levels = sample_levels)

    manual_colors <- parse_manual_colors(input$manual_palette)
    palette_values <- get_palette(input$palette_choice, length(sample_levels), manual_colors)
    color_values <- setNames(palette_values, sample_levels)
    fill_values <- color_values
    point_fill_values <- NULL
    if (paper_ready) {
      fill_values <- setNames(as.character(color_values), sample_levels)
      if (paper_ready_3) {
        fill_values[] <- "#ffffff"
      }
      point_base_colors <- if (paper_ready_3) color_values else fill_values
      point_fill_values <- setNames(darken_hex(point_base_colors, factor = 0.7), sample_levels)
      if (paper_ready_2) {
        point_fill_values[] <- "#ffffff"
      }
      summary$bar_fill <- fill_values[as.character(summary$Sample)]
      points$point_fill <- point_fill_values[as.character(points$Sample)]
    }

    axis_color <- if (paper_ready) "#000000" else "#222222"
    axis_line_width <- (if (paper_ready) 0.7 else 0.4) * size_scale
    bar_width <- if (paper_ready) 0.45 else 0.55
    bar_line_width <- (if (paper_ready) 0.4 else 0.7) * size_scale
    errorbar_line_width <- (if (paper_ready) 0.5 else 0.6) * size_scale
    point_size <- input$point_size * size_scale
    point_stroke <- 0.5 * size_scale
    jitter_width <- input$jitter_width
    x_expand <- input$x_expand
    facet_spacing <- input$facet_spacing
    aspect_ratio <- input$plot_aspect
    base_size <- input$base_font_size * size_scale
    angle <- if (isTRUE(input$rotate_x_labels)) input$x_label_angle else 0

    quote_name <- function(nm) {
      if (is.null(nm) || !nzchar(nm)) {
        return("``")
      }
      if (make.names(nm) != nm || grepl("[^A-Za-z0-9_.]", nm)) {
        paste0("`", nm, "`")
      } else {
        nm
      }
    }
    vec_to_r <- function(x) {
      if (is.null(x) || !length(x)) {
        return("character(0)")
      }
      vals <- paste(sprintf("\"%s\"", x), collapse = ", ")
      paste0("c(", vals, ")")
    }
    named_vec_to_r <- function(x) {
      if (is.null(x) || !length(x)) {
        return("character(0)")
      }
      parts <- mapply(
        function(nm, val) paste0(quote_name(nm), " = \"", val, "\""),
        names(x),
        x,
        SIMPLIFY = TRUE,
        USE.NAMES = FALSE
      )
      paste0("c(", paste(parts, collapse = ", "), ")")
    }

    summary_code <- capture.output(dput(summary))
    points_code <- capture.output(dput(points))

    theme_line <- switch(
      input$theme_choice,
      gray = "theme_gray",
      bw = "theme_bw",
      linedraw = "theme_linedraw",
      light = "theme_light",
      dark = "theme_dark",
      minimal = "theme_minimal",
      classic = "theme_classic",
      void = "theme_void",
      "theme_classic"
    )

    legend_position <- if (isTRUE(input$legend_on)) "right" else "none"
    dodge_width <- 0.55

    code <- c(
      "library(ggplot2)",
      "",
      "summary_df <-",
      summary_code,
      "",
      "points_df <-",
      points_code,
      "",
      "# ---- Tuning parameters ----",
      sprintf("base_size <- %.2f", base_size),
      sprintf("axis_line_width <- %.2f", axis_line_width),
      sprintf("bar_width <- %.2f", bar_width),
      sprintf("bar_line_width <- %.2f", bar_line_width),
      sprintf("errorbar_line_width <- %.2f", errorbar_line_width),
      sprintf("point_size <- %.2f", point_size),
      sprintf("point_stroke <- %.2f", point_stroke),
      sprintf("jitter_width <- %.2f", jitter_width),
      sprintf("dodge_width <- %.2f", dodge_width),
      sprintf("x_expand <- %.2f", x_expand),
      sprintf("facet_spacing <- %.2f", facet_spacing),
      sprintf("aspect_ratio <- %.2f", aspect_ratio),
      sprintf("axis_color <- \"%s\"", axis_color),
      sprintf("x_label_angle <- %d", angle),
      sprintf("legend_position <- \"%s\"", legend_position),
      "# ---------------------------",
      "",
      paste0("sample_levels <- ", vec_to_r(sample_levels)),
      "summary_df$Sample <- factor(summary_df$Sample, levels = sample_levels)",
      "points_df$Sample <- factor(points_df$Sample, levels = sample_levels)",
      paste0("color_values <- ", named_vec_to_r(color_values))
    )

    if (paper_ready) {
      code <- c(
        code,
        paste0("fill_values <- ", named_vec_to_r(fill_values)),
        "summary_df$bar_fill <- fill_values[as.character(summary_df$Sample)]",
        paste0("point_fill_values <- ", named_vec_to_r(point_fill_values)),
        "points_df$point_fill <- point_fill_values[as.character(points_df$Sample)]"
      )
    }

    code <- c(
      code,
      "",
      "p <- ggplot() +"
    )

    if (paper_ready) {
      code <- c(
        code,
        sprintf(
          "  geom_col(data = summary_df, aes(x = Sample, y = mean_fold, fill = bar_fill), width = bar_width, color = \"black\", linewidth = bar_line_width) +"
        ),
        sprintf(
          "  geom_errorbar(data = summary_df, aes(x = Sample, ymin = mean_fold - se_fold, ymax = mean_fold + se_fold), width = 0.2, linewidth = errorbar_line_width, color = \"black\") +"
        )
      )
    } else {
      code <- c(
        code,
        sprintf(
          "  geom_col(data = summary_df, aes(x = Sample, y = mean_fold, color = Sample), width = bar_width, fill = \"white\", linewidth = bar_line_width) +"
        ),
        sprintf(
          "  geom_errorbar(data = summary_df, aes(x = Sample, ymin = mean_fold - se_fold, ymax = mean_fold + se_fold, color = Sample), width = 0.2, linewidth = errorbar_line_width) +"
        )
      )
    }

    if (isTRUE(input$show_points)) {
      if (paper_ready) {
        if (paper_ready_2) {
          code <- c(
            code,
            sprintf(
              "  geom_point(data = points_df, aes(x = Sample, y = fold_change, shape = Sample), position = position_jitterdodge(jitter.width = jitter_width, dodge.width = dodge_width), size = point_size, stroke = point_stroke, color = \"black\", fill = \"white\") +"
            )
          )
        } else {
          code <- c(
            code,
            sprintf(
              "  geom_point(data = points_df, aes(x = Sample, y = fold_change, fill = point_fill), position = position_jitterdodge(jitter.width = jitter_width, dodge.width = dodge_width), shape = 21, size = point_size, stroke = point_stroke, color = \"black\") +"
            )
          )
        }
      } else {
        code <- c(
          code,
          sprintf(
            "  geom_jitter(data = points_df, aes(x = Sample, y = fold_change, color = Sample), width = jitter_width, height = 0, shape = 16, size = point_size) +"
          )
        )
      }
    }

    x_text_code <- if (isTRUE(input$show_sample_labels)) {
      sprintf("element_text(color = axis_color, angle = x_label_angle, hjust = %s)", if (angle > 0) "1" else "0.5")
    } else {
      "element_blank()"
    }
    x_ticks_code <- if (isTRUE(input$show_sample_labels)) {
      "element_line(color = axis_color, linewidth = axis_line_width)"
    } else {
      "element_blank()"
    }

    code <- c(
      code,
      sprintf("  facet_wrap(~Target, scales = \"free_y\") +"),
      "  scale_x_discrete(expand = expansion(mult = c(x_expand, x_expand))) +",
      "  coord_cartesian(clip = \"off\") +",
      "  labs(x = \"Sample\", y = expression(\"Fold Change (\"*2^{-Delta*Delta*Cq}*\")\")) +",
      sprintf("  %s(base_size = base_size, base_family = \"sans\") +", theme_line),
      sprintf("  theme(legend.position = legend_position, axis.line = element_line(color = axis_color, linewidth = axis_line_width), axis.ticks.x = %s, axis.ticks.y = element_line(color = axis_color, linewidth = axis_line_width), axis.text.x = %s, axis.text.y = element_text(color = axis_color), axis.title = element_text(color = axis_color), panel.spacing = grid::unit(facet_spacing, \"lines\"), aspect.ratio = aspect_ratio, strip.background = element_blank(), strip.text = element_text(face = \"bold\"))",
              x_ticks_code,
              x_text_code
      )
    )

    if (paper_ready) {
      code <- c(code, "p <- p + scale_fill_identity(guide = \"none\")")
    } else {
      code <- c(code, paste0("p <- p + scale_color_manual(values = ", named_vec_to_r(color_values), ", drop = FALSE)"))
    }

    if (paper_ready_2) {
      shape_values <- rep(c(21, 22, 23, 24, 25), length.out = length(sample_levels))
      code <- c(code, paste0("p <- p + scale_shape_manual(values = ", named_vec_to_r(setNames(shape_values, sample_levels)), ")"))
    }

    if (paper_ready) {
      code <- c(code, "p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())")
    }

    if (isTRUE(input$log_scale)) {
      code <- c(code, "p <- p + scale_y_continuous(trans = \"log2\")")
    } else if (identical(input$y_limit_mode, "manual")) {
      code <- c(code, sprintf("p <- p + scale_y_continuous(expand = c(0, 0), limits = c(%.4f, %.4f))", input$y_min, input$y_max))
    } else {
      code <- c(code, "p <- p + scale_y_continuous(expand = c(0, 0), limits = c(0, NA))")
    }

    code <- c(code, "", "p")
    paste(code, collapse = "\n")
  })

  output$summary_table <- renderDT({
    summary_state <- summary_all()
    validate(need(summary_state$ok, summary_state$message))
    summary_data <- summary_state$data
    validate(need(nrow(summary_data) > 0, "No summary rows are available to display."))
    datatable(
      summary_data,
      rownames = FALSE,
      options = list(pageLength = 10, scrollX = TRUE),
      extensions = "Buttons"
    )
  })

  output$plot_message <- renderUI({
    summary_state <- filtered_summary()
    if (!summary_state$ok) {
      return(tags$div(class = "text-danger", summary_state$message))
    }
    if (!nrow(summary_state$data)) {
      return(tags$div(class = "text-danger", "Select at least one target to display."))
    }
    NULL
  })

  output$summary_message <- renderUI({
    summary_state <- summary_all()
    if (!summary_state$ok) {
      return(tags$div(class = "text-danger", summary_state$message))
    }
    if (!nrow(summary_state$data)) {
      return(tags$div(class = "text-danger", "No summary rows are available to display."))
    }
    NULL
  })

  output$data_error <- renderUI({
    msg <- data_error_msg()
    if (is.null(msg) || !nzchar(msg)) {
      return(NULL)
    }
    tags$div(class = "text-danger", msg)
  })

  output$manual_cq_table <- renderRHandsontable({
    table <- rhandsontable(
      manual_cq(),
      rowHeaders = FALSE,
      stretchH = "all",
      height = 350
    )
    hot_col(table, "Well", readOnly = TRUE)
  })

  output$manual_plate_table <- renderRHandsontable({
    table <- rhandsontable(
      manual_plate(),
      rowHeaders = FALSE,
      stretchH = "all",
      height = 320
    )
    table <- hot_col(table, "Row", readOnly = TRUE)
    hot_col(table, "Type", readOnly = TRUE)
  })

  output$full_table <- renderDT({
    datatable(
      processed_wells(),
      rownames = FALSE,
      options = list(pageLength = 12, scrollX = TRUE)
    )
  })

  output$download_summary <- downloadHandler(
    filename = function() {
      paste0("qpcr_summary_", Sys.Date(), ".csv")
    },
    content = function(file) {
      summary_state <- summary_all()
      req(summary_state$ok)
      write.csv(summary_state$data, file, row.names = FALSE, na = "")
    }
  )
}

shinyApp(ui, server)
