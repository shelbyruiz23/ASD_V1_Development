rm(list=ls()); gc()
options(stringsAsFactors = F)

# Describe identified co-regulated modules
####################################################################
####################################################################

# load libraries
require(data.table)           # fast version of read.table
require(readxl)               # used for reading Excel format files
#require(viridis)              # colorblind friendly color palette
require(matrixStats)          # to get statistics on matrices
library(dplyr)
library(RColorBrewer)
library(tidyverse)
library(ggplot2)
library(patchwork)


####################################################################
# load Supplementary Table 22, Sheet "module_proteins" and separate into the 4 datasets
# load the representative ORA term as determined by the script module_representative_term

#load("Resubmission/WorkData/modules.RData")
#load("Resubmission/WorkData/modules_ORA_Representative.RData")
# ============================================================
#  4-Panel Module Summary — 2×2 grid layout
#
#  Top-left:    HOM — Young      Top-right:   SYN — Young
#  Bottom-left: HOM — Old        Bottom-right: SYN — Old
#
#  HOM panels: module | n | q.DX | q.AGE | term
#  SYN panels: module | n | q.DX | q.AGE | q.DX.AGE | term
#
#  All circle sizes on same global scale.
#  Rows spaced by circle radius so top (large) rows never overlap.
#
#  install.packages(c("tidyverse", "patchwork"))
# ============================================================



# ============================================================
# A.  DATA PREPARATION
# ============================================================

flatten_modules <- function(module_list, age_label) {
  bind_rows(lapply(names(module_list), function(mod_name) {
    df        <- module_list[[mod_name]]
    df$module <- mod_name
    df$age    <- age_label
    df
  }))
}

clean_term <- function(x) {
  ifelse(is.na(x), "—",
         x |>
           str_remove("^(GOBP_|GOCC_|REACTOME_|KEGG_|GOMF_)") |>
           str_replace_all("_", " ") |>
           str_to_sentence() |>
           str_trunc(40, ellipsis = "…")
  )
}

hom_yng_flat <- flatten_modules(HOM.YNG, "young")
hom_old_flat <- flatten_modules(HOM.OLD, "old")
syn_yng_flat <- flatten_modules(SYN.YNG, "young")
syn_old_flat <- flatten_modules(SYN.OLD, "old")

get_enrich <- function(rep_table) {
  rep_table |>
    rownames_to_column("module") |>
    mutate(term = clean_term(GSterm)) |>
    select(module, term)
}

enrich_hom_yng <- get_enrich(HOM.YNG_Representative)
enrich_hom_old <- get_enrich(HOM.OLD_Representative)
enrich_syn_yng <- get_enrich(SYN.YNG_Representative)
enrich_syn_old <- get_enrich(SYN.OLD_Representative)

# HOM: no q.DX.AGE column
build_hom <- function(flat_df, enrich_df) {
  flat_df |>
    group_by(module) |>
    summarise(
      n        = n(),
      hex_code = first(hex_code),
      n_DE_dx  = sum(q.DX  < 0.10, na.rm = TRUE),
      n_DE_age = sum(q.AGE < 0.05, na.rm = TRUE),
      .groups  = "drop"
    ) |>
    left_join(enrich_df, by = "module") |>
    replace_na(list(term = "—")) |>
    arrange(desc(n)) |>
    mutate(row = row_number())
}

# SYN: includes q.DX.AGE interaction column
build_syn <- function(flat_df, enrich_df) {
  flat_df |>
    group_by(module) |>
    summarise(
      n           = n(),
      hex_code    = first(hex_code),
      n_DE_dx     = sum(q.DX     < 0.10, na.rm = TRUE),
      n_DE_age    = sum(q.AGE    < 0.05, na.rm = TRUE),
      n_DE_dx_age = sum(q.DX.AGE < 0.10, na.rm = TRUE),
      .groups     = "drop"
    ) |>
    left_join(enrich_df, by = "module") |>
    replace_na(list(term = "—")) |>
    arrange(desc(n)) |>
    mutate(row = row_number())
}

hom_yng <- build_hom(hom_yng_flat, enrich_hom_yng)
hom_old <- build_hom(hom_old_flat, enrich_hom_old)
syn_yng <- build_syn(syn_yng_flat, enrich_syn_yng)
syn_old <- build_syn(syn_old_flat, enrich_syn_old)

####################################################################
# ============================================================
#  Fisher's Exact Test — OR and p-value per module per DE test
#  Background population sizes are fixed from the full dataset
# ============================================================

# ── Fixed background totals ───────────────────────────────────
HOM_TOTAL     <- 4727
HOM_DE_DX     <- 230
HOM_DE_AGE    <- 1261

SYN_TOTAL     <- 4287
SYN_DE_DX     <- 304
SYN_DE_AGE    <- 912
SYN_DE_DX_AGE <- 305

# ── Fisher function: returns OR and p-value ───────────────────
run_fisher <- function(n_de_in, n_total_in,
                       n_de_total, n_total) {
  a <- n_de_in
  b <- n_de_total  - n_de_in          # DE outside module
  c <- n_total_in  - n_de_in          # not DE inside module
  d <- (n_total - n_total_in) - b     # not DE outside module
  
  # Guard: negative cells = impossible table, return NA
  if (any(c(a, b, c, d) < 0)) {
    return(tibble(OR = NA_real_, p_fisher = NA_real_))
  }
  
  mat    <- matrix(c(a, c, b, d), nrow = 2,
                   dimnames = list(c("DE", "Not DE"),
                                   c("In module", "Outside")))
  result <- fisher.test(mat, alternative = "greater")
  
  tibble(
    OR       = result$estimate,   # odds ratio
    p_fisher = result$p.value
  )
}

# ── Apply across all modules in a summary data frame ─────────
add_fisher_tests <- function(df, dataset_type) {
  
  if (dataset_type == "HOM") {
    total     <- HOM_TOTAL
    de_dx     <- HOM_DE_DX
    de_age    <- HOM_DE_AGE
    de_dx_age <- NA_integer_
  } else {
    total     <- SYN_TOTAL
    de_dx     <- SYN_DE_DX
    de_age    <- SYN_DE_AGE
    de_dx_age <- SYN_DE_DX_AGE
  }
  
  results <- pmap_dfr(df, function(module, n, n_DE_dx, n_DE_age,
                                   n_DE_dx_age = 0, ...) {
    
    f_dx  <- run_fisher(n_DE_dx,  n, de_dx,  total)
    f_age <- run_fisher(n_DE_age, n, de_age, total)
    
    out <- tibble(
      module       = module,
      n            = n,
      # DX
      n_DE_dx      = n_DE_dx,
      OR_dx        = f_dx$OR,
      p_dx         = f_dx$p_fisher,
      sig_dx       = !is.na(f_dx$p_fisher) & f_dx$p_fisher < 0.05,
      # AGE
      n_DE_age     = n_DE_age,
      OR_age       = f_age$OR,
      p_age        = f_age$p_fisher,
      sig_age      = !is.na(f_age$p_fisher) & f_age$p_fisher < 0.05
    )
    
    # DX:AGE only for SYN
    if (dataset_type == "SYN") {
      f_dxa <- run_fisher(n_DE_dx_age, n, de_dx_age, total)
      out <- out |>
        mutate(
          n_DE_dx_age  = n_DE_dx_age,
          OR_dx_age    = f_dxa$OR,
          p_dx_age     = f_dxa$p_fisher,
          sig_dx_age   = !is.na(f_dxa$p_fisher) & f_dxa$p_fisher < 0.05
        )
    }
    
    out
  })
  
  results
}

# ── Run for all four datasets ─────────────────────────────────
fisher_hom_yng <- add_fisher_tests(hom_yng, "HOM") |> mutate(dataset = "HOM_YNG")
fisher_hom_old <- add_fisher_tests(hom_old, "HOM") |> mutate(dataset = "HOM_OLD")
fisher_syn_yng <- add_fisher_tests(syn_yng, "SYN") |> mutate(dataset = "SYN_YNG")
fisher_syn_old <- add_fisher_tests(syn_old, "SYN") |> mutate(dataset = "SYN_OLD")

# ── Combine into one table ────────────────────────────────────
fisher_results <- bind_rows(
  fisher_hom_yng,
  fisher_hom_old,
  fisher_syn_yng |> select(dataset, module, n,
                           n_DE_dx, OR_dx, p_dx, sig_dx,
                           n_DE_age, OR_age, p_age, sig_age,
                           n_DE_dx_age, OR_dx_age, p_dx_age, sig_dx_age),
  fisher_syn_old |> select(dataset, module, n,
                           n_DE_dx, OR_dx, p_dx, sig_dx,
                           n_DE_age, OR_age, p_age, sig_age,
                           n_DE_dx_age, OR_dx_age, p_dx_age, sig_dx_age)
)

# ── Preview ───────────────────────────────────────────────────
fisher_results |>
  select(dataset, module, n,
         OR_dx, p_dx, sig_dx,
         OR_age, p_age, sig_age) |>
  print(n = 20)

# ── Save ──────────────────────────────────────────────────────
#write.csv(fisher_results, "module_fisher_results.csv", row.names = FALSE)

# ============================================================
# B.  GLOBAL SCALES  (shared across all four panels)
# ============================================================

all_data       <- bind_rows(
  hom_yng |> mutate(n_DE_dx_age = 0L),
  hom_old |> mutate(n_DE_dx_age = 0L),
  syn_yng,
  syn_old
)

GLOBAL_MAX_N       <- max(all_data$n)
GLOBAL_MAX_DX      <- max(all_data$n_DE_dx,     1)
GLOBAL_MAX_AGE     <- max(all_data$n_DE_age,    1)
GLOBAL_MAX_DX_AGE  <- max(all_data$n_DE_dx_age, 1)

# Circle radii in data units
R_MOD_MAX <- 0.42   # largest module circle radius
R_DE_MAX  <- 0.28   # largest DE dot radius

r_mod <- function(n) sqrt(n / GLOBAL_MAX_N)  * R_MOD_MAX
r_de  <- function(n, mx) ifelse(n > 0,
                                pmax(0.05, sqrt(n / mx) * R_DE_MAX),
                                NA_real_)

# Convert data-unit radius → geom_point size (diameter in mm)
# Panel width = half of 6.5 in = 3.25 in; data x range ≈ 12 units
# 1 data unit ≈ 3.25/12 * 25.4 ≈ 6.88 mm
DATA_TO_MM <- (3.25 / 12) * 25.4
r_to_size  <- function(r) r * DATA_TO_MM * 2


# ============================================================
# C.  ROW Y POSITIONS  (shared within each paired row)
#
#  Within each horizontal pair (YNG pair, OLD pair) the y
#  positions are computed from the UNION of both panels so
#  row i in HOM.YNG sits at exactly the same y as row i in
#  SYN.YNG.  The taller panel simply extends further down.
# ============================================================

GAP <- 0.18   # minimum gap between circle edges (data units)

# Compute shared y positions for a pair of datasets.
# Uses the LARGER radius of the two panels at each row to set
# spacing, so neither panel overflows its neighbour's rows.
shared_y_pos <- function(df_left, df_right) {
  # Max rows in either panel
  n_shared <- max(nrow(df_left), nrow(df_right))
  
  # Pad the shorter one with dummy rows so we can take row-wise max(r)
  pad <- function(df) {
    if (nrow(df) < n_shared) {
      extra <- tibble(n = rep(0, n_shared - nrow(df)))
      df    <- bind_rows(df, extra)
    }
    df |> mutate(r = r_mod(pmax(n, 0)))
  }
  dl <- pad(df_left  |> mutate(r = r_mod(n)))
  dr <- pad(df_right |> mutate(r = r_mod(n)))
  
  # Row spacing driven by the larger circle at each row
  r_max <- pmax(dl$r, dr$r)
  
  y_pos <- numeric(n_shared)
  y_pos[1] <- 0
  if (n_shared > 1) {
    for (i in 2:n_shared) {
      y_pos[i] <- y_pos[i-1] - (r_max[i-1] + r_max[i] + GAP)
    }
  }
  y_pos
}

attach_y <- function(df, y_vec) {
  df |>
    mutate(r = r_mod(n)) |>
    mutate(y_pos = y_vec[row_number()])
}

# Compute shared y positions for each paired row
yng_y <- shared_y_pos(hom_yng |> mutate(r = r_mod(n)),
                      syn_yng |> mutate(r = r_mod(n)))
old_y <- shared_y_pos(hom_old |> mutate(r = r_mod(n)),
                      syn_old |> mutate(r = r_mod(n)))

hom_yng <- attach_y(hom_yng, yng_y)
syn_yng <- attach_y(syn_yng, yng_y)
hom_old <- attach_y(hom_old, old_y)
syn_old <- attach_y(syn_old, old_y)


# ============================================================
# D.  COLUMN X POSITIONS
#     Two layouts: HOM (no DX.AGE) and SYN (with DX.AGE)
#     Both panels are ~3.25 in wide (half page width)
# ============================================================

CHAR_W   <- 0.17   # estimated data-unit width per character
NAME_PAD <- 0.12

# Longest name across all datasets
max_name_chars <- max(
  sapply(list(hom_yng, hom_old, syn_yng, syn_old),
         function(d) max(nchar(d$module)))
)
name_col_w <- max_name_chars * CHAR_W + NAME_PAD

X_MOD  <- 0
X_NAME <- 0.42       # right edge of circle + small gap

# HOM layout
HOM_X_N    <- X_NAME + name_col_w
HOM_X_DX   <- HOM_X_N   + 0.50
HOM_X_AGE  <- HOM_X_DX  + 0.62
HOM_X_TERM <- HOM_X_AGE + 0.55
HOM_X_MAX  <- HOM_X_TERM + 5.8

# SYN layout (one extra column)
SYN_X_N      <- X_NAME + name_col_w
SYN_X_DX     <- SYN_X_N     + 0.50
SYN_X_AGE    <- SYN_X_DX    + 0.62
SYN_X_DX_AGE <- SYN_X_AGE   + 0.62
SYN_X_TERM   <- SYN_X_DX_AGE + 0.55
SYN_X_MAX    <- SYN_X_TERM + 5.2


# ============================================================
# E.  PANEL BUILDER
# ============================================================

make_panel <- function(df,
                       panel_title,
                       is_syn          = FALSE,
                       show_col_headers = TRUE) {
  
  # Unpack x positions for this layout
  if (is_syn) {
    xN <- SYN_X_N; xDX <- SYN_X_DX; xAGE <- SYN_X_AGE
    xDXAGE <- SYN_X_DX_AGE; xTERM <- SYN_X_TERM; xMAX <- SYN_X_MAX
  } else {
    xN <- HOM_X_N; xDX <- HOM_X_DX; xAGE <- HOM_X_AGE
    xDXAGE <- NA;  xTERM <- HOM_X_TERM; xMAX <- HOM_X_MAX
  }
  
  df <- df |>
    mutate(
      r_dx     = r_de(n_DE_dx,  GLOBAL_MAX_DX),
      r_age    = r_de(n_DE_age, GLOBAL_MAX_AGE),
      r_dx_age = if (is_syn) r_de(n_DE_dx_age, GLOBAL_MAX_DX_AGE)
      else NA_real_
    )
  
  y_top    <- df$y_pos[1] + df$r[1] + 0.06
  y_bottom <- min(df$y_pos) - max(df$r) - 0.06
  
  # ── dot helper layers (returns a list of geoms) ───────────
  de_layers <- function(r_col, x_col, fill_col, text_col,
                        count_col, threshold_size = 2.2) {
    d_has  <- df |> filter(!is.na(.data[[r_col]]))
    d_none <- df |> filter(is.na(.data[[r_col]]))
    d_lbl  <- df |> filter(!is.na(.data[[r_col]]),
                           r_to_size(.data[[r_col]]) >= threshold_size)
    list(
      geom_point(
        data = d_has,
        aes(x = .env$x_col, y = y_pos,
            size = I(r_to_size(.data[[r_col]]))),
        color = fill_col, alpha = 0.85, shape = 16
      ),
      geom_text(
        data = d_lbl,
        aes(x = .env$x_col, y = y_pos,
            label = .data[[count_col]]),
        size = 1.25, color = "white", fontface = "bold"
      ),
      geom_text(
        data = d_none,
        aes(x = .env$x_col, y = y_pos),
        label = "—", size = 1.35, color = "grey78"
      )
    )
  }
  
  p <- ggplot(df) +
    
    # Alternating row shading
    geom_rect(
      aes(xmin = X_MOD - 0.38, xmax = xMAX,
          ymin = y_pos - r - GAP/2,
          ymax = y_pos + r + GAP/2,
          fill = row %% 2 == 0),
      show.legend = FALSE
    ) +
    scale_fill_manual(values = c("TRUE" = "#f7f7f7", "FALSE" = "white")) +
    
    # Module circle
    geom_point(
      aes(x = X_MOD, y = y_pos, size = I(r_to_size(r))),
      color = df$hex_code, alpha = 0.88, shape = 16
    ) +
    
    # Module name
    geom_text(
      aes(x = X_NAME, y = y_pos, label = module),
      hjust = 0, size = 1.75, fontface = "bold", color = "grey12"
    ) +
    
    # n proteins
    geom_text(
      aes(x = xN, y = y_pos, label = n),
      hjust = 1, size = 1.65, color = "grey45"
    ) +
    
    # q.DX dots
    de_layers("r_dx",  xDX,  "#D4537E", "white", "n_DE_dx") +
    
    # q.AGE dots
    de_layers("r_age", xAGE, "#378ADD", "white", "n_DE_age") +
    
    # q.DX.AGE dots (SYN only)
    {if (is_syn)
      de_layers("r_dx_age", xDXAGE, "#1D9E75", "white", "n_DE_dx_age")
      else list()} +
    
    # Enrichment term
    geom_text(
      aes(x = xTERM, y = y_pos, label = term),
      hjust = 0, size = 1.55, color = "grey38", fontface = "italic"
    ) +
    
    # Panel title
    annotate("text",
             x = X_MOD - 0.3,
             y = y_top + (if (show_col_headers) 0.60 else 0.22),
             label = panel_title,
             hjust = 0, vjust = 1,
             size = 2.1, fontface = "bold", color = "grey18") +
    
    # Column headers
    {if (show_col_headers) {
      hy1 <- y_top + 0.44; hy2 <- y_top + 0.22; div_y <- y_top + 0.09
      hdr <- list(
        annotate("text", x = X_NAME, y = hy1, label = "Module",
                 hjust = 0, size = 1.65, fontface = "bold", color = "grey35"),
        annotate("text", x = xN,     y = hy1, label = "n",
                 hjust = 1, size = 1.65, fontface = "bold", color = "grey35"),
        annotate("text", x = xDX,    y = hy1, label = "q.DX",
                 hjust = 0.5, size = 1.65, fontface = "bold", color = "#993556"),
        annotate("text", x = xDX,    y = hy2, label = "< 0.1",
                 hjust = 0.5, size = 1.30, color = "#993556"),
        annotate("text", x = xAGE,   y = hy1, label = "q.AGE",
                 hjust = 0.5, size = 1.65, fontface = "bold", color = "#185FA5"),
        annotate("text", x = xAGE,   y = hy2, label = "< 0.05",
                 hjust = 0.5, size = 1.30, color = "#185FA5"),
        annotate("text", x = xTERM,  y = hy1, label = "Enrichment term",
                 hjust = 0,   size = 1.65, fontface = "bold", color = "grey35"),
        annotate("segment",
                 x = X_MOD - 0.38, xend = xMAX,
                 y = div_y, yend = div_y,
                 color = "grey45", linewidth = 0.28)
      )
      if (is_syn) {
        hdr <- c(hdr, list(
          annotate("text", x = xDXAGE, y = hy1, label = "q.DX.AGE",
                   hjust = 0.5, size = 1.65, fontface = "bold", color = "#0F6E56"),
          annotate("text", x = xDXAGE, y = hy2, label = "< 0.1",
                   hjust = 0.5, size = 1.30, color = "#0F6E56")
        ))
      }
      hdr
    } else {
      list(annotate("segment",
                    x = X_MOD - 0.38, xend = xMAX,
                    y = y_top + 0.08, yend = y_top + 0.08,
                    color = "grey78", linewidth = 0.18))
    }} +
    
    coord_cartesian(
      xlim = c(X_MOD - 0.4, xMAX),
      ylim = c(y_bottom,
               y_top + (if (show_col_headers) 0.72 else 0.30)),
      clip = "off"
    ) +
    theme_void(base_size = 6) +
    theme(
      plot.background = element_rect(fill = "white", colour = NA),
      plot.margin     = margin(1, 3, 1, 3)
    )
  
  p
}


# ============================================================
# F.  LEGEND PANEL (spans full width below 2×2 grid)
# ============================================================

make_legend <- function() {
  ns     <- c(GLOBAL_MAX_N,
              round(GLOBAL_MAX_N * 0.5),
              round(GLOBAL_MAX_N * 0.25),
              round(GLOBAL_MAX_N * 0.1))
  ns     <- ns[ns > 0]
  leg_df <- tibble(
    x   = cumsum(c(0.6, rep(1.4, length(ns)-1))),
    y   = 0,
    n   = ns,
    r   = r_mod(ns),
    lbl = as.character(ns)
  )
  
  # DE dot examples at half max
  de_x_start <- leg_df$x[length(ns)] + 1.8
  
  ggplot() +
    geom_point(data = leg_df,
               aes(x = x, y = y, size = I(r_to_size(r))),
               color = "grey55", alpha = 0.7, shape = 16) +
    geom_text(data = leg_df,
              aes(x = x, y = -0.30, label = lbl),
              size = 1.4, color = "grey50", hjust = 0.5) +
    annotate("text", x = 0.1, y = 0,
             label = "n =", hjust = 1, size = 1.55, color = "grey40") +
    
    annotate("point", x = de_x_start, y = 0.05,
             size = r_to_size(R_DE_MAX * 0.6),
             color = "#D4537E", alpha = 0.85, shape = 16) +
    annotate("text", x = de_x_start + 0.45, y = 0.05,
             label = "q.DX < 0.1", hjust = 0,
             size = 1.55, color = "#993556") +
    
    annotate("point", x = de_x_start + 3.0, y = 0.05,
             size = r_to_size(R_DE_MAX * 0.6),
             color = "#378ADD", alpha = 0.85, shape = 16) +
    annotate("text", x = de_x_start + 3.45, y = 0.05,
             label = "q.AGE < 0.05", hjust = 0,
             size = 1.55, color = "#185FA5") +
    
    annotate("point", x = de_x_start + 6.2, y = 0.05,
             size = r_to_size(R_DE_MAX * 0.6),
             color = "#1D9E75", alpha = 0.85, shape = 16) +
    annotate("text", x = de_x_start + 6.65, y = 0.05,
             label = "q.DX.AGE < 0.1 (SYN only)", hjust = 0,
             size = 1.55, color = "#0F6E56") +
    
    annotate("text", x = 0.1, y = -0.5,
             label = "All circle sizes on same scale across panels",
             hjust = 0, size = 1.4, color = "grey55",
             fontface = "italic") +
    
    coord_cartesian(xlim = c(0, SYN_X_MAX),
                    ylim = c(-0.65, 0.38), clip = "off") +
    theme_void() +
    theme(plot.background = element_rect(fill = "white", colour = NA),
          plot.margin     = margin(3, 3, 5, 3))
}


# ============================================================
# G.  ASSEMBLE 2×2 GRID + LEGEND
# ============================================================

p_hom_yng <- make_panel(hom_yng, "HOM — Young",
                        is_syn = FALSE, show_col_headers = TRUE)
p_hom_old <- make_panel(hom_old, "HOM — Old",
                        is_syn = FALSE, show_col_headers = FALSE)
p_syn_yng <- make_panel(syn_yng, "SYN — Young",
                        is_syn = TRUE,  show_col_headers = TRUE)
p_syn_old <- make_panel(syn_old, "SYN — Old",
                        is_syn = TRUE,  show_col_headers = FALSE)
p_legend  <- make_legend()

# Height of each panel row = span of the shared y positions
# (both panels in a pair share the same y grid, so the taller
#  one simply has more rows; use the taller one's y extent)
panel_height <- function(df) {
  abs(min(df$y_pos) - df$r[nrow(df)]) +
    (df$y_pos[1]    + df$r[1])
}

h_yng <- max(panel_height(hom_yng), panel_height(syn_yng))
h_old <- max(panel_height(hom_old), panel_height(syn_old))

# 2x2 grid: left = HOM, right = SYN
# Panels in the same row share y positions but may have different
# numbers of modules — the shorter panel clips at its last module
# while the taller one extends further down naturally.
grid_top <- (p_hom_yng | p_syn_yng) +
  plot_layout(widths = c(1, 1))

grid_bot <- (p_hom_old | p_syn_old) +
  plot_layout(widths = c(1, 1))

full_fig <- (grid_top / grid_bot / p_legend) +
  plot_layout(heights = c(h_yng, h_old, 0.7))

# ── Save ─────────────────────────────────────────────────────
FIG_W <- 6.5     # full printable width (in)
FIG_H <- 4.75    # half of 9.5 in printable height

ggsave("module_summary_4panel.pdf",
       plot = full_fig, device = cairo_pdf,
       width = FIG_W, height = FIG_H)

ggsave("module_summary_4panel.png",
       plot = full_fig, device = "png",
       width = FIG_W, height = FIG_H,
       dpi = 600, bg = "white")

