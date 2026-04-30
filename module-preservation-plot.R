rm(list=ls()); gc()
options(stringsAsFactors = F)

# Plot a heatmap plot of the module preservation status between all identified modules
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


####################################################################
# load Supplementary Table 22, Sheet "module_map" 
# load Supplementary Table 22, Sheet "module_proteins" and separate into 4 datasets
# load Supplementary Table 22, Sheet "module_comparisons" and separate into 6 comparisons

#RESULTS=readRDS(file="Resubmission/Results/ASD-DEV-clustering-comparisons-2025-07-30.RDS")
#load("Resubmission/WorkData/module-renaming.RData")
#load("Resubmission/WorkData/module_colors.RData")

big_list=list(cluster.HOM.OLD, cluster.HOM.YNG, cluster.SYN.OLD, cluster.SYN.YNG)


# ============================================================
#  Build comparison-aware module_map
# ============================================================

# Step 1: Normalise dataset names
module_map <- module_map |>
  mutate(dataset = str_replace_all(dataset, "\\.", "_"))

# Step 2: Define all 6 comparisons and which position each dataset appears in
comparisons <- tibble(
  comparison = names(RESULTS),
  ds1 = str_split(comparison, ":", simplify = TRUE)[,1],
  ds2 = str_split(comparison, ":", simplify = TRUE)[,2]
)

# Step 3: For each comparison, create rows for ds1 (prefix=1) and ds2 (prefix=2)
# The module_num in RESULTS uses the prefix to encode position:
# e.g. in HOM_YNG:HOM_OLD, HOM_YNG modules appear as 1.x, HOM_OLD as 2.x
# The extension (.1, .10, .83 etc.) comes from the original module_num decimal

module_map_comp <- bind_rows(
  # ds1 gets prefix 1
  comparisons |>
    left_join(
      module_map |>
        mutate(
          extension  = str_extract(as.character(module_num), "\\..*"),
          comp_num   = paste0("1", extension)
        ) |>
        select(dataset, module, hex_code, module_num, extension, comp_num),
      by = c("ds1" = "dataset")
    ) |>
    mutate(position = "cl1", comp_module_num = comp_num),
  
  # ds2 gets prefix 2
  comparisons |>
    left_join(
      module_map |>
        mutate(
          extension  = str_extract(as.character(module_num), "\\..*"),
          comp_num   = paste0("2", extension)
        ) |>
        select(dataset, module, hex_code, module_num, extension, comp_num),
      by = c("ds2" = "dataset")
    ) |>
    mutate(position = "cl2", comp_module_num = comp_num)
) |>
  select(comparison, position, dataset = ds1,
         module, hex_code, module_num, comp_module_num) |>
  mutate(dataset = if_else(position == "cl2",
                           str_extract(comparison, "[^:]+$"),
                           str_extract(comparison, "^[^:]+")))

# Verify — check one comparison
module_map_comp |>
  filter(comparison == names(RESULTS)[1]) |>
  select(comparison, position, dataset, module, comp_module_num) |>
  print(n = 10)

# Step 4: Remap RESULTS using comparison-specific lookup
remap_results <- function(cmp) {
  df <- RESULTS[[cmp]]$results
  
  lk <- module_map_comp |>
    filter(comparison == cmp) |>
    select(position, comp_module_num, module, hex_code)
  
  lk1 <- lk |> filter(position == "cl1")
  lk2 <- lk |> filter(position == "cl2")
  
  df |>
    mutate(
      cl1 = coalesce(
        lk1$module[match(as.character(cl1), lk1$comp_module_num)],
        as.character(cl1)
      ),
      cl2 = coalesce(
        lk2$module[match(as.character(cl2), lk2$comp_module_num)],
        as.character(cl2)
      )
    )
}

RESULTS_NAMED <- RESULTS
for (cmp in names(RESULTS)) {
  RESULTS_NAMED[[cmp]]$results <- remap_results(cmp)
}

# Verify all comparisons
for (cmp in names(RESULTS_NAMED)) {
  df        <- RESULTS_NAMED[[cmp]]$results
  still_num <- sum(str_detect(as.character(df$cl1), "^[0-9]")) +
    sum(str_detect(as.character(df$cl2), "^[0-9]"))
  message(cmp, ": unmapped=", still_num,
          " | cl1: ", paste(head(df$cl1, 3), collapse=", "),
          " | cl2: ", paste(head(df$cl2, 3), collapse=", "))
}

# ============================================================
#  Module Preservation Heatmap
#
#  Assumes RESULTS_NAMED already exists in your environment
#  with cl1/cl2 remapped to module names like "coral4".
#  If not, run the remapping block at the top first.
#
#  install.packages(c("tidyverse","ggplot2","scales","patchwork"))
# ============================================================

library(tidyverse)
library(ggplot2)
library(scales)
library(patchwork)


# ============================================================
# STEP 1.  REMAP RESULTS TO MODULE NAMES  (run if needed)
#          Skip if RESULTS_NAMED already exists
# ============================================================

if (!exists("RESULTS_NAMED")) {
  
  module_map <- module_map |>
    mutate(dataset = str_replace_all(dataset, "\\.", "_"))
  
  comparisons <- tibble(
    comparison = names(RESULTS),
    ds1 = str_split(comparison, ":", simplify = TRUE)[,1],
    ds2 = str_split(comparison, ":", simplify = TRUE)[,2]
  )
  
  module_map_comp <- bind_rows(
    comparisons |>
      left_join(
        module_map |>
          mutate(
            extension  = str_extract(as.character(module_num), "\\..*"),
            comp_num   = paste0("1", extension)
          ) |>
          select(dataset, module, hex_code, module_num, extension, comp_num),
        by = c("ds1" = "dataset")
      ) |>
      mutate(position = "cl1", comp_module_num = comp_num,
             dataset = ds1),
    
    comparisons |>
      left_join(
        module_map |>
          mutate(
            extension  = str_extract(as.character(module_num), "\\..*"),
            comp_num   = paste0("2", extension)
          ) |>
          select(dataset, module, hex_code, module_num, extension, comp_num),
        by = c("ds2" = "dataset")
      ) |>
      mutate(position = "cl2", comp_module_num = comp_num,
             dataset = ds2)
  ) |>
    select(comparison, position, dataset, module, hex_code, comp_module_num)
  
  remap_results <- function(cmp) {
    df  <- RESULTS[[cmp]]$results
    lk  <- module_map_comp |> filter(comparison == cmp)
    lk1 <- lk |> filter(position == "cl1")
    lk2 <- lk |> filter(position == "cl2")
    df |>
      mutate(
        cl1 = coalesce(
          lk1$module[match(as.character(cl1), lk1$comp_module_num)],
          as.character(cl1)
        ),
        cl2 = coalesce(
          lk2$module[match(as.character(cl2), lk2$comp_module_num)],
          as.character(cl2)
        )
      )
  }
  
  RESULTS_NAMED <- RESULTS
  for (cmp in names(RESULTS)) {
    RESULTS_NAMED[[cmp]]$results <- remap_results(cmp)
  }
  
  # Verify
  for (cmp in names(RESULTS_NAMED)) {
    df        <- RESULTS_NAMED[[cmp]]$results
    still_num <- sum(str_detect(as.character(df$cl1), "^[0-9]")) +
      sum(str_detect(as.character(df$cl2), "^[0-9]"))
    message(cmp, ": unmapped=", still_num,
            " | cl1: ", paste(head(df$cl1, 2), collapse=", "),
            " | cl2: ", paste(head(df$cl2, 2), collapse=", "))
  }
}


# ============================================================
# STEP 2.  BUILD HEX COLOUR LOOKUP
#          module name → hex code, per dataset
# ============================================================

hex_lookup <- module_map |>
  mutate(dataset = str_replace_all(dataset, "\\.", "_")) |>
  select(dataset, module, hex_code) |>
  group_by(dataset, module) |>
  slice(1) |>
  ungroup()

get_hex <- function(mod_names, ds_name) {
  lk <- hex_lookup |> filter(dataset == ds_name)
  coalesce(lk$hex_code[match(mod_names, lk$module)], "#888888")
}


# ============================================================
# STEP 3.  FLATTEN TO LONG DATA FRAME
# ============================================================

comparison_names <- names(RESULTS_NAMED)

plong <- map_dfr(comparison_names, function(cmp) {
  df    <- RESULTS_NAMED[[cmp]]$results
  parts <- str_split(cmp, ":", simplify = TRUE)
  ds1   <- parts[1]
  ds2   <- parts[2]
  
  tibble(
    comparison = cmp,
    ds1        = ds1,
    ds2        = ds2,
    mod1       = as.character(df$cl1),
    mod2       = as.character(df$cl2),
    hex1       = get_hex(as.character(df$cl1), ds1),
    hex2       = get_hex(as.character(df$cl2), ds2),
    n1         = df$n.prt1,
    n2         = df$n.prt2,
    n_shared   = df$n.prt12,
    p          = df$p
  )
})


# ============================================================
# STEP 4.  COLOUR SCALE  (p = 0 → darkest colour)
# ============================================================

plong <- plong |>
  mutate(
    is_zero    = p == 0,
    log10p_raw = if_else(p == 0, NA_real_, -log10(p))
  )

scale_max <- quantile(plong$log10p_raw, 0.99, na.rm = TRUE)

plong <- plong |>
  mutate(
    log10p_plot = pmin(log10p_raw, scale_max, na.rm = TRUE),
    log10p_plot = if_else(is_zero, scale_max * 1.05, log10p_plot)
  )


# ============================================================
# STEP 5.  MODULE SIZE ORDERING
# ============================================================

sizes <- bind_rows(
  plong |> select(ds = ds1, mod = mod1, n = n1),
  plong |> select(ds = ds2, mod = mod2, n = n2)
) |>
  group_by(ds, mod) |>
  summarise(n = max(n), .groups = "drop")

order_by_size <- function(mod_vec, ds_name) {
  sz  <- sizes |> filter(ds == ds_name)
  ord <- sz$mod[order(-sz$n)]
  ord[ord %in% mod_vec]
}


# ============================================================
# STEP 6.  HEATMAP FUNCTION
# ============================================================

FILL_MID  <- "#5DCAA5"
FILL_HIGH <- "#04342C"
FILL_NA   <- "#E0E0E0"

make_heatmap <- function(df, title) {
  
  x_lvls <- order_by_size(unique(df$mod2), unique(df$ds2))
  y_lvls <- order_by_size(unique(df$mod1), unique(df$ds1))
  x_hex  <- get_hex(x_lvls, unique(df$ds2))
  y_hex  <- get_hex(y_lvls, unique(df$ds1))
  
  df <- df |>
    mutate(
      mod1 = factor(mod1, levels = y_lvls),
      mod2 = factor(mod2, levels = x_lvls)
    )
  
  # Colour swatch data frames — placed just outside the panel via clip="off"
  # x swatches sit one unit below y=0 (below the bottom row)
  x_swatch <- tibble(
    mod2 = factor(x_lvls, levels = x_lvls),
    hex  = x_hex
  )
  # y swatches sit one unit left of x=0 (left of the leftmost column)
  y_swatch <- tibble(
    mod1 = factor(y_lvls, levels = y_lvls),
    hex  = y_hex
  )
  
  ggplot(df, aes(x = mod2, y = mod1, fill = log10p_plot)) +
    
    # ── Main heatmap ──────────────────────────────────────
    geom_tile(colour = "grey70", linewidth = 0.15) +
    
    # ── x-axis colour swatches ────────────────────────────
    geom_tile(
      data        = x_swatch,
      aes(x = mod2, y = 0),
      fill        = x_swatch$hex,
      colour      = "white",
      height      = 0.8,
      width       = 0.8,
      linewidth   = 0.15,
      inherit.aes = FALSE
    ) +
    
    # ── y-axis colour swatches ────────────────────────────
    geom_tile(
      data        = y_swatch,
      aes(x = 0, y = mod1),
      fill        = y_swatch$hex,
      colour      = "white",
      height      = 0.8,
      width       = 0.8,
      linewidth   = 0.15,
      inherit.aes = FALSE
    ) +
    
    scale_fill_gradientn(
      colours  = c("white", FILL_MID, FILL_HIGH),
      values   = rescale(c(0, scale_max * 0.5, scale_max * 1.05)),
      limits   = c(0, scale_max * 1.05),
      na.value = FILL_NA,
      name     = "-log10(p)",
      guide    = guide_colorbar(
        barwidth       = 0.35,
        barheight      = 3.5,
        ticks          = FALSE,
        title.position = "top",
        title.hjust    = 0.5
      )
    ) +
    
    scale_x_discrete(limits = x_lvls) +
    scale_y_discrete(limits = y_lvls) +
    
    labs(title = title,
         x = unique(df$ds2),
         y = unique(df$ds1)) +
    
    # clip="off" lets swatches render outside the panel area
    coord_cartesian(clip = "off") +
    
    theme_minimal(base_size = 6) +
    theme(
      plot.title      = element_text(size = 6.5, face = "bold",
                                     hjust = 0.5, margin = margin(b = 2)),
      axis.title.x    = element_text(size = 5.5, face = "bold",
                                     colour = "grey30"),
      axis.title.y    = element_text(size = 5.5, face = "bold",
                                     colour = "grey30", angle = 90),
      # Plain black text — colour identity shown by swatches
      axis.text.x     = element_text(size = 7, angle = 90,
                                     hjust = 1, vjust = 0.5,
                                     colour = "grey10"),
      axis.text.y     = element_text(size = 7, colour = "grey10"),
      axis.ticks      = element_blank(),
      panel.grid      = element_blank(),
      legend.position = "right",
      legend.title    = element_text(size = 5),
      legend.text     = element_text(size = 4.5),
      plot.background = element_rect(fill = "white", colour = NA),
      # Extra bottom/left margin so swatches are not clipped by plot edge
      plot.margin     = margin(4, 4, 16, 16)
    )
}


# ============================================================
# STEP 7.  BUILD ALL 6 PLOTS
# ============================================================

grid_order <- c(
  "HOM_YNG:HOM_OLD", "HOM_YNG:SYN_YNG", "HOM_YNG:SYN_OLD",
  "HOM_OLD:SYN_YNG", "HOM_OLD:SYN_OLD", "SYN_YNG:SYN_OLD"
)
grid_order <- grid_order[grid_order %in% comparison_names]
if (length(grid_order) == 0) grid_order <- comparison_names

plots <- lapply(grid_order, function(cmp) {
  df    <- plong |> filter(comparison == cmp)
  title <- paste(unique(df$ds1), "vs", unique(df$ds2))
  make_heatmap(df, title)
})

# Legend only on last plot
plots <- lapply(seq_along(plots), function(i) {
  plots[[i]] + theme(legend.position = if (i == length(plots)) "right" else "none")
})

full_fig <- wrap_plots(plots, ncol = 3)


# ============================================================
# STEP 8.  SAVE
# ============================================================

ggsave("module_preservation_heatmap.pdf",
       plot = full_fig, device = cairo_pdf,
       width = 14, height = 10)

ggsave("module_preservation_heatmap.png",
       plot = full_fig, device = "png",
       width = 14, height = 10, dpi = 300, bg = "white")
