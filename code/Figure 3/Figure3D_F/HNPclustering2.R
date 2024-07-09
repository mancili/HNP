library(tidyverse)
library(DOSE)
library(enrichplot)
library(readxl)
library(tidyr)

# Data readin ----
cthigh = read_xlsx('plot_ct_high.xlsx')

# Extract high-level term data ----
hi.cat = unique(cthigh$HighLevel)

cthigh.plot <- cthigh |>
  group_by(HighLevel) |>
  nest() |>
  mutate(genes = map(
    data,
    \(x) str_split(x$`matching proteins in your network (labels)`,
                   pattern = ',') |>
      unlist()
  )) |> 
  pluck('genes') |> 
  set_names(hi.cat)

## Plot ----
png(
  filename = "goplot_highlevel.png",
  height = 5,
  width = 5,
  units = "in",
  res = 600
)
cnetplot(
  cthigh.plot,
  layout = 'kk',
  node_label = "none",
  color_category = "#006d77"
)
dev.off()

# Extract genes for categories of interest ----

## Custom function ----
select_cthigh <- function(x) {
  cthigh.sel <- cthigh |>
    filter(str_detect(HighLevel, x)) |>
    dplyr::select(`term description`, 
                  `matching proteins in your network (labels)`) |>
    mutate(genes = map(`matching proteins in your network (labels)`, str_split_1, pattern = ',') |> 
             set_names(`term description`)) |> 
    pluck("genes")
  
  return(cthigh.sel)
}

## Membrane ----

cthigh.sel.mem <- select_cthigh("Membrane")

### Plot ----
png(
  filename = "goplot_membrane.png",
  height = 4,
  width = 5,
  units = "in",
  res = 600
)
cnetplot(
  cthigh.sel.mem,
  showCategory = 8,
  layout = 'dh',
  node_label = "none",
  color_category = "#83c5be"
)
dev.off()

## Translation ----
cthigh.sel.trans <- select_cthigh("Translation")

### Plot with category ----
png(
  filename = "goplot_translation_categoryname.png",
  height = 4,
  width = 5,
  units = "in",
  res = 600
)
cnetplot(
  cthigh.sel.trans,
  showCategory = 8,
  layout = 'dh',
  node_label = "category",
  color_category = "#83c5be"
)
dev.off()

### Plot without category ----
png(
  filename = "goplot_translation_none.png",
  height = 4,
  width = 5,
  units = "in",
  res = 600
)
cnetplot(
  cthigh.sel.trans,
  showCategory = 8,
  layout = 'dh',
  node_label = 'none',
  color_category = "#83c5be"
)
dev.off()
