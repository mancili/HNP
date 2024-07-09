library(GOplot)
library(tidyverse)
library(viridis)

# Read in data ----

mempo = read.delim('MembranePotentialCluster.tsv')
mito = read.delim('MitochondriaCluster.tsv')
atp = read.delim('ATPaseCluster.tsv')
trans = read.delim('TranslationCluster.tsv')
mhc = read.delim('MHCCluster.tsv')
func = read.delim('MemoryCluster.tsv')

gene = read.csv('0.de.markers.ctneu_low.dec.csv')

# Get gene terms ----

## Custom function ----
GetGeneTerm <- function(x) {
  gene_df <- x |>
    group_by(term.description) |>
    nest() |>
    mutate(genes = map(
      data,
      \(x) str_split(x$matching.proteins.in.your.network..labels., 
                     pattern = ',') |>
        unlist() |> 
        data.frame() |> 
        set_names("gene") |> 
        rowid_to_column() |>  
        mutate(term = term.description) |> 
        dplyr::select(-rowid)
    )) |> 
    pluck('genes') |> 
    list_rbind()
  return(gene_df)
}

## Run ----

term_objects <- list(mempo, mito, atp, trans, mhc, func) |>
  set_names(c('mempo', 'mito', 'atp', 'trans', 'mhc', 'func'))

GTall <- imap(term_objects,
    \(x,y) GetGeneTerm(x) |> 
      mutate(group = names(term_objects[y]))) |> 
  list_rbind() |> 
  unique()

write.csv(GTall, 'GTall.csv')