library(BRETIGEA)
library(Seurat) # version 4
library(dplyr)
library(stringr)
library(p)

#cluster ID #seurat4
# Extract cluster assignments from Seurat object ----

clusters <- data.big@meta.data$seurat_clusters %>%
  as.character() %>%
  as.numeric() %>%
  data.frame()
  
.rowNamesDF(clusters, make.names = FALSE) <- data.big@assays$RNA$data@Dimnames[[2]]

colnames(clusters) <- 'cluster_num'

# Convert numeric cluster assignments to cell type name ----

type_cluster_ids <- data.frame(cluster_num = 0:9,
                               cluster = c(rep('oli', 2),
                                           'ast',
                                           'opc',
                                           'ast',
                                           'neu',
                                           'mic',
                                           'neu',
                                           'hybrid',
                                           'end'))

clusters <- clusters %>%
  left_join(type_cluster_ids, by = "cluster_num") %>%
  dplyr::select(-cluster_num)

# png(file = 'newdim.celltype.png', width = 5, height = 5, units= 'in', res = 500)
# DimPlot(data.big, reduction = "umap", group.by = 'cell.type',
#         cols = c( 'neu' = '#F5CAC3',
#                   'ast' = '#F28482',
#                   'mic' = '#84A59D',
#                   'oli' = '#6C584C',
#                   'opc' = '#F6BD60',
#                   'end' = '#DEDEDE',
#                   'hybrid' = '#ECECEC'))
# dev.off()

# Estimate brain cell type proportions ----

dfpv <- GetAssayData(data.big, slot = "scale.data")

ct_res = brainCells(dfpv, nMarker = 50)

ct_res <- as.data.frame(ct_res)

rm(dfpv)

# Calculate module scores ----

## Custom function

add_marker_score <- function(x) {
  marker <- dplyr::filter(markers_df_brain, cell == x)
  
  marker_use <- list(c(marker[1:50, 1]))
  
  AddModuleScore(data.big, 
                 marker_use, 
                 name = paste(x, 'score', sep = '.'))
}

## Run ----

### Distinct marker names (ast, neu, etc.) from BRETIGEA
marker_names <- distinct(markers_df_brain, cell) %>%
  pluck(1)

### Plotting cutoff
plot_cutoff <- 'q10'

### add marker scores
for (x in marker_names) {
  data.big <- add_marker_score(x)
}

### write plots to objects
map(paste0(marker_names, ".score1"),
    \(x) FeaturePlot(data.big,
                     features = x,
                     min.cutoff = plot_cutoff)) %>%
  set_names(paste0("plot.", marker_names)) %>%
  list2env(.GlobalEnv)

# Find hybrids ----

f_id <- data.frame(matrix(ncol = 1, 
                          nrow = ncol(data.big)))

rownames(f_id) <- rownames(ct_res)

colnames(f_id) <- "ID"

for (i in c(1:ncol(data.big))) {
  max1 <- max(ct_res[i, ])
  
  max2 <- max(ct_res[i, ][ct_res[i, ] != max(ct_res[i, ])])
  
  mu <- (max1 - max2) / max1
  
  if (mu >= 0.2) {
    
    f_id[i, 1] <- rownames(data.frame(which.max(ct_res[i, ])))
    
  } else {
    
    f_id[i, 1] <- 'hybrid'
  }
}

ct_res_in <- cbind(ct_res, f_id)

ct_res_in <- cbind(ct_res_in, clusters)

# Cell number per phenotype ----

sel <- grep('^ct', rownames(ct_res))

## hybrid cells disproportionately from AD

length(grep('hybrid', f_id[sel, 1]))

length(grep('hybrid', f_id[-sel, 1]))

length(grep('^ct', rownames(ct_res_in))) #6946

length(grep('^AD', rownames(ct_res_in))) #7333

# Cell type-phenotype data exploration ----

## Hybrid ----

hybrid_type <- subset(ct_res_in, ID == 'hybrid' |
                        cluster == 'hybrid') %>% 
  select('ID', 'cluster')

## Custom function for standard cell types ----

subset_marker <- function(x) {
  type <- subset(ct_res_in, ID == x, 
                     select = c(x, 'ID'))
  
  nct <- length(grep('^ct', rownames(type)))
  nad <- length(grep('^AD', rownames(type)))
  
  type <- dplyr::arrange(type, desc(!!as.symbol(x)))
  type <- type[1:round((nct + nad) * 0.95), ]
  
  return(type)
}

## Run subset ----
map(marker_names,
   subset_marker) %>%
  set_names(paste0(marker_names, "_type")) %>%
  list2env(.GlobalEnv)

## Bind rows ----

cell_type <- bind_rows(hybrid_type,
                       ast_type,
                       end_type,
                       mic_type,
                       neu_type,
                       oli_type,
                       opc_type) %>% 
  select(ID)

# Unidentified ----

unID <- setdiff(rownames(ct_res_in), rownames(cell_type))

noID <- c(rep('no_id', times = length(unID))) %>% data.frame()

rownames(noID) <- unID

colnames(noID) <- 'ID'

all(rownames(noID) == unID)

sel <- which(!is.na(match(
  rownames(ct_res_in), rownames(hybrid_type)
)))

ct_res_in$cluster[sel] <- 'hybrid'

sel <- which(!is.na(match(rownames(ct_res_in), rownames(noID))))

ct_res_in$cluster[sel] <- 'no_id'

#create a new metadata for cell ID then color 
write.csv(ct_res_in, 'cell_type.csv')

data.big[['cell.type']] <- ct_res_in$cluster
