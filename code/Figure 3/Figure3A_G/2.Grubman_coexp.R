library(Seurat) # requires Seurat version 5
library(patchwork)
library(dplyr)
library(tidyverse)
library(readxl)
library(viridis)
library(ggridges)
library(cowplot)
library(pracma)
library(pheatmap)
library(BRETIGEA)
library(tidyr)
library(Matrix)
library(knitr)
library(scCustomize)
library(scuttle)

# Load Seurat object and add metadata ----
data.big <- Read10X_Multi_Directory(
  "../adnp_supplemental/data_grubman",
  default_10X_path = FALSE,
  sample_list = c("ct_2", "ct_4", "ct_6", 
                  "AD_2", "AD_4", "AD_6"),
  merge = T,
) %>%
  # read in filters
  # only features found in at least 3 cells
  # only nuclei with at least 200 features
  CreateSeuratObject(., min.cells = 3, min.features = 200) %>%
  # Disease Group
  AddMetaData(., str_extract(.@meta.data$orig.ident, 
                             "^(AD|ct)"), 
              col.name = "type") %>%
  # Hemoglobin gene content
  AddMetaData(., PercentageFeatureSet(., "^HB[^(P)]"), 
              col.name = "percent_hb") %>%
  # MT, ribo genes;
  Add_Cell_QC_Metrics(seurat_object = ., species = "human")

# Define project ----
Project(seurat_big) <- "AD_TCX"

# Format ID metadata ----
seurat_big@meta.data$type <- factor(seurat_big@meta.data$type, 
                                    levels = c("ct", "AD"))

data.big$orig.ident <- factor(data.big$orig.ident,
                              levels = c('ct1_2', 
                                         'ct3_4', 
                                         'ct5_6', 
                                         'AD1_2', 
                                         'AD3_4', 
                                         'AD5_6'))

## Check dimensions
dim(data.big)
table(data.big$orig.ident)

# Filtering, normalization, scaling ----
all_genes <- rownames(seurat_big)

data.big <- data.big %>%
  subset(., subset = nFeature_RNA < 2500 &
           nFeature_RNA > 200 &
           percent_mito < 5
  ) %>%
  NormalizeData(.) %>%
  FindVariableFeatures(., 
                       selection.method = "vst",
                       nfeatures = 2000
  ) %>%
  ScaleData(., features = all_genes) %>% 
  RunPCA(., features = VariableFeatures(object = .))

## Elbow plot ----
ElbowPlot(seurat_filter)

# Clustering ----
data.big <- data.big %>%
  FindNeighbors(., dims = 1:9, seed.use = 100) %>%
  FindClusters(., resolution = 0.2) %>%
  RunUMAP(., dims = 1:9, seed.use = 100) 

# Add cell types ----

ct_res_in = read.csv('cell_type.csv') #this file was obtained by running above and 1. AssignCellType in R 4.0 to export the cell type file

## Add cluster to metadata
data.big <- AddMetadata(data.big, ct_res_in$cluster, 'cluster')

## Add combined phenotype-cell type to meta
data.big <- AddMetadata(data.big,
                        paste(data.big@meta.data$type, 
                              data.big@meta.data$cell.type, 
                              sep = '_'), 
                        col.name = "phenotype_cell.type")

other <- read_xlsx('allnp.xlsx')

c <- intersect(toupper(other$NPs), 
               rownames((data.big@assays[["RNA"]]@features)))

allother <- FetchData(object = data.big,
                      vars = c,
                      layer = "data") |>
  mutate(
    phenotype = data.big$type,
    cell_type = factor(
      data.big$cell.type,
      levels = c('neu', 'ast', 
                 'mic', 'oli', 
                 'opc', 'end', 
                 'hybrid', 'no_id')
    ),
    phenotype_cell.type = data.big$phenotype_cell.type,
    origin_cell.type = paste(data.big$orig.ident, 
                             data.big$cell.type, 
                             sep = '_'),
    total = rowSums(allother[, 1:60]),
    no.coexpression = (apply(allother[, 1:60], 1, 
                             function(x)
      sum(x != 0)) %>% 
        data.frame())$. 
  )

# Assign NP coexpression groups ----
allother <- allother |> 
  mutate(coexp = if_else(cell_type == "neu",
                         case_when(no.coexpression >= 6 ~ 'high',
                                   no.coexpression <= 1 ~ 'low',
                                   TRUE ~ 'mid'),
                         case_when(no.coexpression >= 3 ~ 'high',
                                   no.coexpression <= 1 ~ 'low',
                                   TRUE ~ 'mid')))

data.big <- AddMetaData(data.big, allother$coexp, 'coexp')

for (i in c(1:nrow(allother))) {
  if (grepl('neu', allother[i, 'cell_type']) == TRUE) {
    if (allother[i, 'no.coexpression'] >= 6) {
      allother[i, 67] <- ('high')
    } else if (allother[i, 'no.coexpression'] <= 1) {
      allother[i, 67] <- ('low')
    } else {
      allother[i, 67] <- ('mid')
    }
  } else {
    if (allother[i, 'no.coexpression'] >= 3) {
      allother[i, 67] <- ('high')
    } else if (allother[i, 'no.coexpression'] <= 1) {
      allother[i, 67] <- ('low')
    } else {
      allother[i, 67] <- ('mid')
    }
  }
}

# Coexpression correlation testing ----
total.ct.neu <- allother %>%
  filter(cell_type == 'neu' & phenotype == 'ct')

total.AD.neu <- allother %>%
  filter(cell_type == 'neu' & phenotype == 'AD')

cor.test(total.ct.neu$total, 
         total.ct.neu$no.coexpression, 
         method = "spearman")

cor.test(total.AD.neu$total, 
         total.AD.neu$no.coexpression, 
         method = "spearman")

## Plotting ----

### Ct
png(
  filename = 'CT-cor.png',
  height = 5,
  width = 5,
  res = 600,
  units = 'in'
)
total.ct.neu %>%
  ggplot(aes(x = no.coexpression, y = total)) +
  geom_jitter(size = 1, col = '#83c5be') +
  geom_smooth(method = lm,
              se = FALSE,
              col = '#006d77') +
  xlim(0, 15) +
  ylim(0, 25) + #20000
  theme_bw()
dev.off()

### AD
png(
  filename = 'AD-cor.png',
  height = 5,
  width = 5,
  res = 600,
  units = 'in'
)
total.AD.neu %>%
  ggplot(aes(x = no.coexpression, y = total)) +
  geom_jitter(size = 1, col = '#ffddd2') +
  geom_smooth(method = lm,
              se = FALSE,
              col = '#e29578') +
  xlim(0, 15) +
  ylim(0, 25) + #20000
  theme_bw()
dev.off()

# Differential gene expression analysis ----
allother <- allother |>
  mutate(phenotype_cell.type_coexp = paste0(allother$phenotype_cell.type, '_', allother$coexp))

data.big <- AddMetaData(data.big, allother |> 
                          dplyr::select(phenotype_cell.type_coexp,
                                        phenotype_cell.type_no.coexp))

Idents(data.big) <- 'phenotype_cell.type_coexp'
levels(data.big)

## Interested in both increased and decrease expression
de.markers.neu_mid <- FindMarkers(
  data.big,
  ident.1 = "AD_neu_mid",
  ident.2 = "ct_neu_mid",
  min.cells.feature = 3,
  min.cells.group = 3,
  logfc.threshold = 0.25
)

de.genes.neu_mid.inc = filter(de.markers.neu_mid, 
                              p_val_adj < 0.05 &
                                avg_log2FC > 0)
write.csv(de.genes.neu_mid.inc, '02.de.genes.neu_mid.inc.csv')
#https://version-11-5.string-db.org/cgi/network?networkId=bQJkLvjoGu2q
de.genes.neu_mid.dec = filter(de.markers.neu_mid, 
                              p_val_adj < 0.05 &
                                avg_log2FC < 0)
write.csv(de.genes.neu_mid.dec, 
          '02.de.genes.neu_mid.dec.csv')
#https://version-11-5.string-db.org/cgi/network?networkId=blIJ9nBJKdj9

## Interested in decreased expression only 
de.markers.ctneu_low <- FindMarkers(
  data.big,
  ident.1 = "ct_neu_low",
  ident.2 = "ct_neu_high",
  min.cells.feature = 3,
  min.cells.group = 3,
  logfc.threshold = 0.25
)

de.markers.ctneu_low.inc <- filter(de.markers.ctneu_low, 
                                  p_val_adj < 0.05 &
                                    avg_log2FC > 0)
write.csv(de.markers.ctneu_low.inc, 
          '0.de.markers.ctneu_low.inc.csv')


de.markers.ctneu_low.dec <- filter(de.markers.ctneu_low, 
                                  p_val_adj < 0.05 &
                                    avg_log2FC < 0)
write.csv(de.markers.ctneu_low.dec, 
          '0.de.markers.ctneu_low.dec.csv')
#https://version-11-5.string-db.org/cgi/network?networkId=bOsb73IchYQd



ct_neu_dec <- read.csv('GTall.csv')
ct_neu_dec <- filter(ct_neu_dec, group != 'func')

ct_neu <- subset(data.big, 
                idents = c('ct_neu_low', 
                           'ct_neu_mid', 
                           'ct_neu_high'))

sel <- unique(ct_neu_dec$gene)

sel <- sel[which(sel %in% rownames(de.markers.ctneu_low.dec))]

ct_neu <- ct_neu[sel, ]

Idents(ct_neu) <- "phenotype_cell.type_coexp"

levels(ct_neu) <- c('ct_neu_low', 'ct_neu_mid', 'ct_neu_high')

png(
  filename = '0.ct_neu-low.mid.high.png',
  height = 10,
  width = 8,
  res = 600,
  units = 'in'
)
DoHeatmap(subset(ct_neu, downsample = 100),
          features = sel,
          size = 3) + #,
  # group.colors = c('#1d3557', '#fb8500','#ffd60a'))+ #'#e29578',,'#fca311'
  scale_fill_viridis(option = 'inferno')
dev.off()

ave.exp <- AverageExpression(object = ct_neu, layer = "data")

ct_neu_exp <- t(GetAssayData(object = ct_neu, 
                             layer = "data")) %>% 
  data.frame()

ct_neu_exp$no.coexp <- ct_neu@meta.data$phenotype_cell.type_no.coexp

sel <- grep('^[1-9]$|^10$', ct_neu_exp$no.coexp)

ct_neu_exp <- ct_neu_exp[sel, ]

lm_ct_neu_exp <- data.frame()

for (i in c(1:(ncol(ct_neu_exp) - 1))) {
  sel = which(ct_neu_exp[, i] >= 0)
  relation <- glm(ct_neu_exp[sel, i] ~ ct_neu_exp[sel, 'no.coexp'])
  lm_ct_neu_exp[i, 1] = relation$coefficients['ct_neu_exp[sel, "no.coexp"]']
  deviance <- summary(relation)$deviance
  null_deviance <- summary(relation)$null.deviance
  rsquared <- 1 - (deviance / null_deviance)
  lm_ct_neu_exp[i, 2] = rsquared
  lm_ct_neu_exp[i, 3] = coef(summary(relation))[2, 4]
}

print(summary(relation))

rownames(lm_ct_neu_exp) <- colnames(ct_neu_exp) [1:(ncol(ct_neu_exp) - 1)]
colnames(lm_ct_neu_exp) <- c('Coefficient', 'R_squared', 'P_value')

lm_ct_neu_exp$P_value <- p.adjust(lm_ct_neu_exp$P_value, 
                                  method = 'BH')
write.csv(lm_ct_neu_exp, 'glm.csv')

png(
  filename = 'ERBB4.png',
  width = 3,
  height = 5,
  res = 600,
  units = 'in'
)
ggplot(data = ct_neu_exp, aes(x = no.coexp, y = ERBB4)) +
  geom_jitter(colour = 'gray') +
  geom_smooth(method = 'lm', colour = 'black') +
  scale_x_continuous(breaks = scales::pretty_breaks()) +
  theme_bw()
dev.off()

sessionInfo()
