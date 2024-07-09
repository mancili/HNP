library(GOplot)
library(tidyverse)
library(viridis)
library(readxl)

# Data readin ----

## STRING output
en.in <- read_xlsx('plot_vs_mid.xlsx')
en.de <- read_xlsx('plot_vs_mid_dec.xlsx')
en.dat <- rbind(en.in, en.de)

## STRING input
dat.in <- read.csv('02.de.genes.neu_mid.inc.csv') 
dat.de <- read.csv('02.de.genes.neu_mid.dec.csv')
gene = rbind(dat.in, dat.de)

## Define color palette
mycol <- plasma(n = (nrow(en.dat)), 
                direction = -1)

dat <- en.dat

# Format genes and descriptions ----

gene_list <- map(dat$`matching proteins in your network (labels)`,
                 str_split,
                 pattern = ",") %>%
  flatten()

names(gene_list) = c(dat$`term description`)

## Unique gene names
sel = unique(unlist(gene_list))

gene_sel = which(gene$X %in% sel)

gene_sel = gene[gene_sel, ]

## Generate gene matrix ----

gene_matrix <- matrix(nrow = nrow(gene_sel), 
                     ncol = length(gene_list)) %>%
  data.frame()

colnames(gene_matrix) <- names(gene_list)
rownames(gene_matrix) <- gene_sel$X

for (i in c(1:ncol(gene_matrix))) {
  for (j in c(1:nrow(gene_matrix))) {
    if (rownames(gene_matrix)[j] %in% gene_list[[i]]) {
      gene_matrix[j,i] = 1
    } else {
      gene_matrix[j,i] = 0
    }
  }
}

gene_matrix$logFC = gene_sel$avg_log2FC

# Chord diagram ----

png(
  filename = 'Goplot_mid.png',
  width = 13,
  height = 15,
  units = 'in',
  res = 600
)

GOChord(
  gene_matrix,
  space = 0.02,
  gene.order = 'logFC',
  ribbon.col = mycol,
  lfc.min = -10,
  lfc.max = 10,
  #lfc.col	= mycol2,
  gene.space = 0.25,
  gene.size = 5
)
dev.off()
