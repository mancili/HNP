GetMatchedNP = function (np) {
  cano$`Neuropeptide auto-annotation` = strsplit(cano$`Neuropeptide auto-annotation`, " ") 
  result = vector(mode ='list')
  count_np = c()
  length(result) = nrow(cano)
  for (i in 1:nrow(cano)) {
    dat = cano$`Neuropeptide auto-annotation`[i]
    if (is.na(dat) == TRUE){
      nnp = 0
    } else {
      nnp = length(c((dat)[[1]]))
      for (j in 1:nnp){
        matches <- dat[[1]][j] %in% np
        # Use `ifelse` to keep matches and remove non-matches
        result[[i]][j] <- ifelse(matches, dat[[1]][j], NA)
      }
      result[[i]] <- result[[i]][!is.na(result[[i]])]
      count_np[i] = length(result[[i]])
    }
    count_np[is.na(count_np)] <-0
  }
  return(count_np)
}

GetSectionCellCount = function(data, region) {
  df = data %>% 
    select(region, ncell_total) %>%
    separate(region, into = c('region', 'percentage'), sep = ': ') %>% 
    separate(percentage, into = c('percent'), sep = '%') 
  df$percent = as.numeric(df$percent)
  df$ncell = round(df$percent/100 * df$ncell_total)
  df[,1] = gsub("^\\s", "", df[,1])
  return(df)
}

library(tidyverse)
library(readxl)
library(ggplot2)

cano = read_xlsx('cluster_annotation.xlsx')
cano = cano[-nrow(cano),]
np = read.csv('NPs.csv')
np = np$Origin

adnp_neu = c('PENK','PDYN','CCK','SST','NPY','PYY','NMB','VIP',
             'TAC1','CHGA','CHGB','SCG2','SCG3','SCG5','VGF','NPB',
             'NPW','NXPH2','IGF1','PTHLH','APLN','CBLN4','CRH')
adnp_neu_count = GetMatchedNP(adnp_neu)
np_count = GetMatchedNP(np)
nonadnp = np_count-adnp_neu_count

sel = which(adnp_neu_count>=5)
antisel= which(nonadnp >= 3)

sasel = setdiff(sel, antisel)
view(cano[sasel,])
cano = cano[sasel,]

topregions = cano$`Top three regions` %>% data.frame()
colnames(topregions) = 'top_regions'
topregions = separate(topregions, top_regions, into = c('R1','R2','R3'), sep = ',')
topregions$ncell_total = cano$`Number of cells`-cano$H18.30.001

topdissections = cano$`Top three dissections` %>% data.frame()
colnames(topdissections) = 'top_dissections'
topdissections = separate(topdissections, top_dissections, into = c('D1','D2','D3'), sep = ',')
topdissections$ncell_total = cano$`Number of cells`-cano$H18.30.001


R1 = GetSectionCellCount(topregions, 'R1')
R2 = GetSectionCellCount(topregions, 'R2')
R3 = GetSectionCellCount(topregions, 'R3')
uregions = unique(c(R1$region,R2$region, R3$region)) %>% data.frame()
for (i in 1:nrow(uregions)) {
  sel1 = which(R1$region == uregions[i,1])
  sum1 = sum(R1$ncell[sel1])
  sel2 = which(R2$region == uregions[i,1])
  sum2 = sum(R2$ncell[sel2])
  sel3 = which(R3$region == uregions[i,1])
  sum3 = sum(R3$ncell[sel3])
  sum = sum(sum1, sum2, sum3)
  uregions[i,2] = sum
}
colnames(uregions) = c('region',"ncell")



D1 = GetSectionCellCount(topdissections, 'D1')
D2 = GetSectionCellCount(topdissections, 'D2')
D3 = GetSectionCellCount(topdissections, 'D3')
udissections = unique(c(D1$region,D2$region, D3$region)) %>% data.frame()
for (i in 1:nrow(udissections)) {
  sel1 = which(D1$region == udissections[i,1])
  sum1 = sum(D1$ncell[sel1])
  sel2 = which(D2$region == udissections[i,1])
  sum2 = sum(D2$ncell[sel2])
  sel3 = which(D3$region == udissections[i,1])
  sum3 = sum(D3$ncell[sel3])
  sum = sum(sum1, sum2, sum3)
  udissections[i,2] = sum
}

udissections = separate(udissections, '.', into = c('species','region'), sep = ' ')
udissections= udissections[,-1]

samples = read_excel('samples.xlsx')
regions = samples[,c("Dissection","Dissection abbreviation")]
regions = separate(regions, Dissection, into = c('big', 'rest'), sep = ' - ') %>% 
  data.frame()

sel = grep("Cerebral cortex", regions$big)
topcortex= unique(regions[sel,3])
sel = which(udissections$region %in% topcortex)
topcortex_dissections = udissections[sel, ]

write.csv(uregions, 'Table2.csv')
write.csv(topcortex_dissections, 'Table3.csv')

sel1 = grep("MEC", D1$region)
sel2 = grep("MEC", D2$region)
sel3 = grep("MEC", D3$region)
sel = c(sel1,sel2,sel3)

mec_hnp = cano[sel,]
mec_hnp = mec_hnp %>% select('Cluster name', 'Supercluster', 'Class auto-annotation',
                              'Neurotransmitter auto-annotation',"Neuropeptide auto-annotation",
                               "Subtype auto-annotation", "Transferred MTG Label", "Top Enriched Genes",
                               "Number of cells")
colnames(mec_hnp) = c('Cluster_name', 'Supercluster', 'Class', 'Neurotransmitter',
                      'Neuropeptide', 'Subtype', "MTG_label", "DEG", "ncell")
mec_hnp_NP = mec_hnp %>% group_by(Neuropeptide) %>% summarise(sum(ncell))
mec_hnp_NT = mec_hnp %>% group_by(Neurotransmitter) %>% summarise(sum(ncell))
mec_hnp_ST = mec_hnp %>% group_by(Subtype) %>% summarise(sum(ncell))
mec_hnp_MTG = mec_hnp %>% group_by(MTG_label) %>% summarise(sum(ncell))

colnames(mec_hnp_NT)[2] = 'ncell'
colnames(mec_hnp_ST)[2] = 'ncell'
colnames(mec_hnp_MTG)[2] = 'ncell'


png(filename = 'mec_hnp_NT.png', height = 5, width = 5, res = 600, units = 'in')
ggplot(mec_hnp_NT, aes(x="", y= ncell , fill= Neurotransmitter)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) +
  scale_fill_manual(values=c("#E29578","#FFDDD2","#006D77", "#83C5BE")) +
  theme_void() +
  theme(legend.position = "none")
dev.off()

png(filename = 'mec_hnp_ST.png', height = 5, width = 5, res = 600, units = 'in')
ggplot(mec_hnp_ST, aes(x="", y= ncell , fill= Subtype)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) +
  scale_fill_manual(values=c("#E29578","#FFDDD2","#006D77","#83C5BE","#9D9D9D")) +
  theme_void() +
  theme(legend.position = "none")
dev.off()

png(filename = 'mec_hnp_MTG.png', height = 5, width = 5, res = 600, units = 'in')
ggplot(mec_hnp_MTG, aes(x="", y= ncell , fill= MTG_label)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) +
  scale_fill_manual(values=c("#E29578",'#F2B342',"#FFDDD2", "#9D9D9D", "#006D77",
                             '#7AB1CC',"#83C5BE")) +
  theme_void() +
  theme(legend.position = "none")
dev.off()

