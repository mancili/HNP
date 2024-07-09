#Functions
#-------------------------------------------------------------------------------------

#function to get correlation between NP transcript and NP
GetResult_Cor_MIT  = function(input_data, expression_threshold) {
  data = input_data
  col = ncol(data) - 4
  gene_cols <- colnames(data)[2:col]
  
  # Get the unique donor IDs for each Braak stage
  ct_donors <- unique(data$donor_id[data$status == 'ct'])
  mci_donors <- unique(data$donor_id[data$status == 'mci'])
  ad_donors <- unique(data$donor_id[data$status == 'ad'])
  
  # Initialize a result_cor data frame to store the test result_cors
  result_cor <- data.frame(
    Gene = character(),
    ct_Prop = numeric(),
    mci_Prop = numeric(),
    ad_Prop = numeric(),
    P_Value = numeric(),
    stringsAsFactors = FALSE
  )
  
  # Iterate over each gene
  for (gene in gene_cols) {
    ct_props <- numeric(length(ct_donors))
    mci_props <- numeric(length(mci_donors))
    ad_props <- numeric(length(ad_donors))
    
    for (i in seq_along(ct_donors)) {
      donor <- ct_donors[i]
      donor_adata.exp <- data[data$donor_id == donor &
                                data$status == 'ct', ]
      ct_props[i] <- sum(donor_adata.exp[[gene]] > expression_threshold) / nrow(donor_adata.exp)
    }
    
    for (i in seq_along(mci_donors)) {
      donor <- mci_donors[i]
      donor_adata.exp <- data[data$donor_id == donor &
                                data$status == 'mci', ]
      mci_props[i] <- sum(donor_adata.exp[[gene]] > expression_threshold) / nrow(donor_adata.exp)
    }
    
    for (i in seq_along(ad_donors)) {
      donor <- ad_donors[i]
      donor_adata.exp <- data[data$donor_id == donor &
                                data$status == 'ad', ]
      ad_props[i] <- sum(donor_adata.exp[[gene]] > expression_threshold) / nrow(donor_adata.exp)
    }
    
    # Create a data frame with the proportions and the corresponding Braak stage
    prop_data <- data.frame(
      Proportion = c(ct_props, mci_props, ad_props),
      Diagnosis = rep(c("ct", "mci", "ad"), times = c(
        length(ct_props), length(mci_props), length(ad_props)
      ))
    )
    
    # Perform Spearman correlation test with "less" alternative hypothesis
    cor_test <- cor.test(
      prop_data$Proportion,
      as.numeric(factor(
        prop_data$Diagnosis, levels = c("ct", "mci", "ad")
      )),
      method = "spearman",
      alternative = "less"
    )
    
    result_cor <- rbind(
      result_cor,
      data.frame(
        Gene = gene,
        ct_Prop = mean(ct_props),
        mci_Prop = mean(mci_props),
        ad_Prop = mean(ad_props),
        P_Value = cor_test$p.value,
        stringsAsFactors = FALSE
      )
    )
  }
  
  result_cor$Adjusted_P_Value <- p.adjust(result_cor$P_Value, method = "BH")
  return (result_cor)
}


#function to get sum of transcript value and coexpression for each cell/neuron
GetResult_Cor_Status = function (input_data, expression_threshold) {
  data = input_data
  col = ncol(data) - 4
  gene_cols <- colnames(data)[2:col]
  result_cor_status <- data.frame(
    Aggregated_Value = numeric(),
    Coexpressed_Count = integer(),
    stringsAsFactors = FALSE
  )
  
  for (i in 1:nrow(data)) {
    cell_values <- data[i, gene_cols]
    
    coexpressed_count <- sum(cell_values > expression_threshold)
    
    aggregated_value <- sum(cell_values[cell_values > expression_threshold])
    
    result_cor_status <- rbind(
      result_cor_status,
      data.frame(
        Aggregated_Value = aggregated_value,
        Coexpressed_Count = coexpressed_count,
        stringsAsFactors = FALSE
      )
    )
  }
  
  result_cor_status$Cell = data$X
  result_cor_status = cbind(result_cor_status, data[, (col + 1):(col + 4)])
  result_cor_status$status = factor(result_cor_status$status, levels = c("ct", "mci", "ad"))
  return (result_cor_status)
}




#function to get correlation between NP transcript and NP
GetResult_Cor_Leng  = function (input_data, expression_threshold) {
  data = input_data
  col = ncol(data) - 2
  gene_cols <- colnames(data)[2:col]
  
  # Get the unique donor IDs for each Braak stage
  ct_donors <- unique(data$donor_id[data$BraakStage == '0'])
  mci_donors <- unique(data$donor_id[data$BraakStage == '2'])
  ad_donors <- unique(data$donor_id[data$BraakStage == '6'])
  
  # Initialize a result_cor data frame to store the test result_cors
  result_cor <- data.frame(
    Gene = character(),
    ct_Prop = numeric(),
    mci_Prop = numeric(),
    ad_Prop = numeric(),
    P_Value = numeric(),
    stringsAsFactors = FALSE
  )
  
  # Iterate over each gene
  for (gene in gene_cols) {
    ct_props <- numeric(length(ct_donors))
    mci_props <- numeric(length(mci_donors))
    ad_props <- numeric(length(ad_donors))
    
    for (i in seq_along(ct_donors)) {
      donor <- ct_donors[i]
      donor_adata.exp <- data[data$donor_id == donor &
                                data$BraakStage == '0', ]
      ct_props[i] <- sum(donor_adata.exp[[gene]] > expression_threshold) / nrow(donor_adata.exp)
    }
    
    for (i in seq_along(mci_donors)) {
      donor <- mci_donors[i]
      donor_adata.exp <- data[data$donor_id == donor &
                                data$BraakStage == '2', ]
      mci_props[i] <- sum(donor_adata.exp[[gene]] > expression_threshold) / nrow(donor_adata.exp)
    }
    
    for (i in seq_along(ad_donors)) {
      donor <- ad_donors[i]
      donor_adata.exp <- data[data$donor_id == donor &
                                data$BraakStage == '6', ]
      ad_props[i] <- sum(donor_adata.exp[[gene]] > expression_threshold) / nrow(donor_adata.exp)
    }
    
    # Create a data frame with the proportions and the corresponding Braak stage
    prop_data <- data.frame(
      Proportion = c(ct_props, mci_props, ad_props),
      Diagnosis = rep(c('0', '2', '6'), times = c(
        length(ct_props), length(mci_props), length(ad_props)
      ))
    )
    
    # Perform Spearman correlation test with "less" alternative hypothesis
    cor_test <- cor.test(
      prop_data$Proportion,
      as.numeric(factor(
        prop_data$Diagnosis, levels = c('0', '2', '6')
      )),
      method = "spearman",
      alternative = "less"
    )
    
    result_cor <- rbind(
      result_cor,
      data.frame(
        Gene = gene,
        ct_Prop = mean(ct_props),
        mci_Prop = mean(mci_props),
        ad_Prop = mean(ad_props),
        P_Value = cor_test$p.value,
        stringsAsFactors = FALSE
      )
    )
  }
  
  result_cor$Adjusted_P_Value <- p.adjust(result_cor$P_Value, method = "BH")
  return (result_cor)
}



GetResult_Test <- function(input_data, expression_threshold) {
  data <- input_data
  col = ncol(data) - 2
  gene_cols <- colnames(data)[2:col]
  
  # Get the unique donor IDs for each Braak stage
  ct_mci_donors <- unique(data$donor_id[data$BraakStage %in% c('0')])
  ad_donors <- unique(data$donor_id[data$BraakStage == '6'])
  
  # Initialize a result_cor data frame to store the test result_cors
  result_cor <- data.frame(
    Gene = character(),
    ct_mci_Prop = numeric(),
    ad_Prop = numeric(),
    P_Value = numeric(),
    stringsAsFactors = FALSE
  )
  
  # Iterate over each gene
  for (gene in gene_cols) {
    ct_mci_props <- numeric(length(ct_mci_donors))
    ad_props <- numeric(length(ad_donors))
    
    for (i in seq_along(ct_mci_donors)) {
      donor <- ct_mci_donors[i]
      donor_adata.exp <- data[data$donor_id == donor &
                                data$BraakStage %in% c('0'), ]
      ct_mci_props[i] <- sum(donor_adata.exp[[gene]] > expression_threshold) / nrow(donor_adata.exp)
    }
    
    for (i in seq_along(ad_donors)) {
      donor <- ad_donors[i]
      donor_adata.exp <- data[data$donor_id == donor &
                                data$BraakStage == '6', ]
      ad_props[i] <- sum(donor_adata.exp[[gene]] > expression_threshold) / nrow(donor_adata.exp)
    }
    
    # Perform Wilcoxon test with "greater" alternative hypothesis
    wilcox_test <- wilcox.test(ct_mci_props, ad_props, alternative = "greater")
    
    result_cor <- rbind(
      result_cor,
      data.frame(
        Gene = gene,
        ct_mci_Prop = mean(ct_mci_props),
        ad_Prop = mean(ad_props),
        P_Value = wilcox_test$p.value,
        stringsAsFactors = FALSE
      )
    )
  }
  
  result_cor$Adjusted_P_Value <- p.adjust(result_cor$P_Value, method = "BH")
  return(result_cor)
}


#function to get sum of transcript value and coexpression for each cell/neuron
GetResult_Cor_BraakStage = function (input_data, expression_threshold) {
  data = input_data
  col = ncol(data) - 2
  gene_cols <- colnames(data)[2:col]
  result_cor_BraakStage <- data.frame(
    Aggregated_Value = numeric(),
    Coexpressed_Count = integer(),
    stringsAsFactors = FALSE
  )
  
  for (i in 1:nrow(data)) {
    cell_values <- data[i, gene_cols]
    
    coexpressed_count <- sum(cell_values > expression_threshold)
    
    aggregated_value <- sum(cell_values[cell_values > expression_threshold])
    
    result_cor_BraakStage <- rbind(
      result_cor_BraakStage,
      data.frame(
        Aggregated_Value = aggregated_value,
        Coexpressed_Count = coexpressed_count,
        stringsAsFactors = FALSE
      )
    )
  }
  
  result_cor_BraakStage$Cell = data$X
  result_cor_BraakStage = cbind(result_cor_BraakStage, data[, (col + 1):(col +
                                                                           2)])
  result_cor_BraakStage$BraakStage = factor(result_cor_BraakStage$BraakStage, levels = c('0', '2', '6'))
  return (result_cor_BraakStage)
}

#------------
#Main Text
#----------------------------------------------------------------------------------------------------

#Figure 3B: MIT Multiomics

library(tidyverse)
library(DescTools)

adata.exp = read.csv('expMatrix.csv')
aadata.exp.meta = read.csv('adata.meta.csv')
cellmeta = read.csv('cell.meta.csv')

adata.exp$status = cellmeta$cogdx
adata.exp$braak = cellmeta$braaksc
adata.exp$majorcelltype = cellmeta$major.celltype

adata.exp$status <- ifelse(adata.exp$status == 1,
                           "ct",
                           ifelse(
                             adata.exp$status == 2,
                             "mci",
                             ifelse(adata.exp$status == 4, "ad", adata.exp$status)
                           ))

mcidonor = cellmeta[which(adata.exp$status == 'mci'), ]
mcidonor %>% select(projid, braaksc) %>% unique()


adata.exp$donor_id = aadata.exp.meta$projid

adata.exp.exc = filter(adata.exp, majorcelltype == 'Exc')
adata.exp.inh = filter(adata.exp, majorcelltype == 'Inh')

np = read.csv('ADNPs.csv')
np = np$.
adata.exp.np = adata.exp[, which(colnames(adata.exp) %in% c('X', np))]
adata.exp.np = cbind(adata.exp.np, adata.exp[, 94:97])

all.cor = GetResult_Cor_MIT(adata.exp.np, 2)
all.coexp = GetResult_Cor_Status(adata.exp.np, 2)


summary_data <- all.coexp %>%
  mutate(
    Coexpressed_Level = case_when(
      Coexpressed_Count %in% 0:1 ~ "low",
      Coexpressed_Count %in% 2:5 ~ "mid",
      Coexpressed_Count >= 6 ~ "high"
    )
  ) %>%
  group_by(donor_id, status, Coexpressed_Level) %>%
  summarise(num_rows = n()) %>%
  ungroup() %>%
  group_by(donor_id) %>%
  mutate(proportion = num_rows / sum(num_rows))

summary_data$Coexpressed_Level = factor(summary_data$Coexpressed_Level, levels = c('low', 'mid', 'high'))


summary_data_with_error <- summary_data %>%
  group_by(Coexpressed_Level, status) %>%
  summarise(
    mean_var = mean(proportion),
    n = sum(num_rows),
    se = sd(proportion) / sqrt(length(proportion))
  )

summary_data_trimmed <- summary_data %>%
  group_by(Coexpressed_Level, status) %>%
  arrange(proportion) %>%
  slice(2:(n() - 1)) %>%
  ungroup()


png(
  filename =  "maintext_cognition.png",
  width = 5,
  height = 8,
  units = 'in',
  res = 500
)
ggplot(summary_data_with_error,
       aes(x = Coexpressed_Level, y = mean_var, fill = status)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(
    aes(ymin = mean_var - se, ymax = mean_var + se),
    position = position_dodge(width = 0.9),
    width = 0.2
  ) +
  geom_point(
    data = summary_data,
    aes(x = Coexpressed_Level, y = proportion, group = status),
    position = position_jitterdodge(dodge.width = 0.9),
    size = 2,
    alpha = 0.4
  ) +
  scale_color_manual(values = c("#34623f", "#34623f", "#34623f")) +
  scale_fill_manual(values = c("#dbfeb8", "#98c9a3", "#14746f")) +
  facet_wrap( ~ Coexpressed_Level, ncol = 1) +
  ylim(0, 1) +
  labs(x = "Coexpressed Level", y = "Proportion of Rows", fill = "Status") +
  theme_bw() +
  theme(strip.background = element_blank(), strip.text = element_blank())
dev.off()


a = filter(summary_data, status == 'ct'  &
             Coexpressed_Level == 'low')
b = filter(summary_data, status == 'mci' &
             Coexpressed_Level == 'low')
c = filter(summary_data, status == 'ad'  &
             Coexpressed_Level == 'low')

wilcox.test(a$proportion, b$proportion, 'less')
wilcox.test(b$proportion, c$proportion, 'less')
wilcox.test(a$proportion, c$proportion, 'less')


a = filter(summary_data, status == 'ct'  &
             Coexpressed_Level == 'mid')
b = filter(summary_data, status == 'mci' &
             Coexpressed_Level == 'mid')
c = filter(summary_data, status == 'ad'  &
             Coexpressed_Level == 'mid')

wilcox.test(a$proportion, b$proportion, 'greater')
wilcox.test(b$proportion, c$proportion, 'greater')
wilcox.test(a$proportion, c$proportion, 'greater')

a = filter(summary_data, status == 'ct'  &
             Coexpressed_Level == 'high')
b = filter(summary_data, status == 'mci' &
             Coexpressed_Level == 'high')
c = filter(summary_data, status == 'ad'  &
             Coexpressed_Level == 'high')

wilcox.test(a$proportion, b$proportion, 'greater')
wilcox.test(b$proportion, c$proportion, 'greater')
wilcox.test(a$proportion, c$proportion, 'greater')




#Figure 3C: Leng dataset
library(readxl)
library(DescTools)
library(tidyverse)

exc_np = read.csv('ExpMatrix_exc_np.csv')
inh_np = read.csv('ExpMatrix_inh_np.csv')

data = rbind(exc_np, inh_np)

np = read.csv('ADNPs.csv')
np = np$.
data.np = data[, which(colnames(data) %in% c('X', np))]
data.np = cbind(data.np, data[, 94:95])


all.cor = GetResult_Cor_Leng(data.np, 2)
all.test = GetResult_Test(data.np, 2)
all.coexp = GetResult_Cor_BraakStage(data.np, 2)

result_cor_BraakStage_ct = filter(all.coexp, BraakStage == '0')
cor.test(
  result_cor_BraakStage_ct$Coexpressed_Count,
  result_cor_BraakStage_ct$Aggregated_Value,
  method = "spearman"
)

result_cor_BraakStage_mci = filter(all.coexp, BraakStage == '2')
cor.test(
  result_cor_BraakStage_mci$Coexpressed_Count,
  result_cor_BraakStage_mci$Aggregated_Value,
  method = "spearman"
)

result_cor_BraakStage_AD = filter(all.coexp, BraakStage == '6')
cor.test(
  result_cor_BraakStage_AD$Coexpressed_Count,
  result_cor_BraakStage_AD$Aggregated_Value,
  method = "spearman"
)


png(
  filename =  "linear_adnp_braak.png",
  width = 5,
  height = 8,
  units = 'in',
  res = 500
)
ggplot(all.coexp,
       aes(x = Coexpressed_Count, y = Aggregated_Value, color = BraakStage)) +
  geom_jitter(alpha = 0.7,
              width = 1,
              height = 1) +
  geom_smooth(method = "lm", se = FALSE, size = 1) +
  scale_color_manual(values = c(
    "0" = "#ffd670",
    "2" = "#f4a261",
    '6' = '#e76f51'
  )) +
  facet_wrap( ~ BraakStage, ncol = 1) +
  ylim(0, 200) +
  xlim(0, 10) +
  scale_x_continuous(breaks = scales::breaks_pretty()) +
  theme_bw() +
  theme(strip.background = element_blank(), strip.text = element_blank())
dev.off()



#------------
#Supplementary Figures
#----------------------------------------------------------------------------------------------------
#MIT

all.cor = GetResult_Cor_MIT(adata.exp, 2)
all.coexp = GetResult_Cor_Status(adata.exp, 2)

#plot these in the same graph
result_cor_status_ct = filter(all.coexp, status == 'ct')
cor.test(
  result_cor_status_ct$Coexpressed_Count,
  result_cor_status_ct$Aggregated_Value,
  method = "spearman"
)

result_cor_status_mci = filter(all.coexp, status == 'mci')
cor.test(
  result_cor_status_mci$Coexpressed_Count,
  result_cor_status_mci$Aggregated_Value,
  method = "spearman"
)

result_cor_status_AD = filter(all.coexp, status == 'ad')
cor.test(
  result_cor_status_AD$Coexpressed_Count,
  result_cor_status_AD$Aggregated_Value,
  method = "spearman"
)

png(
  filename =  "linear_allnp_cognition.png",
  width = 6,
  height = 12,
  units = 'in',
  res = 500
)
ggplot(all.coexp,
       aes(x = Coexpressed_Count, y = Aggregated_Value, color = status)) +
  geom_jitter(alpha = 0.7,
              width = 1,
              height = 1) +
  geom_smooth(method = "lm", se = FALSE, size = 1) +
  scale_color_manual(values = c(
    "ct" = "blue",
    "mci" = "pink",
    'ad' = 'maroon'
  )) +
  facet_wrap( ~ status, ncol = 1) +
  ylim(0, 250) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
  theme_bw() +
  theme(strip.background = element_blank(), strip.text = element_blank())
dev.off()


png(
  filename = "Distribution_coexpression_cognition.png",
  width = 12,
  height = 6,
  units = 'in',
  res = 500
)
ggplot(all.coexp,
       aes(
         x = Coexpressed_Count,
         fill = as.character(status),
         color = as.character(status)
       )) +
  geom_histogram(
    aes(y = after_stat(density)),
    position = "identity",
    binwidth = 1,
    boundary = 0,
    alpha = 0.75
  ) +
  xlim(c(0, 15)) +
  facet_grid(. ~ status) +
  scale_color_manual(values = c("#BD9FC2", "#906498", "#72417B")) +
  scale_fill_manual(values = c("#FFC7E6", "#D7C5D9", "#885990")) +
  scale_y_continuous(labels = scales::percent) +
  theme_bw() +
  theme(
    legend.position = "none",
    strip.text = element_text(size = 20),
    axis.text = element_text(size = 12)
  )
dev.off()


#Leng

all.cor = GetResult_Cor_Leng(data, 2)
all.test = GetResult_Test(data, 2)
all.coexp = GetResult_Cor_BraakStage(data, 2)


result_cor_BraakStage_ct = filter(all.coexp, BraakStage == '0')
cor.test(
  result_cor_BraakStage_ct$Coexpressed_Count,
  result_cor_BraakStage_ct$Aggregated_Value,
  method = "spearman"
)

result_cor_BraakStage_mci = filter(all.coexp, BraakStage == '2')
cor.test(
  result_cor_BraakStage_mci$Coexpressed_Count,
  result_cor_BraakStage_mci$Aggregated_Value,
  method = "spearman"
)

result_cor_BraakStage_AD = filter(all.coexp, BraakStage == '6')
cor.test(
  result_cor_BraakStage_AD$Coexpressed_Count,
  result_cor_BraakStage_AD$Aggregated_Value,
  method = "spearman"
)

png(
  filename =  "linear_allnp_braak.png",
  width = 6,
  height = 12,
  units = 'in',
  res = 500
)
ggplot(all.coexp,
       aes(x = Coexpressed_Count, y = Aggregated_Value, color = BraakStage)) +
  geom_jitter(alpha = 0.7,
              width = 1,
              height = 1) +
  geom_smooth(method = "lm", se = FALSE, size = 1) +
  scale_color_manual(values = c(
    "0" = "blue",
    "2" = "pink",
    '6' = 'maroon'
  )) +
  facet_wrap( ~ BraakStage, ncol = 1) +
  ylim(0, 200) +
  xlim(0, 10) +
  scale_x_continuous(breaks = scales::breaks_pretty()) +
  theme_bw() +
  theme(strip.background = element_blank(), strip.text = element_blank())
dev.off()


png(
  filename = "Distribution_coexpression.png",
  width = 12,
  height = 6,
  units = 'in',
  res = 500
)
ggplot(all.coexp,
       aes(
         x = Coexpressed_Count,
         fill = as.character(BraakStage),
         color = as.character(BraakStage)
       )) +
  geom_histogram(
    aes(y = after_stat(density)),
    position = "identity",
    binwidth = 1,
    boundary = 0,
    alpha = 0.75
  ) +
  xlim(c(0, 15)) +
  facet_grid(. ~ BraakStage) +
  scale_color_manual(values = c("#BD9FC2", "#906498", "#72417B")) +
  scale_fill_manual(values = c("#FFC7E6", "#D7C5D9", "#885990")) +
  scale_y_continuous(labels = scales::percent) +
  theme_bw() +
  theme(
    legend.position = "none",
    strip.text = element_text(size = 20),
    axis.text = element_text(size = 12)
  )
dev.off()