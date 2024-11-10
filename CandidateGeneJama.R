library(tidyverse)
library(tidyr)
library(table1)
library(sesame)
library(sesameData)
library(data.table)
library(GenomicRanges)
library(performance)
library(dplyr)
library(factoextra)
library(flextable)
library(broom)
library(officer)
library(lme4)
library(marginaleffects)


m <- read.csv("annotated_candidate_sites.csv")## annotated mvalue data file

d <- read.csv("data_raw/d.csv")## this is the pheno datafile with the final 91 subjects in it


##Make a long form of the m-value dataset
m_long <- m %>%
  pivot_longer(cols = starts_with("X20"), names_to = "Subject", values_to = "M_Value") %>%
 mutate(subject2 = substr(Subject, 2, nchar(Subject))) %>%
  left_join(d, by = c("subject2" = "Filename"))

unique_subjects <- unique(m_long$Subject)
num_unique_subjects <- length(unique_subjects)
print(num_unique_subjects)

# be sure that make non-exposed the ref group
m_long$group_factor <- factor(m_long$group_factor, levels = c("Non exposed", "Exposed to genocide", "Exposed to genocide and rape"))

write.csv(m_long, "data_cooked/m_long.csv")


## Descriptive and Correlation table


##candidate gene analysis


d <- subset(m_long, select = c("subject2", "Name", "M_Value", "group_factor", "bio_sex", "ace_total", "PC1", "UCSC_RefGene_Group", "deprtotal_promis29", "anxttotal_promis29"))

probes <- unique(d$Name)
n_tests <- length(probes) 

length(probes)


model_contrasts <- data.frame(term = c(),
                              estimate = c(),
                              std.error = c(),
                              statistic = c(),
                              p.value = c(),
                              conf.low = c(),
                              conf.high = c()) # init data frame

###ACEs exposure

for (i in 1:n_tests) {
  m_temp <- lm(M_Value ~ ace_total + bio_sex + PC1, data = dplyr::filter(d, Name == probes[i]))
  #contrast_groups <- avg_comparisons(m_temp, variables = "group_factor")
  summary(m_temp)$coefficients[2,]
  
  contrast_groups <- data.frame(
    term = "ace_total_plus1",
    estimate = summary(m_temp)$coefficients[2,1],
    std.error = summary(m_temp)$coefficients[2,2],
    statistic = summary(m_temp)$coefficients[2,3],
    p.value = summary(m_temp)$coefficients[2,4],
    conf.low = confint(m_temp, "ace_total")[1],
    conf.high = confint(m_temp, "ace_total")[2]
  )
  
  temp_df <- as.data.frame(contrast_groups) %>% mutate(probe = probes[i])
  model_contrasts <- bind_rows(model_contrasts, temp_df)
}

model_contrasts$significance <- ifelse(model_contrasts$p.value < 0.05, "p < 0.05", "p >= 0.05")

m_long_gene <- m_long %>% 
  group_by(Name) %>% 
  summarise(Gene = paste(unique(Gene_Name), collapse = "_"))

sig_contrasts <- model_contrasts %>% 
  left_join(m_long_gene, by = c("probe" = "Name")) %>% 
  filter(p.value < 0.05)

write.csv(sig_contrasts, "figures_tables/Significant_CpG_noFDR_ace.csv")


nrow(sig_contrasts) # total number sig contrasts

# significant contrasts by gene
sig_contrasts %>% 
  group_by(Gene) %>% 
  tally() %>% 
  arrange(desc(n))


## FDR correction

pvalues_adj <- p.adjust(model_contrasts$p.value, method = "fdr")
model_contrasts$significant_FDR <- pvalues_adj < 0.05

sig_contrasts_FDR <- model_contrasts %>% 
  filter(significant_FDR == T) %>% 
  left_join(m_long_gene, by = c("probe" = "Name"))

## save as table

write.csv(sig_contrasts_FDR, "figures_tables/Significant_CpG_FDR_ace.csv")


### War Trauma Exposure
model_contrasts <- data.frame(term = c(),
                              estimate = c(),
                              std.error = c(),
                              statistic = c(),
                              p.value = c(),
                              conf.low = c(),
                              conf.high = c()) # init data frame

for (i in 1:n_tests) {
  m_temp <- lm(M_Value ~ group_factor + bio_sex + PC1, data = dplyr::filter(d, Name == probes[i]))
  contrast_groups <- avg_comparisons(m_temp, variables = "group_factor")
  
  temp_df <- as.data.frame(contrast_groups) %>% mutate(probe = probes[i])
  model_contrasts <- bind_rows(model_contrasts, temp_df)
}

model_contrasts$significance <- ifelse(model_contrasts$p.value < 0.05, "p < 0.05", "p >= 0.05")

m_long_gene <- m_long %>% 
  group_by(Name) %>% 
  summarise(Gene = paste(unique(Gene_Name), collapse = "_"))

sig_contrasts <- model_contrasts %>% 
  left_join(m_long_gene, by = c("probe" = "Name")) %>% 
  filter(p.value < 0.05)

nrow(sig_contrasts) # total number sig contrasts

# significant contrasts by gene
sig_contrasts %>% 
  group_by(Gene) %>% 
  tally() %>% 
  arrange(desc(n))


write.csv(sig_contrasts, "figures_tables/Significant_CpG_noFDR_group.csv")

## FDR correction

pvalues_adj <- p.adjust(model_contrasts$p.value, method = "fdr")
model_contrasts$significant_FDR <- pvalues_adj < 0.05

sig_contrasts_FDR <- model_contrasts %>% 
  filter(significant_FDR == T) %>% 
  left_join(m_long_gene, by = c("probe" = "Name"))

## save as table

write.csv(sig_contrasts_FDR, "figures_tables/Significant_CpG_FDR_group.csv")




### War Trauma Exposure with aces as a covariate
model_contrasts <- data.frame(term = c(),
                              estimate = c(),
                              std.error = c(),
                              statistic = c(),
                              p.value = c(),
                              conf.low = c(),
                              conf.high = c()) # init data frame

for (i in 1:n_tests) {
  m_temp <- lm(M_Value ~ group_factor + ace_total+ bio_sex + PC1, data = dplyr::filter(d, Name == probes[i]))
  contrast_groups <- avg_comparisons(m_temp, variables = "group_factor")
  
  temp_df <- as.data.frame(contrast_groups) %>% mutate(probe = probes[i])
  model_contrasts <- bind_rows(model_contrasts, temp_df)
}

model_contrasts$significance <- ifelse(model_contrasts$p.value < 0.05, "p < 0.05", "p >= 0.05")

m_long_gene <- m_long %>% 
  group_by(Name) %>% 
  summarise(Gene = paste(unique(Gene_Name), collapse = "_"))

sig_contrasts <- model_contrasts %>% 
  left_join(m_long_gene, by = c("probe" = "Name")) %>% 
  filter(p.value < 0.05)

nrow(sig_contrasts) # total number sig contrasts

# significant contrasts by gene
sig_contrasts %>% 
  group_by(Gene) %>% 
  tally() %>% 
  arrange(desc(n))


write.csv(sig_contrasts, "figures_tables/Significant_CpG_noFDR_group_covar_ace.csv")

## FDR correction

pvalues_adj <- p.adjust(model_contrasts$p.value, method = "fdr")
model_contrasts$significant_FDR <- pvalues_adj < 0.05

sig_contrasts_FDR <- model_contrasts %>% 
  filter(significant_FDR == T) %>% 
  left_join(m_long_gene, by = c("probe" = "Name"))

## save as table

write.csv(sig_contrasts_FDR, "figures_tables/Significant_CpG_FDR_group_covar_ace.csv")





### Exploratory analyses on phenotypic correlations with methylation values on sig genes

### FKBP5 probe cg03245912

subset_d <- subset(d, Name == "cg03245912")

depression=lm(M_Value~deprtotal_promis29 + bio_sex +PC1, data= subset_d)
summary(depression1)

anxiety=lm(M_Value~anxttotal_promis29 +bio_sex +PC1, data= subset_d)
summary(anxiety1)

###BDNF probe cg06979684
subset_d <- subset(d, Name == "cg06979684")

depression=lm(M_Value~deprtotal_promis29 + bio_sex +PC1, data= subset_d)
summary(depression2)

anxiety=lm(M_Value~anxttotal_promis29 +bio_sex +PC1, data= subset_d)
summary(anxiety2)

###SLC6A4 probe cg26438554
subset_d <- subset(d, Name == "cg26438554")

depression=lm(M_Value~deprtotal_promis29 + bio_sex +PC1, data= subset_d)
summary(depression3)

anxiety=lm(M_Value~anxttotal_promis29 +bio_sex +PC1, data= subset_d)
summary(anxiety3)

depression1 <- as_flextable(depression1)  
anxiety1 <- as_flextable(anxiety1)
depression2 <- as_flextable(depression2)  
anxiety2 <- as_flextable(anxiety2)
depression3 <- as_flextable(depression3)  
anxiety3 <- as_flextable(anxiety3)

save_as_docx(
  `cg03245912 Methylation ~ Promise 29 Depression Total` = depression1,  `cg03245912 Methylation ~ Promise 29 Anxiety Total` = anxiety1,  `cg06979684 Methylation ~ Promise 29 Depression Total`=depression2,`cg06979684 Methylation ~ Promise 29 AnxietyTotal`=anxiety2, `cg26438554 Methylation ~ Promise 29 Depression Total`=depression3,`cg26438554 Methylation ~ Promise 29 Anxiety Total`=anxiety3, path ="figures_tables/depression_anxiety.docx")

