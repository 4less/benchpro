library(ggplot2)
library(dplyr)
library(tidyr)
library(ggpubr)

msa.summary <- read.csv("/usr/users/QIB_fr017/fritsche/Projects/benchpro/data/SPECIES46/MSSS200C/msa.summary2.tsv", sep='\t', header=FALSE)
colnames(msa.summary) <- c("Tool", "Species", "Genome", "GroupSize", "SNPs", "AlignmentLength", "ErrorRate", "MinCov", "WholeGroup")
msa.summary$WholeGroup <- msa.summary$WholeGroup == "True"
typeof(msa.summary$MinCov)
msa.summary$MinCovInt <- as.integer(msa.summary$MinCov)
msa.summary$MinCovInt <- as.integer(msa.summary$MinCov/3) * 3
msa.summary <- msa.summary %>% filter(WholeGroup)

msa.summary %>% ggplot(aes(x=Species, y=SNPs, fill=Tool)) +
  geom_boxplot() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

msa.summary %>% ggplot(aes(x=MinCov, y=SNPs, fill=Tool, color=Tool)) +
  geom_point() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  facet_wrap(~ Species, nrow=3, scales="free_y")

msa.summary %>% ggplot(aes(x=MinCov, y=SNPs, fill=Tool, color=Tool)) +
  geom_line() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))#

msa.summary %>%
  mutate(MinCovInt=as.integer(msa.summary$MinCov/3) * 3) %>%
  ggplot(aes(x=MinCovInt, y=SNPs, fill=Tool, color=Tool)) +
  stat_summary(fun.y=mean, geom="line") + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
  facet_wrap(~ Species, nrow=3, scales="free_y")

exclude_species <- c("s__Bacteroides_xylanisolvens", "s__Blautia_A_sp003471165")
msa.summary %>%
  filter(!(Species %in% exclude_species)) %>%
  mutate(MinCovInt=as.integer(MinCov/1) * 1) %>%
  ggplot(aes(x=MinCovInt, y=SNPs, fill=Tool, color=Tool)) +
  stat_summary(fun.y=mean, geom="line") + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

msa.summary %>%
  pivot_longer(cols=c("ErrorRate", "SNPs", "GroupSize"), names_to="Metric", values_to="Value") %>% 
  ggplot(aes(x=Species, y=Value, fill=Tool)) +
  geom_boxplot() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
  facet_wrap(~ Metric, nrow=3, scales="free_y") +
  scale_y_continuous(trans=scales::pseudo_log_trans(base = 10)) + #theme_classic()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))


msa.summary$SNPs2 <- msa.summary$SNPs 
msa.summary$Group2 <- msa.summary$Group
msa.summary$ErrorRate2 <- msa.summary$ErrorRate
msa.summary$SNPs2[msa.summary$SNPs == 0] <- 1e-2
msa.summary$ErrorRate2[msa.summary$ErrorRate == 0] <- 1e-6


msa.summary %>%
  ggplot(aes(x=Tool, y=ErrorRate, fill=Tool)) +
  geom_boxplot() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
  facet_wrap(~ Species, nrow=3, scales="free_y") +
  scale_y_continuous(trans=scales::pseudo_log_trans(base = 10)) + #theme_classic()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))


msa.summary %>% 
  pivot_longer(cols=c("AlignmentLength", "GroupSize", "ErrorRate"), names_to="Metric", values_to="Value") %>%
  ggplot(aes(x=Tool, y=Value, fill=Tool)) +
  geom_boxplot() + 
  facet_wrap(~ Metric, ncol=3, scales="free_y") + 
  #scale_y_continuous(trans=scales::pseudo_log_trans(base = 10)) + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
