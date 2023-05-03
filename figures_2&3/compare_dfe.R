library(ggpubr)
library(tidyverse)
library(outliers)


#### Load data
relFitness <- read.csv("relative_fitness.tsv", sep="\t")
relFitness_reps <- read.csv("relative_fitness_reps.csv", sep=";")



#### DFE

relFitness$strain_name <- factor(relFitness$strain_name, levels = relFitness$strain_name)

ggplot(relFitness[relFitness$source == "This work", ], aes(x=strain_name, y=w_norm)) +
  geom_segment(aes(x=strain_name, xend=strain_name, y=1, yend=w_norm, color=species), size=5, alpha=0.5) +
  theme_classic() +
  geom_errorbar(aes(ymin=w_norm-SE_prop, ymax=w_norm+SE_prop), width=0, color="black") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))



#### Normality tests

# checking outliers
grubbs.test(relFitness[relFitness$source == "This work", ]$w_norm)

relFitness_mod <- relFitness
relFitness_mod[relFitness_mod=="C_freundii"] <- "Other"
relFitness_mod[relFitness_mod=="K_ascorbata"] <- "Other"
relFitness_mod[relFitness_mod=="K_variicola"] <- "Other"

# by species
shapiro.test(relFitness_mod[relFitness_mod$source == "Alonso-del Valle et al. 2021" & relFitness_mod$species == "E_coli", ]$w_norm)
shapiro.test(relFitness_mod[relFitness_mod$source == "Alonso-del Valle et al. 2021" & relFitness_mod$species == "Klebsiella_spp", ]$w_norm)
shapiro.test(relFitness_mod[relFitness_mod$source == "This work" & relFitness_mod$species == "E_coli", ]$w_norm)
shapiro.test(relFitness_mod[relFitness_mod$source == "This work" & relFitness_mod$species == "Klebsiella_pneumoniae", ]$w_norm)
shapiro.test(relFitness_mod[relFitness_mod$source == "This work" & relFitness_mod$species == "Other", ]$w_norm)

# all strains
shapiro.test(relFitness_mod[relFitness_mod$source == "Alonso-del Valle et al. 2021", ]$w_norm)
shapiro.test(relFitness_mod[relFitness_mod$source == "This work", ]$w_norm)
relFitness_woN46 <- subset(relFitness, strain_name!="AJ_N46")
shapiro.test(relFitness_woN46[relFitness_woN46$source == "This work", ]$w_norm)

# homoscedasticity test
bartlett.test(w_norm ~ source, relFitness_mod)



#### Comparing means

# between species of the same collection
t.test(relFitness_mod[relFitness_mod$source == "Alonso-del Valle et al. 2021" & relFitness_mod$species == "E_coli", ]$w_norm, relFitness_mod[relFitness_mod$source == "Alonso-del Valle et al. 2021"  & relFitness_mod$species == "Klebsiella_spp", ]$w_norm)
wilcox.test(relFitness_mod[relFitness_mod$source == "This work" & relFitness_mod$species == "E_coli", ]$w_norm, relFitness_mod[relFitness_mod$source == "This work"  & relFitness_mod$species == "Klebsiella_pneumoniae", ]$w_norm)
t.test(relFitness_woN46[relFitness_woN46$source == "This work" & relFitness_woN46$species == "E_coli", ]$w_norm, relFitness_woN46[relFitness_woN46$source == "This work"  & relFitness_woN46$species == "Klebsiella_pneumoniae", ]$w_norm)
kruskal.test(relFitness_mod[relFitness_mod$source == "This work", ]$w_norm, relFitness_mod[relFitness_mod$source == "This work", ]$species)
# between collections
t.test(relFitness_mod[relFitness_mod$source == "Alonso-del Valle et al. 2021", ]$w_norm, relFitness_mod[relFitness_mod$source == "This work", ]$w_norm)
ks.test(relFitness_mod[relFitness_mod$source == "Alonso-del Valle et al. 2021", ]$w_norm, relFitness_mod[relFitness_mod$source == "This work", ]$w_norm)


# among strains of this collection
table_pvalues <- data.frame(Strain=character(0), pvalue=numeric(0))

table_pvalues <- table_pvalues %>% add_row (Strain= "AJ_N46",
                                            pvalue= t.test(relFitness_reps[relFitness_reps$Strain == "AJ_N46" & relFitness_reps$Genotype == "WT", ]$Rel_fitness,
                                                           relFitness_reps[relFitness_reps$Strain == "AJ_N46" & relFitness_reps$Genotype == "pOXA-48_free", ]$Rel_fitness,
                                                           paired=TRUE)$p.value)
table_pvalues <- table_pvalues %>% add_row (Strain= "AJ_C310",
                                            pvalue= t.test(relFitness_reps[relFitness_reps$Strain == "AJ_C310" & relFitness_reps$Genotype == "WT", ]$Rel_fitness,
                                                           relFitness_reps[relFitness_reps$Strain == "AJ_C310" & relFitness_reps$Genotype == "pOXA-48_free", ]$Rel_fitness,
                                                           paired=TRUE)$p.value)
table_pvalues <- table_pvalues %>% add_row (Strain= "AJ_C728",
                                            pvalue= t.test(relFitness_reps[relFitness_reps$Strain == "AJ_C728" & relFitness_reps$Genotype == "WT", ]$Rel_fitness,
                                                           relFitness_reps[relFitness_reps$Strain == "AJ_C728" & relFitness_reps$Genotype == "pOXA-48_free", ]$Rel_fitness,
                                                           paired=TRUE)$p.value)
table_pvalues <- table_pvalues %>% add_row (Strain= "CF12",
                                            pvalue= t.test(relFitness_reps[relFitness_reps$Strain == "CF12" & relFitness_reps$Genotype == "WT", ]$Rel_fitness,
                                                           relFitness_reps[relFitness_reps$Strain == "CF12" & relFitness_reps$Genotype == "pOXA-48_free", ]$Rel_fitness,
                                                           paired=TRUE)$p.value)
table_pvalues <- table_pvalues %>% add_row (Strain= "K163",
                                            pvalue= t.test(relFitness_reps[relFitness_reps$Strain == "K163" & relFitness_reps$Genotype == "WT", ]$Rel_fitness,
                                                           relFitness_reps[relFitness_reps$Strain == "K163" & relFitness_reps$Genotype == "pOXA-48_free", ]$Rel_fitness,
                                                           paired=TRUE)$p.value)
table_pvalues <- table_pvalues %>% add_row (Strain= "AJ_Nh27",
                                            pvalue= t.test(relFitness_reps[relFitness_reps$Strain == "AJ_Nh27" & relFitness_reps$Genotype == "WT", ]$Rel_fitness,
                                                           relFitness_reps[relFitness_reps$Strain == "AJ_Nh27" & relFitness_reps$Genotype == "pOXA-48_free", ]$Rel_fitness,
                                                           paired=TRUE)$p.value)
table_pvalues <- table_pvalues %>% add_row (Strain= "R10",
                                            pvalue= t.test(relFitness_reps[relFitness_reps$Strain == "R10" & relFitness_reps$Genotype == "WT", ]$Rel_fitness,
                                                           relFitness_reps[relFitness_reps$Strain == "R10" & relFitness_reps$Genotype == "pOXA-48_free", ]$Rel_fitness,
                                                           paired=TRUE)$p.value)
table_pvalues <- table_pvalues %>% add_row (Strain= "AJ_C527",
                                            pvalue= t.test(relFitness_reps[relFitness_reps$Strain == "AJ_C527" & relFitness_reps$Genotype == "WT", ]$Rel_fitness,
                                                           relFitness_reps[relFitness_reps$Strain == "AJ_C527" & relFitness_reps$Genotype == "pOXA-48_free", ]$Rel_fitness,
                                                           paired=TRUE)$p.value)
table_pvalues <- table_pvalues %>% add_row (Strain= "AJ_N23",
                                            pvalue= t.test(relFitness_reps[relFitness_reps$Strain == "AJ_N23" & relFitness_reps$Genotype == "WT", ]$Rel_fitness,
                                                           relFitness_reps[relFitness_reps$Strain == "AJ_N23" & relFitness_reps$Genotype == "pOXA-48_free", ]$Rel_fitness,
                                                           paired=TRUE)$p.value)
table_pvalues <- table_pvalues %>% add_row (Strain= "AJ_K219",
                                            pvalue= t.test(relFitness_reps[relFitness_reps$Strain == "AJ_K219" & relFitness_reps$Genotype == "WT", ]$Rel_fitness,
                                                           relFitness_reps[relFitness_reps$Strain == "AJ_K219" & relFitness_reps$Genotype == "pOXA-48_free", ]$Rel_fitness,
                                                           paired=TRUE)$p.value)
table_pvalues <- table_pvalues %>% add_row (Strain= "AJ_K198",
                                            pvalue= t.test(relFitness_reps[relFitness_reps$Strain == "AJ_K198" & relFitness_reps$Genotype == "WT", ]$Rel_fitness,
                                                           relFitness_reps[relFitness_reps$Strain == "AJ_K198" & relFitness_reps$Genotype == "pOXA-48_free", ]$Rel_fitness,
                                                           paired=TRUE)$p.value)
table_pvalues <- table_pvalues %>% add_row (Strain= "AJ_L38",
                                            pvalue= t.test(relFitness_reps[relFitness_reps$Strain == "AJ_L38" & relFitness_reps$Genotype == "WT", ]$Rel_fitness,
                                                           relFitness_reps[relFitness_reps$Strain == "AJ_L38" & relFitness_reps$Genotype == "pOXA-48_free", ]$Rel_fitness,
                                                           paired=TRUE)$p.value)
table_pvalues <- table_pvalues %>% add_row (Strain= "C642",
                                            pvalue= t.test(relFitness_reps[relFitness_reps$Strain == "C642" & relFitness_reps$Genotype == "WT", ]$Rel_fitness,
                                                           relFitness_reps[relFitness_reps$Strain == "C642" & relFitness_reps$Genotype == "pOXA-48_free", ]$Rel_fitness,
                                                           paired=TRUE)$p.value)
table_pvalues <- table_pvalues %>% add_row (Strain= "AJ_S8",
                                            pvalue= t.test(relFitness_reps[relFitness_reps$Strain == "AJ_S8" & relFitness_reps$Genotype == "WT", ]$Rel_fitness,
                                                           relFitness_reps[relFitness_reps$Strain == "AJ_S8" & relFitness_reps$Genotype == "pOXA-48_free", ]$Rel_fitness,
                                                           paired=TRUE)$p.value)
table_pvalues <- table_pvalues %>% add_row (Strain= "C325",
                                            pvalue= t.test(relFitness_reps[relFitness_reps$Strain == "C325" & relFitness_reps$Genotype == "WT", ]$Rel_fitness,
                                                           relFitness_reps[relFitness_reps$Strain == "C325" & relFitness_reps$Genotype == "pOXA-48_free", ]$Rel_fitness,
                                                           paired=TRUE)$p.value)
table_pvalues <- table_pvalues %>% add_row (Strain= "C609",
                                            pvalue= t.test(relFitness_reps[relFitness_reps$Strain == "C609" & relFitness_reps$Genotype == "WT", ]$Rel_fitness,
                                                           relFitness_reps[relFitness_reps$Strain == "C609" & relFitness_reps$Genotype == "pOXA-48_free", ]$Rel_fitness,
                                                           paired=TRUE)$p.value)
table_pvalues <- table_pvalues %>% add_row (Strain= "AJ_K25",
                                            pvalue= t.test(relFitness_reps[relFitness_reps$Strain == "AJ_K25" & relFitness_reps$Genotype == "WT", ]$Rel_fitness,
                                                           relFitness_reps[relFitness_reps$Strain == "AJ_K25" & relFitness_reps$Genotype == "pOXA-48_free", ]$Rel_fitness,
                                                           paired=TRUE)$p.value)
table_pvalues <- table_pvalues %>% add_row (Strain= "K147",
                                            pvalue= t.test(relFitness_reps[relFitness_reps$Strain == "K147" & relFitness_reps$Genotype == "WT", ]$Rel_fitness,
                                                           relFitness_reps[relFitness_reps$Strain == "K147" & relFitness_reps$Genotype == "pOXA-48_free", ]$Rel_fitness,
                                                           paired=TRUE)$p.value)
table_pvalues <- table_pvalues %>% add_row (Strain= "C662",
                                            pvalue= t.test(relFitness_reps[relFitness_reps$Strain == "C662" & relFitness_reps$Genotype == "WT", ]$Rel_fitness,
                                                           relFitness_reps[relFitness_reps$Strain == "C662" & relFitness_reps$Genotype == "pOXA-48_free", ]$Rel_fitness,
                                                           paired=TRUE)$p.value)
table_pvalues <- table_pvalues %>% add_row (Strain= "CF13",
                                            pvalue= t.test(relFitness_reps[relFitness_reps$Strain == "CF13" & relFitness_reps$Genotype == "WT", ]$Rel_fitness,
                                                           relFitness_reps[relFitness_reps$Strain == "CF13" & relFitness_reps$Genotype == "pOXA-48_free", ]$Rel_fitness,
                                                           paired=TRUE)$p.value)
table_pvalues <- table_pvalues %>% add_row (Strain= "H53",
                                            pvalue= t.test(relFitness_reps[relFitness_reps$Strain == "H53" & relFitness_reps$Genotype == "WT", ]$Rel_fitness,
                                                           relFitness_reps[relFitness_reps$Strain == "H53" & relFitness_reps$Genotype == "pOXA-48_free", ]$Rel_fitness,
                                                           paired=TRUE)$p.value)
table_pvalues <- table_pvalues %>% add_row (Strain= "AJ_K88",
                                            pvalue= t.test(relFitness_reps[relFitness_reps$Strain == "AJ_K88" & relFitness_reps$Genotype == "WT", ]$Rel_fitness,
                                                           relFitness_reps[relFitness_reps$Strain == "AJ_K88" & relFitness_reps$Genotype == "pOXA-48_free", ]$Rel_fitness,
                                                           paired=TRUE)$p.value)
table_pvalues <- table_pvalues %>% add_row (Strain= "AJ_Z29",
                                            pvalue= t.test(relFitness_reps[relFitness_reps$Strain == "AJ_Z29" & relFitness_reps$Genotype == "WT", ]$Rel_fitness,
                                                           relFitness_reps[relFitness_reps$Strain == "AJ_Z29" & relFitness_reps$Genotype == "pOXA-48_free", ]$Rel_fitness,
                                                           paired=TRUE)$p.value)
table_pvalues <- table_pvalues %>% add_row (Strain= "AJ_K127",
                                            pvalue= t.test(relFitness_reps[relFitness_reps$Strain == "AJ_K127" & relFitness_reps$Genotype == "WT", ]$Rel_fitness,
                                                           relFitness_reps[relFitness_reps$Strain == "AJ_K127" & relFitness_reps$Genotype == "pOXA-48_free", ]$Rel_fitness,
                                                           paired=TRUE)$p.value)
table_pvalues <- table_pvalues %>% add_row (Strain= "AJ_K293",
                                            pvalue= t.test(relFitness_reps[relFitness_reps$Strain == "AJ_K293" & relFitness_reps$Genotype == "WT", ]$Rel_fitness,
                                                           relFitness_reps[relFitness_reps$Strain == "AJ_K293" & relFitness_reps$Genotype == "pOXA-48_free", ]$Rel_fitness,
                                                           paired=TRUE)$p.value)
table_pvalues <- table_pvalues %>% add_row (Strain= "AJ_K308",
                                            pvalue= t.test(relFitness_reps[relFitness_reps$Strain == "AJ_K308" & relFitness_reps$Genotype == "WT", ]$Rel_fitness,
                                                           relFitness_reps[relFitness_reps$Strain == "AJ_K308" & relFitness_reps$Genotype == "pOXA-48_free", ]$Rel_fitness,
                                                           paired=TRUE)$p.value)
table_pvalues <- table_pvalues %>% add_row (Strain= "AJ_K244",
                                            pvalue= t.test(relFitness_reps[relFitness_reps$Strain == "AJ_K244" & relFitness_reps$Genotype == "WT", ]$Rel_fitness,
                                                           relFitness_reps[relFitness_reps$Strain == "AJ_K244" & relFitness_reps$Genotype == "pOXA-48_free", ]$Rel_fitness,
                                                           paired=TRUE)$p.value)
table_pvalues <- table_pvalues %>% add_row (Strain= "AJ_K318",
                                            pvalue= t.test(relFitness_reps[relFitness_reps$Strain == "AJ_K318" & relFitness_reps$Genotype == "WT", ]$Rel_fitness,
                                                           relFitness_reps[relFitness_reps$Strain == "AJ_K318" & relFitness_reps$Genotype == "pOXA-48_free", ]$Rel_fitness,
                                                           paired=TRUE)$p.value)
table_pvalues <- table_pvalues %>% add_row (Strain= "AJ_C646",
                                            pvalue= t.test(relFitness_reps[relFitness_reps$Strain == "AJ_C646" & relFitness_reps$Genotype == "WT", ]$Rel_fitness,
                                                           relFitness_reps[relFitness_reps$Strain == "AJ_C646" & relFitness_reps$Genotype == "pOXA-48_free", ]$Rel_fitness,
                                                           paired=TRUE)$p.value)
table_pvalues <- table_pvalues %>% add_row (Strain= "K153",
                                            pvalue= t.test(relFitness_reps[relFitness_reps$Strain == "K153" & relFitness_reps$Genotype == "WT", ]$Rel_fitness,
                                                           relFitness_reps[relFitness_reps$Strain == "K153" & relFitness_reps$Genotype == "pOXA-48_free", ]$Rel_fitness,
                                                           paired=TRUE)$p.value)
table_pvalues <- table_pvalues %>% add_row (Strain= "J57",
                                            pvalue= t.test(relFitness_reps[relFitness_reps$Strain == "J57" & relFitness_reps$Genotype == "WT", ]$Rel_fitness,
                                                           relFitness_reps[relFitness_reps$Strain == "J57" & relFitness_reps$Genotype == "pOXA-48_free", ]$Rel_fitness,
                                                           paired=TRUE)$p.value)
table_pvalues <- table_pvalues %>% add_row (Strain= "AJ_J61",
                                            pvalue= t.test(relFitness_reps[relFitness_reps$Strain == "AJ_J61" & relFitness_reps$Genotype == "WT", ]$Rel_fitness,
                                                           relFitness_reps[relFitness_reps$Strain == "AJ_J61" & relFitness_reps$Genotype == "pOXA-48_free", ]$Rel_fitness,
                                                           paired=TRUE)$p.value)
table_pvalues <- table_pvalues %>% add_row (Strain= "AJ_C163",
                                            pvalue= t.test(relFitness_reps[relFitness_reps$Strain == "AJ_C163" & relFitness_reps$Genotype == "WT", ]$Rel_fitness,
                                                           relFitness_reps[relFitness_reps$Strain == "AJ_C163" & relFitness_reps$Genotype == "pOXA-48_free", ]$Rel_fitness,
                                                           paired=TRUE)$p.value)
table_pvalues <- table_pvalues %>% add_row (Strain= "AJ_C164",
                                            pvalue= t.test(relFitness_reps[relFitness_reps$Strain == "AJ_C164" & relFitness_reps$Genotype == "WT", ]$Rel_fitness,
                                                           relFitness_reps[relFitness_reps$Strain == "AJ_C164" & relFitness_reps$Genotype == "pOXA-48_free", ]$Rel_fitness,
                                                           paired=TRUE)$p.value)
table_pvalues <- table_pvalues %>% add_row (Strain= "C288",
                                            pvalue= t.test(relFitness_reps[relFitness_reps$Strain == "C288" & relFitness_reps$Genotype == "WT", ]$Rel_fitness,
                                                           relFitness_reps[relFitness_reps$Strain == "C288" & relFitness_reps$Genotype == "pOXA-48_free", ]$Rel_fitness,
                                                           paired=TRUE)$p.value)

table_pvalues$p.adj <- p.adjust(table_pvalues$pvalue, method="bonferroni", n=35)
table_pvalues[table_pvalues$p.adj < 0.05, ]


#### Q-Q plots

ggqqplot(relFitness_mod, x = "w_norm",
         color = "source", 
         palette = c("gray", "red3")
)



#### Histogram and density plots

phist <- gghistogram(relFitness_mod, x="w_norm", fill="source", add="mean", add.params=list(size=1, linetype=1),
                     palette = c("gray70", "red3"), bins=15
)
phist




#### Cumulative distribution function plots

ggplot(relFitness_mod[relFitness_mod$source == "Alonso-del Valle et al. 2021", ], aes(x=w_norm)) + 
  stat_ecdf(aes(x=w_norm, col=species)) + stat_ecdf(size=0.7) + 
  theme_classic()

ggplot(relFitness_mod[relFitness_mod$source == "This work", ], aes(x=w_norm)) + 
  stat_ecdf(aes(x=w_norm, col=species)) + stat_ecdf(size=0.7) +
  theme_classic()
 
