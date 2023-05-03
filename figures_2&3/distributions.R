library(tidyverse)
library(ggplot2)
library(ggpubr)

# Load data
total_collection <- read.table("total_collection.csv", header=T, sep="\t")
total_collection$Species_category <- factor(
  total_collection$Species_category, 
  levels=rev(sort(unique(total_collection$Species_category))), 
  ordered=TRUE
)

selected <- read.table("selected.csv", header=T, sep="\t")
selected$Species_category <- factor(
  selected$Species_category, 
  levels=rev(sort(unique(selected$Species_category))), 
  ordered=TRUE
)

cured <- read.table("cured.csv", header=T, sep="\t")
cured$Species_category <- factor(
  cured$Species_category, 
  levels=rev(sort(unique(cured$Species_category))), 
  ordered=TRUE
)


# Set color palette
Enterobacteria_n_ST <- total_collection %>% group_by(ST) %>% summarise(n_ST = n())
palette_ST = setNames(object = scales::hue_pal()(49), nm = Enterobacteria_n_ST$ST)
print(palette_ST)


# Plots

all_enterobacteria_distribution <-
  ggplot(total_collection, aes(fill=ST, x=Species_category)) +
  geom_bar(position="fill", stat="count", color="black", alpha=0.6) +
  theme(legend.position="none", axis.text.x=element_text(color="black", size=8, angle=90, vjust=0.5, hjust=1)) +
  theme_classic() +
  theme(legend.position="none") +
  labs(y="Frequency", x="Species", title="Total enterobacteria (n = 225)") +
  coord_flip() +
  scale_fill_manual(values=palette_ST)

selected_enterobacteria_distribution <-
  ggplot(selected, aes(fill=ST, x=Species_category)) +
  geom_bar(position="fill", stat="count", color="black", alpha=0.6) +
  theme(legend.position="none", axis.text.x=element_text(color="black", size=8, angle=90, vjust=0.5, hjust=1)) +
  theme_classic() +
  theme(legend.position="none") +
  labs(y="Frequency", x="Species", title="Selected enterobacteria (n = 60)") +
  coord_flip() +
  scale_fill_manual(values=palette_ST)

cured_enterobacteria_distribution <-
  ggplot(cured, aes(fill=ST, x=Species_category)) +
  geom_bar(position="fill", stat="count", color="black", alpha=0.6) +
  theme(legend.position="none", axis.text.x=element_text(color="black", size=8, angle=90, vjust=0.5, hjust=1)) +
  theme_classic() +
  theme(legend.position="none") +
  labs(y="Frequency", x="Species", title="Cured enterobacteria (n = 35)") +
  coord_flip() +
  scale_fill_manual(values=palette_ST)


# Final plot
ggarrange(all_enterobacteria_distribution, selected_enterobacteria_distribution, cured_enterobacteria_distribution,
          labels = c("(a)", "(b)", "(c)"),
          ncol = 1, nrow = 3, align = "v")

