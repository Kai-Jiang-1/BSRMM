# Construct Tree from GTBD
rm(list = ls())
setwd("~/Desktop/BSRMM/realdata/data_setup")

# required pacakge
library(ape)
library(data.table)
library(tidyverse)

# import archaea
# tree 
ar_tree <- read.tree("ar53_r207.tree")

# taxonomy 
ar_tax <- fread("ar53_taxonomy_r207.tsv", sep = "\t", header = FALSE, data.table = FALSE)
ar_tax <- ar_tax %>%
  filter(V1 %in% ar_tree$tip.label)

# reorder ar_tax based on the ar_tree$tip.label
all(ar_tree$tip.label == ar_tax$V1)
ar_tax <- ar_tax[match(ar_tree$tip.label, ar_tax$V1), ]
all(ar_tree$tip.label == ar_tax$V1)

# rename the tree
ar_tree$tip.label <- ar_tax$V2

# import bacteria
# tree
bac_tree <- read.tree("bac120_r207.tree")

# taxonomy
bac_tax <- fread("bac120_taxonomy_r207.tsv", sep = "\t", header = FALSE, data.table = FALSE)
bac_tax <- bac_tax %>% 
  filter(V1 %in% bac_tree$tip.label)

# reorder bac_tax based on the bac_tree$tip.label
all(bac_tree$tip.label == bac_tax$V1)
bac_tax <- bac_tax[match(bac_tree$tip.label, bac_tax$V1), ]
all(bac_tree$tip.label == bac_tax$V1)

# rename the tree 
bac_tree$tip.label <- bac_tax$V2

# merge ar_tree and bac_tree together
tree <- bind.tree(ar_tree, bac_tree)
write.tree(tree, "gtbd.tree")
