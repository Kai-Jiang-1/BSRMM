# data processing and cleaning 
rm(list=ls())

# set working directory
setwd("~/Desktop/BSRMM/realdata/data_setup")

# packages
library(phyloseq)
library(ape)
library(ggtree)
library(pheatmap)
library(data.table)
library(stringr)
library(ppcor)
library(phyloseq)
library(microbiome)

# import data set 
# meta data 
meta <- fread("metadata.tsv", data.table = F)

# microbiome data 
otu.tab <- fread("species.tsv", data.table = F) 

# metabolite data 
metab <- fread("mtb.tsv", data.table = F)
sum(metab$C00378_Thiamine ==0)/nrow(metab)
# 0.3342939

# check whether the data is in the same order
# meta vs otu.tab 
all(meta$Sample == otu.tab$Sample)
all(sort(meta$Sample) == sort(otu.tab$Sample))
meta <- meta[order(meta$Sample),]
otu.tab <- otu.tab[order(otu.tab$Sample),]
all(meta$Sample == otu.tab$Sample)

# meta vs metab 
all(meta$Sample == metab$Sample)
all(meta$Sample == sort(metab$Sample))
metab <- metab[order(metab$Sample),]
all(meta$Sample == metab$Sample)

# otu.tab vs metab
all(otu.tab$Sample == metab$Sample)
 
# convert the sample to rownames 
rownames(meta) <- meta$Sample
rownames(otu.tab) <- otu.tab$Sample
rownames(metab) <- metab$Sample

# remove the sample column for otu.tab
otu.tab$Sample <- NULL

# import tree information 
tree <- read.tree("gtbd.tree")
tree$tip.label[1:10]

# construct taxanomic table 
otu.name <- data.frame(name = colnames(otu.tab), stringsAsFactors = F)
otu.name <- str_split(otu.name$name, ";")
otu.name <- do.call(rbind, otu.name)
colnames(otu.name) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")
otu.name <- data.frame(otu.name, stringsAsFactors = F)
otu.name$Species <- str_replace_all(otu.name$Species, " ", "_")

# match the name between OTU vs referecne tree
name_list <- apply(otu.name, 1, paste, collapse = "-")

# check if all name from name_list are contained in the tree$tip.label
all(name_list %in% tree$tip.label)

# covert to colnames for otu.tab
colnames(otu.tab) <- name_list

# covert to rownames for otu.name 
rownames(otu.name) <- name_list

# filter tree based on otu name 
tree <- drop.tip(tree, tree$tip.label[!(tree$tip.label %in% name_list)])

# create the pholyseq object 
meta <- sample_data(meta)
otu.tab <- otu_table(as.matrix(otu.tab), taxa_are_rows = FALSE)
otu.name <- tax_table(as.matrix(otu.name))

myphylo <- phyloseq(meta, otu.tab, otu.name,tree)

# check the order to match 
all(rownames(myphylo@otu_table) == rownames(myphylo@sam_data))
all(colnames(myphylo@otu_table) == rownames(myphylo@tax_table))
all(rownames(myphylo@tax_table) == myphylo@phy_tree$tip.label)
all(myphylo@sam_data$Sample == metab$Sample)

# data cleaning
# proportion of presence for each otu
prop_nonzero <- colMeans(myphylo@otu_table > 0)

# keep otu which are present in at least 30% of the samples
length(which(prop_nonzero < 0.3))
myphylo <- prune_taxa(prop_nonzero > 0.3, myphylo)
dim(myphylo@otu_table)

# drop if average relative abundance less than 0.03% 
myphylo <- prune_taxa(colMeans(myphylo@otu_table) > 0.0003, myphylo)
myphylo <- transform(myphylo, "compositional")
dim(myphylo@otu_table)
prop_zero <- colMeans(myphylo@otu_table == 0)
quantile(prop_zero)

# construct Q matrix.
# step 1:
# phylogenetic correlation matrix 
q_corr <- vcv(myphylo@phy_tree, corr = T)

# step 2: 
# inverse phylogenetic correlation matrix
q_inv <- solve(q_corr)

# step 3:
# partial correlation matrix
q_corr_par <- -q_inv / sqrt(outer(diag(q_inv), diag(q_inv)))

# heatmap
pheatmap(q_corr_par,
         cluster_rows = F,
         cluster_cols = F,
         show_rownames = F,
         main = " ",
         show_colnames = F,
         fontsize = 12,
         border_color = NA,
         color = colorRampPalette(c("#f0f0f0", "#85C0F9", "#003f5c"))(100))

# set diagonal to 0 
diag(q_corr_par) <- 0

# check the negative values 
sum(q_corr_par < 0)
sum(q_corr_par)

# we set all negative values to 0 
q_corr_par[q_corr_par < 0] <- 0
sum(q_corr_par)

# re-do the heatmap 
pheatmap(q_corr_par,
         filename = "Partial_correlation.pdf",
         cluster_rows = F,
         cluster_cols = F,
         show_rownames = F,
         main = " ",
         show_colnames = F,
         fontsize = 12,
         border_color = NA,
         color = colorRampPalette(c("#f0f0f0", "#85C0F9", "#003f5c"))(100),
         breaks = seq(0, 1, length.out = 101), 
         legend = TRUE,
         legend_breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1),
         width = 10,
         height = 10)

# obtain the X matrix
X <- as.matrix(myphylo@otu_table)
X <- as(X, "matrix")

X <- apply(X, 2, function(x) {
  min_nonzero <- min(x[x > 0])
  impute_nonzero <- 0.5 * min_nonzero
  x[x == 0] <- impute_nonzero
  return(x)
})

# redo the compostional
X <- t(apply(X, 1, function(x) {
  x <- x / sum(x)
  return(x)
}))

sum(rowSums(X))

# save the dataset
mydata <- list(Y_thi = metab$`C00378_Thiamine`,
               X = X,
               group = meta$Study.Group,
               q_corr_par = q_corr_par, # partial correlation matrix
               q_corr = q_corr, # phylogenetic correlation matrix
               q_zero = diag(0, nrow = ncol(X)), # zero matrix
               tax_info = myphylo@tax_table)

saveRDS(mydata, "mydata_relative_abundance.rds")
