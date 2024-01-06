library(phyloseq)
library(ggplot2)
library(plyr)
library(vegan)
library(agricolae)
library(FSA)
library(rcompanion)
library(boot)

otus <- read.csv("otu_table.csv", row.names=1)
otumat <- sample_data(otus)
otumat
taxmat <- read.csv("tax_table.csv", row.names=1)
colnames(taxmat) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus")
taxmat

meta_data = read.csv("sample_info.csv", header=T, row.names=1)
meta_data = sample_data(meta_data)
meta_data

#check that rownames for each table match
all(rownames(otumat) == rownames(taxmat))
#check that colnames match meta data
all(colnames(otumat) == unlist(meta_data[,'Sample_ID']))

OTU = otu_table(as.matrix(otumat), taxa_are_rows=T)
TAX = tax_table(as.matrix(taxmat))
OTU
TAX

physeq = phyloseq(OTU, TAX)
physeq

#combine meta data with physeq object
sample_data(physeq) = meta_data
physeq

metadata <- as(sample_data(physeq), "data.frame")

data_phylo_filt = filter_taxa(physeq, function(x) sum(x > 2) > (0.11 * length(x)), T)
set.seed(1425)
OTU_filt_rar = rarefy_even_depth(otu_table(data_phylo_filt), rngseed = T, replace = F)
data_otu_filt_rar = data.frame(otu_table(OTU_filt_rar))
data_phylo_filt_rar <- phyloseq(OTU_filt_rar, TAX, meta_data)

dist_bc <- as.matrix(vegdist(data_otu_filt_rar, method = 'bray'))
dist_bc[1:5, 1:5]

pcoa_bc = ordinate(data_phylo_filt_rar, "PCoA", "bray")
plot_ordination(data_phylo_filt_rar, pcoa_bc, color = "Location") +
  geom_point(size=10)

perm <- adonis2(data_otu_filt_rar~Location, data = meta_data, permutations = 9999, method = 'bray')

###PCoA Plot, calculate median

distmat <- phyloseq::distance(physeq, method = "bray")

pcoa_res <- cmdscale(distmat, eig = T, k = 2)

df <- as.data.frame(pcoa_res$points)
names(df) <- c("PC1", "PC2")

df$Sample_ID <- rownames(df)
df$Location <- physeq@sam_data$Location[match(df$Sample_ID, rownames(physeq@sam_data))]

med_pc1 <- median(df$PC1)
med_pc2 <- median(df$PC2)

p4 <- ggplot(df, aes(x = PC1, y = PC2)) +
  geom_point(aes(color = Location), size = 10, alpha = 0.5) +
  geom_text(aes(label = Sample_ID), size = 4, nudge_y = 0.05) +
  geom_hline(yintercept = med_pc2, linetype = "dashed", color = "red") +
  geom_vline(xintercept = med_pc1, linetype = "dashed", color = "red") +
  scale_color_manual(values = c('violet', 'cyan')) +
  theme(
    axis.text.x = element_text(size = 25),
    axis.text.y = element_text(size = 25),
    axis.title.x = element_text(size = 35),
    axis.title.y = element_text(size = 35),
    legend.title = element_text(size = 35),
    legend.text = element_text(size = 25),
    legend.key.size = unit(1.5, 'cm')
  )
p4

###PCoA Simple

pc <- pcoa_phyloseq(physeq, c("Location"), method = 'bray',
                    circle = 0.75, colors = c('violet', 'cyan')) + 
  geom_point(size=10, color = "black", alpha=0.5) +
  geom_text(aes(label=Sample_ID), size=5, nudge_y=0.06) +
  geom_hline(yintercept = med_pc2, linetype = "dashed", color = "red",
             linewidth=1) +
  geom_vline(xintercept = med_pc1, linetype = "dashed", color = "red",
             linewidth = 1) +
  theme(
    panel.border = element_rect(color = 'black', fill = NA, size = 2) 
  ) +
  theme(
    axis.text.x = element_text(size=20),
    axis.ticks.x = element_line(size=1),
    axis.text.y = element_text(size=20),
    axis.ticks.y = element_line(size=1),
    axis.title.x = element_text(size=25),
    axis.title.y = element_text(size=25),
    legend.title = element_text(size=25),
    legend.text = element_text(size=20),
    legend.key.size = (unit(1.5, 'cm'))
  )
pc
