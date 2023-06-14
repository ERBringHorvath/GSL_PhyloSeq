library(phyloseq)
library(ggplot2)
library(plyr)

theme_set(theme_bw())

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

###Plot
nmds.horn <- ordinate(physeq, method='NMDS', distance='horn')
p <- plot_ordination(physeq, nmds.horn, type='samples', color="Location") +
  theme_bw() +
  coord_fixed(ratio=1) +
  geom_point(size=10) +
  geom_text(aes(label=Sample_ID), size=4, nudge_y=0.05) +
  theme(
    axis.text.x = element_text(size=15),
    axis.text.y = element_text(size=15),
    axis.title.x = element_text(size=15),
    axis.title.y = element_text(size=15),
    legend.title = element_text(size=15),
    legend.text = element_text(size=15),
    legend.key.size = (unit(1.5, 'cm'))
  )
p
