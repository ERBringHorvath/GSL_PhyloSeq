library(phyloseq)
library(ggplot2)
library(plyr)
library(vegan)

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
physeq

metadata <- as(sample_data(physeq), "data.frame")

#Bray is good for abundance data
perm <- adonis2(distance(physeq, method="bray", permutations = 999) ~ Location,
               data = metadata)
perm
write.csv(perm, file = "PERMANOVA_bray.csv")

#Euclidean is good for continuous data
perm.s <- adonis2(distance(physeq, method = "euclidean", permutations = 999) ~ Location,
                  data = metadata)
perm.s

###If we wanted to manually select locations to compare
physeq.subs <- subset_samples(physeq, Location %in% c("BRB", "Marina"))

