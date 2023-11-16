setwd()

library(phyloseq)
library(plyr)
library(ggplot2)

###Import datasets
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

###Create OTU objects for physeq
OTU = otu_table(as.matrix(otumat), taxa_are_rows=T)
TAX = tax_table(as.matrix(taxmat))
OTU
TAX

###Create physeq object
physeq = phyloseq(OTU, TAX)
physeq

#merge meta data with physeq object, phyloseq method
physeq1 = merge_phyloseq(physeq, meta_data)
physeq1

head(sample_data(physeq1)$Location, 8)

library(DESeq2)

loc = phyloseq_to_deseq2(physeq1, ~ Location)
loc = DESeq(loc, test="Wald", fitType='parametric')

res = results(loc, cooksCutoff=F)
alpha = 0.01
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(physeq1)[rownames(sigtab), ], "matrix"))
head(sigtab)
#dimensions of matrix/array/dataframe
dim(sigtab)
# write.csv(sigtab, "sigtab.csv")
# write.csv(res, "res.csv")
summary(res)

ytitle <- expression(paste("Log"[2]*" Fold Change"))

# phylum order
x = tapply(sigtab$log2FoldChange, sigtab$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab$Phylum = factor(as.character(sigtab$Phylum), levels=names(x))
# Genus order
x = tapply(sigtab$log2FoldChange, sigtab$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtab$Genus = factor(as.character(sigtab$Genus), levels=names(x))

p <- ggplot(sigtab, aes(x=Genus, y=log2FoldChange, color=Phylum)) +
  theme_minimal() +
  geom_point(size=8) +
  labs(y=ytitle, x="Genus") +
  theme(
    axis.text.x = element_text(angle=-90, hjust=0, vjust=0.5),
    axis.title.x = element_text(size=20),
    axis.title.y = element_text(size=20),
    legend.position='bottom')
p + ggplot2::scale_color_manual(values=alpha(c(
  "gold2", "antiquewhite3", "cornsilk4", "darkgoldenrod", 
  "black", "dodgerblue",
  "chartreuse", "chartreuse3", "chartreuse4", 
  "brown", "firebrick1", "navy",
  "darkolivegreen1", "coral2", "coral4", "cyan", "cyan3", "cyan4",
  "darkorchid1", "darkorchid3", "darkorchid4", 
  "mediumspringgreen", "magenta"), 0.8))
