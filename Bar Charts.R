setwd()

library(phyloseq)
library(plyr)

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

physeq1 = merge_phyloseq(physeq, meta_data)
physeq1

library(ggplot2)
library(ggbreak)
theme_set(theme_bw())

###Actinomycetota
a = subset_taxa(physeq1, Phylum=="Actinobacteria")
p <- plot_bar(a, "Phylum", fill="Genus", facet_grid = ~Location) +
  labs(x="Phylum", y="Abundance (OTU)") +
  theme(
    axis.title.x = element_text(size=25),
    axis.title.y = element_text(size=25),
    axis.text.y = element_text(size=12),
    axis.text.x = element_text(size=12))
p

###Pseudomonadota
pro = subset_taxa(physeq1, Family=="Enterobacteriaceae")
pro1 <- plot_bar(pro, "Family", fill="Genus", facet_grid = ~Location) +
  theme(legend.position='right') +
  labs(x="Family", y="Abundance (OTU)") +
  theme(
    axis.title.x = element_text(size=25),
    axis.title.y = element_text(size=25),
    axis.text.y = element_text(size=12),
    axis.text.x = element_text(size=12, angle = 0, hjust = 0.5))
pro1

###Bar charts
phy = subset_taxa(physeq1, )
p = plot_bar(physeq1, "Phylum", fill="Phylum")
p + geom_bar(aes(color=Phylum, fill=Phylum), stat="identity") +
  theme(legend.position = 'none') +
  scale_y_break(c(6.5e+04, 1.5e+05)) + scale_y_break(c(1.7e+05, 3.5e+05)) +
  theme(
    axis.title.x = element_text(size=25),
    axis.title.y = element_text(size=25, angle=90),
    axis.text.x = element_text(size=15),
    axis.text.y = element_text(size=12)) +
  labs(y = "Abundance (ASV)") +
  coord_flip()

phydf <- psmelt(physeq1)

custom_col42 = c("#781156","#A51876","#D21E96","#E43FAD",
                 "#EA6CC0","#F098D3","#114578","#185EA5","#1E78D2",
                 "#3F91E4","#6CABEA","#98C4F0","#117878","#18A5A5",
                 "#3FE4E4","#6CEAEA","#98F0F0", "#117845","#18A55E",
                 "#1ED278","#3FE491","#6CEAAB","#98F0C4","#787811",
                 "#A5A518","#D2D21E","#E4E43F","#EAEA6C","#F0F098",
                 "#F7F7C5","#784511","#A55E18","#D2781E","#E4913F",
                 "#EAAB6C","#F0C498","#781122","#A5182F","#D21E2C",
                 "#E43F5B","#EA6C81","#F098A7", "#771155", 
                 "#AA4488", "#EA6CC0", "#CC99BB", "#114477",
                 "#4477AA","#1E78D2", "#77AADD", "#117777", 
                 "#44AAAA", "#3FE4E4", "#77CCCC", "#117744","#44AA77", 
                 "#1ED278", "#88CCAA", "#771122", "#AA4455", "#D21E2C",
                 "#DD7788","#777711", "#AAAA44", "#D2D21E", "#DDDD77",
                 "#774411", "#AA7744", "#D2781E", "#DDAA77", "#771155",
                 "#AA4488", "#CC99BB", "#114477", "#4477AA", "#77AADD", 
                 "#117777", "#44AAAA", "#77CCCC", "#117744", "#44AA77", 
                 "#88CCAA", "#777711", "#AAAA44", "#DDDD77", "#774411", 
                 "#AA7744", "#DDAA77", "#771122", "#AA4455", "#DD7788")

ggplot(phydf[order(phydf$Abundance),], aes(x=reorder(Phylum, Abundance), y=Abundance, fill=Phylum)) +
  geom_bar(stat="identity") +
  scale_fill_manual(values=custom_col42,
                    na.value="black") +
  theme(legend.position='none') +
  scale_y_break(c(6.5e+04, 1.5e+05)) + scale_y_break(c(1.7e+05, 3.5e+05)) +
  theme(
    axis.title.x = element_text(size=25),
    axis.title.y = element_text(size=25, angle=90),
    axis.text.x = element_text(size=15, angle=-90,
                               hjust=0.5),
    axis.text.y = element_text(size=12)) +
  labs(x = "Phylum", y = "Abundance (ASV)") +
  coord_flip()

ggplot(phydf[order(phydf$Abundance),], aes(x=reorder(Class, Abundance), y=Abundance, fill=Class)) +
  geom_bar(stat="identity") +
  scale_fill_manual(values=custom_col42,
                    na.value="black") +
  scale_y_break(c(4e+04, 5.5e+04)) +
  scale_y_break(c(62000,127500)) + scale_y_break(c(1.4e+05, 2.15e+05)) +
  scale_y_break(c(2.2e+05, 2.4e+05)) +
  theme(legend.position='none') +
  theme(
    axis.title.x = element_text(size=25),
    axis.title.y = element_text(size=25, angle=90),
    axis.text.x = element_text(size=15, angle=-90),
    axis.text.y = element_text(size=8)) +
  labs(x = "Class", y = "Abundance (ASV)") +
  coord_flip()
