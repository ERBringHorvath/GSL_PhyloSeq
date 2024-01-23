library(dplyr)
library(readr)
library(stringr)
library(ggplot2)
library(purrr)
library(forcats)
library(ggbreak)

genes_df <- read_csv('')

genes_df <- genes_df %>%
  mutate(
    COG_First_Letter = substr(COG_category, 1, 1)
  )

#Define COG categories and broader category groups
cog_definitions <- c(
  'J' = 'Transl, rib struc and biogen',
  'A' = 'RNA proc and mod',
  'K' = 'Transc',
  'L' = 'Repl, recomb and repair',
  'B' = 'Chroma struc and dyna',
  'D' = 'Cell cycle ctrl, cell div, chromo part',
  'Y' = 'Nuc struc',
  'V' = 'Defense mech',
  'T' = 'Sig transd mech',
  'M' = 'Cell wall/memb/env biogen',
  'N' = 'Cell motil',
  'Z' = 'Cytoskeleton',
  'W' = 'Extracell struc',
  'U' = 'Intracell traf, sec, and vesic trans',
  'O' = 'Posttrans mod, prot t/o, chaperones', #t/o = turnover
  'C' = 'Energy prod and conv',
  'G' = 'Carbo trans and metab',
  'E' = 'AA transp and metab',
  'F' = 'Nucl transp and metab',
  'H' = 'Coenzyme transp and metab',
  'I' = 'Lipid transp and metab',
  'P' = 'Inorg ion transp and metab',
  'Q' = 'Sec metabol biosyn, transp and catab',
  'R' = 'General funct pred',
  'S' = 'Function unknown'
)

broader_categories <- list(
  'Information Storage and Processing' = c('J', 'A', 'K', 'L', 'B'),
  'Cellular Processes and Signaling' = c('D', 'Y', 'V', 'T', 'M', 'N', 'Z', 'W', 'U', 'O'),
  'Metabolism' = c('C', 'G', 'E', 'F', 'H', 'I', 'P', 'Q'),
  'Poorly Characterized' = c('R', 'S')
)

get_broader_category <- function(cog_code) {
  for (category in names(broader_categories)) {
    if (cog_code %in% broader_categories[[category]]) {
      return(category)
    }
  }
  return('Poorly Characterized')
}

# Apply function to create new columns for broader categories
genes_df <- genes_df %>%
  mutate(
    Broad_Category = map_chr(COG_category, get_broader_category)
  )

# For Figure A: Count the occurrences of each broader category per org
broad_counts <- genes_df %>%
  count(org, Broad_Category)

# Plot Figure A
p1 <- ggplot(broad_counts, aes(x = org, y = n, fill = Broad_Category)) +
  geom_bar(stat = "identity", position = "stack",
           color="lightgrey", linewidth=0.25) +
  scale_fill_manual(values = c('Information Storage and Processing' = 'aquamarine3',
                               'Cellular Processes and Signaling' = 'darkmagenta',
                               'Metabolism' = 'turquoise1',
                               'Poorly Characterized' = 'coral2',
                               'Unknown' = 'coral2')) +
  theme_minimal() +
  labs(y="Gene Count") +
  theme(
    axis.text.x = element_text(size = 20),
    axis.text.y = element_text(size = 20),
    axis.title.y = element_text(size = 35),
    legend.key.size = unit(1, 'cm'),
    legend.text = element_text(size = 20),
    legend.title = element_text(size = 30))
p1

genes_df <- genes_df %>%
  mutate(
    COG_First_Letter = substr(COG_category, 1, 1)
  )

# Count the occurrences of each COG subcategory per org and include the first letter for coloring
sub_counts <- genes_df %>%
  count(org, COG_category, name = "count") %>%
  mutate(COG_First_Letter = substr(COG_category, 1, 1))

# Define your custom color palette
my_colors <- c(
  'J' = '#e6194B', 'A' = '#3cb44b', 'K' = '#ffe119', 'L' = '#0082c8', 
  'B' = '#f58231', 'M' = '#911eb4', 'Y' = '#46f0f0', 'V' = '#f032e6', 
  'T' = '#d2f53c', 'D' = '#fabebe', 'N' = '#008080', 'Z' = '#e6beff', 
  'W' = '#aa6e28', 'U' = '#fffac8', 'O' = '#800000', 'C' = '#aaffc3', 
  'G' = '#808000', 'E' = '#ffd8b1', 'F' = '#000080', 'H' = '#808080', 
  'S' = 'turquoise', 'P' = '#000000', 'Q' = '#0A5F38', 'R' = '#4363d8', 
  'I' = '#fabebe'
)

# Apply these colors to your plot
p2 <- ggplot(sub_counts, aes(x = org, y = count, fill = COG_First_Letter)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = my_colors) +
  theme_minimal() +
  labs(y="Gene Count") +
  theme(
    axis.text.x = element_text(size = 20),
    axis.text.y = element_text(size = 20),
    axis.title.y = element_text(size = 35),
    legend.key.size = unit(1, 'cm'),
    legend.text = element_text(size = 20),
    legend.title = element_text(size = 30)
  )

p2
