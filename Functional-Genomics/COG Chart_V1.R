library(dplyr)
library(readr)
library(stringr)
library(ggplot2)
library(purrr)
library(forcats)
library(ggbreak)

genes_df <- read_csv("GSL17_113_Soil.csv")

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
  'S' = 'Function unknown',
  '_' = 'Hypothetical',
  '-' = 'Hypothetical'
)

broader_categories <- list(
  'Information Storage and Processing' = c('J', 'A', 'K', 'L', 'B'),
  'Cellular Processes and Signaling' = c('D', 'Y', 'V', 'T', 'M', 'N', 'Z', 'W', 'U', 'O'),
  'Metabolism' = c('C', 'G', 'E', 'F', 'H', 'I', 'P', 'Q'),
  'Poorly Characterized' = c('R', 'S'),
  'Unknown' = c('Hypothetical', 'Unknown')
)

##Function to map individual COG codes to their full names
get_cog_full_name <- function(cog_code) {
  unlist(str_split(cog_code, "")) %>%
    map_chr(~cog_definitions[.]) %>%
    paste(collapse = ' | ')
}

##Function to determine broader category
get_broader_category <- function(cog_code) {
  for (category in names(broader_categories)) {
    if (any(str_detect(cog_code, broader_categories[[category]]))) {
      return(category)
    }
  }
  return('Unknown')
}

##Apply function to create new columns
genes_df <- genes_df %>%
  mutate(
    COG_Full_Name = map_chr(COG_category, get_cog_full_name),
    Broad_Category = map_chr(COG_category, get_broader_category)
  )

# Count the occurrences of each combination of broader category and full name
category_counts <- genes_df %>%
  count(Broad_Category, COG_Full_Name) %>%
  arrange(Broad_Category, desc(n))

# Create ordered factors for consistent coloring and arranging
category_counts$Broad_Category <- factor(category_counts$Broad_Category, 
                                         levels = names(broader_categories))

p <- ggplot(category_counts, aes(x = fct_inorder(COG_Full_Name), y = n, fill = Broad_Category)) +
  geom_col(show.legend = FALSE) +
  labs(x = "COG Full Name", y = "Count") +
  scale_fill_manual(values = c("aquamarine3", "darkmagenta", "turquoise1", "coral2", "seagreen1")) +
  theme_minimal() +
  #scale_y_break(c(250,350), expand=c(0,0)) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1,
                               size = 10),
    axis.text.y = element_text(size=15),
    axis.title.y = element_text(size=20),
    axis.title.x = element_text(size=20)) 
p
