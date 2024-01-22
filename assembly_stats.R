###Calculate average coverage depth
contigs_data <- read.csv("111_Cov.csv")

weighted_coverage <- sum(contigs_data$Coverage * contigs_data$Length)

total_genome_length <- sum(contigs_data$Length)

average_coverage_depth <- weighted_coverage / total_genome_length

average_coverage_depth


###Calculate N50
library(data.table)

data <- fread("111_Cov.csv")

##Ensure Length column is numeric
data[, Length := as.numeric(Length)]

##Order data by contig length in descending order
data <- data[order(-Length)]

##Calculate the total length of all contigs
total_length <- sum(data$Length)

##Calculate N50
running_total <- 0
for (length in data$Length) {
  running_total <- running_total + length
  if (running_total >= total_length / 2) {
    N50 <- length
    break
  }
}

N50

##Calculate L50
running_total2 <- 0
L50 <- 0
for (length in data$Length) {
  running_total2 <- running_total2 + length
  L50 <- L50 + 1
  if (running_total2 >= total_length / 2) {
    break
  }
}

L50
