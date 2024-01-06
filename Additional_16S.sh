screen -S GSS
srun --pty R
library('dada2')
library('ShortRead')
library('ggplot2')
library('grid')
library('gridExtra')

path <- getwd()
reads <- list.files(path)
fastqs <- sort(reads[grepl('.fastq+', reads)])
ff <- fastqs[grepl('forward', fastqs)]
rf <- fastqs[grepl('reverse', fastqs)]
samples <- sapply(strsplit(ff, '[.]'), '[', 1)
ff <- paste0(path, '/', ff)
rf <- paste0(path, '/', rf)

n <- sample(length(ff), 4)
forward_plots <- list()
reverse_plots <- list()
for (i in 1:length(n)) {
  sample_index <- n[i]
  fp <- plotQualityProfile(ff[sample_index])
  rp <- plotQualityProfile(rf[sample_index])
  fp <- fp + ggtitle(samples[sample_index])
  rp <- rp + ggtitle(samples[sample_index])
  forward_plots[[i]] <- fp
  reverse_plots[[i]] <- rp
}

png(filename="forward_plots.png")
grid.arrange(grobs=forward_plots)
dev.off()
png(filename="reverse_plots.png")
grid.arrange(grobs=reverse_plots)
dev.off()

######
run remaining commands with r script:

nano dada.r

#!/usr/bin/Rscript
library('dada2')
library('ShortRead')
library('ggplot2')
library('grid')
library('gridExtra')

path <- getwd()
reads <- list.files(path)
fastqs <- sort(reads[grepl('.fastq+', reads)])
ff <- fastqs[grepl('forward', fastqs)]
rf <- fastqs[grepl('reverse', fastqs)]
samples <- sapply(strsplit(ff, '[.]'), '[', 1)
ff <- paste0(path, '/', ff)
rf <- paste0(path, '/', rf)

ff.filt <- paste0(path, '/', samples, '.forward.trimmed.filtered.fastq.gz')
rf.filt <- paste0(path, '/', samples, '.reverse.trimmed.filtered.fastq.gz')
filtered <- filterAndTrim(ff, ff.filt, rf, rf.filt, maxN=0, maxEE=3, truncQ=2, compress=TRUE, verbose=TRUE, matchIDs=TRUE, rm.phix=TRUE)
ff.derep <- derepFastq(ff.filt, verbose=TRUE)
rf.derep <- derepFastq(rf.filt, verbose=TRUE)
names(ff.derep) <- samples
names(rf.derep) <- samples
n <- sample(length(ff), 5)
ff.err <- learnErrors(ff.derep[n], errorEstimationFunction=loessErrfun, randomize=FALSE, multithread=TRUE)
rf.err <- learnErrors(rf.derep[n], errorEstimationFunction=loessErrfun, randomize=FALSE, multithread=TRUE)
ff.dada <- dada(ff.derep, err=ff.err, pool=TRUE, multithread=TRUE)
rf.dada <- dada(rf.derep, err=rf.err, pool=TRUE, multithread=TRUE)
merged <- mergePairs(ff.dada, ff.derep, rf.dada, rf.derep, verbose=TRUE, maxMismatch=0, trimOverhang=FALSE, justConcatenate=FALSE)
seqtable <- makeSequenceTable(merged, orderBy='abundance')
seqtable.nochim <- removeBimeraDenovo(seqtable, method="consensus", verbose=TRUE)
save.image()

srun -p highmem --mem 200G Rscript dada.r &


1 Aug 2019
screen -S GSS
srun -p highmem --mem 100G --pty R
library('dada2')
library('ShortRead')
library('ggplot2')
library('grid')
library('gridExtra')

dim(seqtable)
[1]    8 3373

table(nchar(colnames(seqtable)))
250  251  252  253  254  255  256  257  261  264  292  293  357  365  377  395
   9    2  127 2850  331   30    8    1    3    1    1    5    1    1    1    1
 419
   1
  
getN <- function(x) sum(getUniques(x))
track <- cbind(sapply(ff.dada, getN), sapply(merged, getN), rowSums(seqtable.nochim))
colnames(track) <- c("denoised", "merged", "chimera-checked")
rownames(track) <- samples
track
