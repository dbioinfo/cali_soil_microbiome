library(tidyverse)
library(dada2)
library(phyloseq)
library(Biostrings)

#Tailoring dada2

#import metadata
setwd('/scratch/')
meta <- read_csv('data/tailoring/metadata.csv')


##### Start with 16S 
fnFs <- paste0('raw_data/Pietrasiak2_Project_006/',meta$V4_R1) 
fnRs <- paste0('raw_data/Pietrasiak2_Project_006/',meta$V4_R2)

#plot qc
for (i in 1:5){
    png(paste0('figs/tailoring/dada2_fwd_qc',i,'.png'))
    print(plotQualityProfile(fnFs[(i*10):(i*10 +3)]))
    dev.off()
    png(paste0('figs/tailoring/dada2_rev_qc',i,'.png'))
    print(plotQualityProfile(fnRs[(i*10):(i*10+3)]))
    dev.off()
}


#set params according to qc
trunc_f <- 275
trunc_r <- 250
maxEE_f <- 2
maxEE_r <- 2

##Filter by params
filtFs <- file.path("raw_data/tailoring/filt_fastq/", meta$V4_R1)
filtRs <- file.path("raw_data/tailoring/filt_fastq/", meta$V4_R2)
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen = c(trunc_f, trunc_r),
                     maxN=0, maxEE=c(maxEE_f,maxEE_r), rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE)
write.csv(out, 'data/tailoring/16Sdada2/dada2_reads_filt.csv', row.names = T)

#filter samples for low reads after trimming (<25,000 reads.out)
lsamp <- data.frame(out) %>% filter(reads.out<25000) %>% rownames()
ctrls <- c("C6_1-V4_S164_R1_001.fastq.gz", "C6_2-V4_S166_R1_001.fastq.gz", "C6_3-V4_S168_R1_001.fastq.gz", 
           "C6_4-V4_S170_R1_001.fastq.gz", "C6_5-V4_S172_R1_001.fastq.gz", "C6_6-V4_S174_R1_001.fastq.gz",
           "ZS1PRE_1-V4_S998_R1_001.fastq.gz","ZS1PRE_2-V4_S1000_R1_001.fastq.gz",
           "ZS2POST_1-V4_S1002_R1_001.fastq.gz","ZS2POST_2-V4_S1004_R1_001.fastq.gz")
for (i in lsamp){
   if (!(i %in% ctrls)){
     if (file.exists(paste0("raw_data/tailoring/filt_fastq/",i))){
       file.remove(paste0("raw_data/tailoring/filt_fastq/",i))
     }
   }
}

#prepare for error rate learning
filt_fs <- c()
filt_rs <- c()
for (i in 1:nrow(meta)){
  if ((file.exists(paste0("raw_data/tailoring/filt_fastq/",meta[i,'V4_R1']))) & (file.exists(paste0("raw_data/tailoring/filt_fastq/",meta[i,'V4_R2'])))){
    filt_fs <- c(filt_fs, paste0("raw_data/tailoring/filt_fastq/",meta[i,'V4_R1']))
    filt_rs <- c(filt_rs, paste0("raw_data/tailoring/filt_fastq/",meta[i,'V4_R2']))
    }
}

##Main dada2 computation
#      #learn error rates 
err_f <- learnErrors(filt_fs, multithread = T)
err_r <- learnErrors(filt_rs, multithread = T)

#run dada2
dadaFs <- dada(filt_fs, err_f, multithread = T)
dadaRs <- dada(filt_rs, err_r, multithread = T)

intout <- list(err_f=err_f, err_r=err_r, dadaFs=dadaFs, dadaRs=dadaRs)
saveRDS(intout, 'data/tailoring/16Sdada2/dada2_errs.rds')

#intout <- readRDS('data/tailoring/16Sdada2/dada2_errs.rds')
#list2env(intout, .GlobalEnv)

##Merge fwd-rev reads
mergers <- mergePairs(dadaFs, filt_fs, dadaRs, filt_rs, verbose=TRUE)
saveRDS(mergers,'data/tailoring/16Sdada2/merged_reads_dada2.rds')

##Chimera removal
seqtab <- makeSequenceTable(mergers)
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
1 - sum(seqtab.nochim)/sum(seqtab) #what percent of reads were chimera/bimera
tax <- assignTaxonomy(seqtab.nochim, "raw_data/silva_ref/silva_nr99_v138.2_toGenus_trainset.fa.gz", multithread=TRUE)
#tax <- addSpecies(tax, "raw_data/silva_ref/silva_v138.2_assignSpecies.fa.gz")

#Prep phyloseq obj
rownames(seqtab.nochim) <- gsub('-V4_S.*_R1_001.fastq.gz','',rownames(seqtab.nochim)) 
rownames(seqtab.nochim)
imeta <- meta %>% as.data.frame() %>% column_to_rownames('SampleID')
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows = F), sample_data(imeta), tax_table(tax))
#fix ASV names, exact seqs are annoying
dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))

saveRDS(ps, 'data/tailoring/16Sdada2/phyloseq_dada2.rds')

