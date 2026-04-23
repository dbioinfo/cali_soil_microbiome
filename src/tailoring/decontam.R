
library(tidyverse)
library(vegan)
library(phyloseq)
library(decontam)

setwd('/scratch/')


#Start with 16S 
ps <- readRDS('data/tailoring/16Sdada2/phyloseq_dada2.rds')
sample_data(ps)$is.neg <- !(sample_data(ps)$Control=='Pos'|is.na(sample_data(ps)$Control) )

#run decontam
decdf <- isContaminant(ps, method='prevalence', neg='is.neg')
cindx <- which(decdf$contaminant)

#check abundances
gdat <- as.data.frame(otu_table(ps))[,cindx] 
gdat <- gdat %>% rownames_to_column("SampleID") %>% 
  pivot_longer(!SampleID, names_to = "ASV", values_to = "Abundance") %>% 
  left_join(., data.frame(sample_data(ps)) %>% rownames_to_column("SampleID") , by="SampleID") %>% 
  left_join(., data.frame(tax_table(ps)) %>% rownames_to_column("ASV") , by="ASV")

g1<-ggplot(gdat)+
  geom_tile(mapping=aes(x=SampleID, y=ASV, fill=Abundance))+
  scale_fill_gradient(low='black',high='purple')+
  theme_bw()+
  theme(axis.text=element_blank(), axis.ticks = element_blank())+
  xlab("Sample")+
  ggtitle("16S Contamination Matrix")
ggsave('figs/tailoring/contamination_matrix_16S.png',g1)

tdat <- gdat %>% group_by(Control, Phylum) %>% summarize(Abundance=sum(Abundance))
g2<-ggplot(tdat)+
  geom_bar(mapping=aes(x=Control, y=Abundance, fill=Phylum), stat='identity')+
  theme_bw()+
  scale_y_continuous(expand=c(0,0))+
  scale_fill_manual(values=c("#F8A582","#9ED080","#A37BA7","#88F1E1"))
ggsave('figs/tailoring/contamination_phylum_16S.png',g2)


#filter and save results
csamps <- data.frame(sample_data(ps)) %>% rownames_to_column("SampleID") %>% filter(!is.na(Control)) %>% pull(SampleID)
ps <- prune_taxa(!(taxa_names(ps) %in% unique(gdat$ASV)), ps)
ps <- prune_samples(!(rownames(data.frame(sample_data(ps))) %in% csamps),ps)
saveRDS(ps, 'data/tailoring/16Sdada2/phyloseq_dada2.filtered.rds')


#Then do ITS2 
ps <- readRDS('data/tailoring/ITS2dada2/phyloseq_dada2.rds')
sample_data(ps)$is.neg <- !(sample_data(ps)$Control=='Pos'|is.na(sample_data(ps)$Control) )

#run decontam
decdf <- isContaminant(ps, method='prevalence', neg='is.neg')
cindx <- which(decdf$contaminant)

#check abundances -- only 1 ASV detected 
cname <- colnames(as.data.frame(otu_table(ps)))[cindx] 
gdat <- as.data.frame(otu_table(ps)) %>% rownames_to_column("SampleID") %>% select(SampleID, cname) %>% 
  pivot_longer(!SampleID, names_to = "ASV", values_to = "Abundance") %>% 
  left_join(., data.frame(sample_data(ps)) %>% rownames_to_column("SampleID") , by="SampleID") %>% 
  left_join(., data.frame(tax_table(ps)) %>% rownames_to_column("ASV") , by="ASV")

g1<-ggplot(gdat)+
  geom_tile(mapping=aes(x=SampleID, y=ASV, fill=Abundance))+
  scale_fill_gradient(low='black',high='purple')+
  scale_y_discrete(expand = c(0,0))+
  theme_minimal()+
  theme(axis.text=element_blank(), axis.ticks = element_blank())+
  xlab("Sample")+
  ggtitle("ITS2 Contamination Matrix")
g1
ggsave('figs/tailoring/contamination_matrix_ITS2.png',g1)

tdat <- gdat %>% group_by(Control, ASV) %>% summarize(Abundance=sum(Abundance))
g2<-ggplot(tdat)+
  geom_bar(mapping=aes(x=Control, y=Abundance, fill=ASV), stat='identity')+
  theme_bw()+
  scale_y_continuous(expand=c(0,0))
  #scale_fill_manual(values=c("#A37BA7"))
g2
ggsave('figs/tailoring/contamination_asv_ITS2.png',g2)


#filter and save results
csamps <- data.frame(sample_data(ps)) %>% rownames_to_column("SampleID") %>% filter(!is.na(Control)) %>% pull(SampleID)
ps <- prune_taxa(!(taxa_names(ps) %in% unique(gdat$ASV)), ps)
ps <- prune_samples(!(rownames(data.frame(sample_data(ps))) %in% csamps),ps)
saveRDS(ps, 'data/tailoring/ITS2dada2/phyloseq_dada2.filtered.rds')
