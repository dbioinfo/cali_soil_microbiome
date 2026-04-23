library(tidyverse)
library(vegan)
library(phyloseq)

setwd("/scratch/")

#16S
ps <- readRDS("data/tailoring/16Sdada2/phyloseq_dada2.filtered.rds")
meta <- data.frame(sample_data(ps)) %>% rownames_to_column("SampleID")

#check out rarefaction curve
otu <- data.frame(otu_table(ps))

g1 <- rarecurve(otu, step=1000, tidy=T)
g2 <- g1 %>% group_by(Site) %>% summarize(Sample=max(Sample), Species=max(Species))
gg<-ggplot(g1)+
  geom_line(aes(x=Sample, y=Species, group=Site))+
  geom_label(data=g2, aes(x=Sample, y=Species, label=Site))+
  geom_vline(aes(xintercept=10000), color='salmon')+
  theme_bw()+
  xlab("Sample Size")+
  ggtitle("16S Rarefaction Curve")
ggsave("figs/tailoring/alpha_rarefaction_curve_16S.png", gg)

#rarefy 
irare <- rrarefy(otu, 10000)

#calculate alphas
rich <- specnumber(irare)
shan <- diversity(irare, index = 'shannon')
simp <- diversity(irare, index = 'simpson')
peven <- shan / log(rich)

gdat <- data.frame(SampleID=names(rich),
                   Richness=rich, 
                   Shannon=shan, 
                   Simpson=simp, 
                   Pielou_Even=peven) %>% 
      pivot_longer(!SampleID, names_to = "Metric", values_to = "Alpha") %>% 
      left_join(., meta, by='SampleID')

sitepal <- c("#A37BA7","#80D0B2","#EEA84D","#9ED080","#F8A582","#F8E082","#E694B5","#F88299")
gg<-ggplot(gdat)+
  geom_boxplot(aes(x=Site, y=Alpha, fill=Site), outliers = F)+
  geom_jitter(aes(x=Site, y=Alpha), width=0.1, height=0, size=.75, alpha=0.8)+
  scale_fill_manual(values=sitepal)+
  facet_grid(Metric~MicrositeType, scales = 'free' )+
  theme_bw()+
  guides(fill='none')+
  ggtitle("16S Alpha Diversity By Site")
ggsave("figs/tailoring/alpha_panel_site_16S.png",gg)

#ITS2
ps <- readRDS("data/tailoring/ITS2dada2/phyloseq_dada2.filtered.rds")
meta <- data.frame(sample_data(ps)) %>% rownames_to_column("SampleID")

#check out rarefaction curve
otu <- data.frame(otu_table(ps))

g1 <- rarecurve(otu, step=1000, tidy=T)
g2 <- g1 %>% group_by(Site) %>% summarize(Sample=max(Sample), Species=max(Species))
gg<-ggplot(g1)+
  geom_line(aes(x=Sample, y=Species, group=Site))+
  geom_label(data=g2, aes(x=Sample, y=Species, label=Site))+
  geom_vline(aes(xintercept=20000), color='salmon')+
  theme_bw()+
  xlab("Sample Size")+
  ggtitle("ITS2 Rarefaction Curve")
ggsave("figs/tailoring/alpha_rarefaction_curve_ITS2.png", gg)

#rarefy 
irare <- rrarefy(otu, 20000)

#calculate alphas
rich <- specnumber(irare)
shan <- diversity(irare, index = 'shannon')
simp <- diversity(irare, index = 'simpson')
peven <- shan / log(rich)

gdat <- data.frame(SampleID=names(rich),
                   Richness=rich, 
                   Shannon=shan, 
                   Simpson=simp, 
                   Pielou_Even=peven) %>% 
  pivot_longer(!SampleID, names_to = "Metric", values_to = "Alpha") %>% 
  left_join(., meta, by='SampleID')

sitepal <- c("#A37BA7","#80D0B2","#EEA84D","#9ED080","#F8A582","#F8E082","#E694B5","#F88299")
gg<-ggplot(gdat)+
  geom_boxplot(aes(x=Site, y=Alpha, fill=Site), outliers = F)+
  geom_jitter(aes(x=Site, y=Alpha), width=0.1, height=0, size=.75, alpha=0.8)+
  scale_fill_manual(values=sitepal)+
  facet_grid(Metric~MicrositeType, scales = 'free' )+
  theme_bw()+
  guides(fill='none')+
  ggtitle("ITS2 Alpha Diversity By Site")
ggsave("figs/tailoring/alpha_panel_site_ITS2.png",gg)


