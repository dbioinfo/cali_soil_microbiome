library(tidyverse)
library(vegan)
library(phyloseq)

setwd("/scratch/")
sitepal <- c("#A37BA7","#80D0B2","#EEA84D","#9ED080","#F8A582","#F8E082","#E694B5","#F88299")

#16S 
ps <- readRDS("data/tailoring/16Sdada2/phyloseq_dada2.filtered.rds")
meta <- data.frame(sample_data(ps)) %>% rownames_to_column("SampleID")

mat <- data.frame(otu_table(ps))
X <- vegdist(mat, method='bray')
saveRDS(X, 'data/tailoring/bray.dist.16S.rds')
out <- metaMDS(X,
               k=3,
               wascores = T,
               weakties=T,
               try=50,
               trymax=300,
               parallel=80,
               maxit=300)
saveRDS(out, 'data/tailoring/16S_nmds.rds')

gdat <- data.frame(out$points)
gdat$goodness <- goodness(out)
gdat$SampleID <- rownames(gdat)
gdat <- left_join(gdat, meta, by='SampleID')

stressplot(out)

centroids <- gdat %>% group_by(Site) %>% summarize(MDS1=mean(MDS1),MDS2=mean(MDS2))
gg<- ggplot(gdat)+
  geom_point(aes(x=MDS1,y=MDS2, color=Site))+
  geom_point(data=centroids, mapping=aes(x=MDS1, y=MDS2, color=Site), shape=4, size=5, stroke=2)+
  stat_ellipse(aes(x=MDS1, y=MDS2, group=Site, color=Site))+
  scale_color_manual(values=sitepal)+
  theme_bw()+
  ggtitle("16S Bray-Curtis NMDS")
gg
ggsave("figs/tailoring/beta_bc_nmds_16S.png",gg)


#ITS2
ps <- readRDS("data/tailoring/ITS2dada2/phyloseq_dada2.filtered.rds")
meta <- data.frame(sample_data(ps)) %>% rownames_to_column("SampleID")

mat <- data.frame(otu_table(ps))
X <- vegdist(mat, method='bray')
saveRDS(X, 'data/tailoring/bray.dist.ITS2.rds')
out <- metaMDS(X,
               k=3,
               wascores = T,
               weakties=T,
               try=50,
               trymax=100,
               parallel=80,
               maxit=300)
saveRDS(out, 'data/tailoring/ITS2_nmds.rds')

gdat <- data.frame(out$points)
gdat$goodness <- goodness(out)
gdat$SampleID <- rownames(gdat)
gdat <- left_join(gdat, meta, by='SampleID')

stressplot(out)

centroids <- gdat %>% group_by(Site) %>% summarize(MDS1=mean(MDS1),MDS2=mean(MDS2))
gg<- ggplot(gdat)+
  geom_point(aes(x=MDS1,y=MDS2, color=Site))+
  geom_point(data=centroids, mapping=aes(x=MDS1, y=MDS2, color=Site), shape=4, size=5, stroke=2)+
  stat_ellipse(aes(x=MDS1, y=MDS2, group=Site, color=Site))+
  scale_color_manual(values=sitepal)+
  theme_bw()+
  ggtitle("ITS2 Bray-Curtis NMDS")
gg
ggsave("figs/tailoring/beta_bc_nmds_ITS2.png",gg)

