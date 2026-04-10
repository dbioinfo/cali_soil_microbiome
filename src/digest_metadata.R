library(tidyverse)
library(readxl)

setwd('~/WorkForaging/Academia/Nicole/lindsay/')
xldat <- read_xlsx('raw_data/BLM-CA Tailoring Soil microbiome Workbook Edit 2026 04 06.xlsx', sheet=2, range="A5:M581")

colnames(xldat) <- c("Project","SampleID",
                     "Site","TransectNo","MicrositeType","MicrositeNo",
                     "Control",
                     "Volume","Concentration","Mass",
                     "PlateNo","Well","PlateID")
#write full data before filtering
write_csv(xldat, 'data/full_metadata.csv')

#subset for Stim Regen
tmp <- xldat %>% filter(Project=="Stimulating Regeneration 2024") %>% select(-c(Control, TransectNo)) %>% mutate(SampleID=gsub('^0','',SampleID))
fnames <- list.files('raw_data/Pietrasiak2_Project_006/', pattern = '*gz')
ttmp <- tibble(filename=fnames) %>%   
  mutate(mat = str_match(filename,   "^(.+?)-(ITS2|V4)_S(\\d+)_R(\\d+)_001\\.fastq\\.gz$")) %>%
  mutate(SampleID = mat[, 2],
    ftype = paste0(mat[, 3],'_R',mat[, 5])) %>%
  select(-mat)%>% 
  pivot_wider(id_cols = SampleID, names_from = ftype, values_from = filename)
ttmp <- left_join(tmp, ttmp, by="SampleID")
write_csv(ttmp, 'data/stim_regen/metadata.csv')

#subset for Tailoring
tmp <- xldat %>% filter(Project=="Tailoring Restoration 2025") %>% select(-Control)
fnames <- list.files('raw_data/Pietrasiak2_Project_006/', pattern = '*gz')
ttmp <- tibble(filename=fnames) %>%   
  mutate(mat = str_match(filename,   "^(.+?)-(ITS2|V4)_S(\\d+)_R(\\d+)_001\\.fastq\\.gz$")) %>%
  mutate(SampleID = mat[, 2],
         ftype = paste0(mat[, 3],'_R',mat[, 5])) %>%
  select(-mat)%>% 
  pivot_wider(id_cols = SampleID, names_from = ftype, values_from = filename)
ttmp <- left_join(tmp, ttmp, by="SampleID")
write_csv(ttmp, 'data/tailoring/metadata.csv')

