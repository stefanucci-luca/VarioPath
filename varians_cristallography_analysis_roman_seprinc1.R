
library(tidyverse)

set.seed(1)

all_serpin_in_ukb = read.delim("Desktop/VarioPath/serpinC1_cristallography/all_variants_in_UKB_in_serpinc_VEP_annotation.txt") %>% 
  filter(Protein_position != "-") %>% 
  select(Consequence, IMPACT, Protein_position, Amino_acids, gnomAD_AF) %>% 
  filter(Consequence != 'stop_gained' & Consequence != 'synonymous_variant')
all_serpin_in_ukb$source = "UKB_not_reported"

mdt_variants = xlsx::read.xlsx("Desktop/VarioPath/MDT_variants/V20210211/MDT_for_thrombosis.xls", sheetIndex = 1) %>% 
  filter(GENE == 'SERPINC1') %>% 
  filter(Protein_position != "-") %>% 
  select(Consequence, IMPACT, Protein_position, Amino_acids, gnomAD_AF) %>% 
  filter(Consequence != 'stop_gained')
mdt_variants$source = "MDT_spreadsheet"

variants_from_karyn = data.table::fread("/Volumes/GoogleDrive/My Drive/Phd/VarioPath/variants/VarioPath_VEP_annotation_only_variopath_transcripts_variants_20200615.tab") %>% 
  filter(SYMBOL == 'SERPINC1') %>% 
  filter(Protein_position != "-") %>% 
  select(Consequence, IMPACT, Protein_position, Amino_acids, gnomAD_AF) %>% 
  filter(Consequence != 'stop_gained') # %>% 
  filter(IMPACT == "HIGH")
variants_from_karyn$source = "karyn_list"

variants_from_karyn_pathogenicity = data.table::fread("/Volumes/GoogleDrive/My Drive/Phd/VarioPath/variants/VarioPath_variants_20200615.tsv") %>% 
  filter(QUAL == 'DM|Pathogenic' ) %>% 
  filter(FILTER == "ClinVar; curatedPub") 
variants_from_karyn_pathogenicity_serpin = grep( grep(pattern = "PROT=", unlist(str_split(variants_from_karyn_pathogenicity$INFO, pattern = "\\||\\,")),value = T), pattern = 'NP_000479.1|NP_001351981.1|NP_001373231.1', value = T)

var_kar_not_ukb = data_frame(
  "Protein_position" = c(23, 425, 424, 425, 381),
  "Amino_acids" = c("L/P", "R/H", "G/D", "R/C", "S/P"),
  "source" = "karyn_list"
           )

random_25_ukb = all_serpin_in_ukb[sample(x = 1:dim(all_serpin_in_ukb)[1], replace = F, size = 25),] %>% 
  select(Protein_position, Amino_acids, source )

rbind(random_25_ukb, var_kar_not_ukb)

write_csv(random_25_ukb, file = 'Desktop/VarioPath/serpinC1_cristallography/random_25_ukb.csv', quote_escape = 'none' )
write_csv(var_kar_not_ukb, file = 'Desktop/VarioPath/serpinC1_cristallography/var_kar_not_ukb.csv', quote_escape = 'none' )
write_csv(mdt_variants[,c("Protein_position","Amino_acids","source")], file = 'Desktop/VarioPath/serpinC1_cristallography/mdt_variants.csv', quote_escape = 'none' )

