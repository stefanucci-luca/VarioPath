library(tidyverse)
# Clean dataframe
# Read in the data.frame
df = data.table::fread("/Users/luca/Desktop/VarioPath/MDT_master_spreadsheet_to_clean.tab")
# create a working copy
df_clean = df
dim(df_clean)
# Select relevant colums
# colnames(df)
# All columns
#   "NEW_ID","CHROM", "ORIG_ID","REF","ALT","EUR_ALT_FREQS", "EUR_OBS_CT","AFR_ALT_FREQS","AFR_OBS_CT","AMR_ALT_FREQS","AMR_OBS_CT",
#   "EAS_ALT_FREQS","EAS_OBS_CT","SAS_ALT_FREQS","SAS_OBS_CT","Location",  "Allele","Gene",  "Feature","Feature_type", "Consequence",  "cDNA_position",
#   "CDS_position", "Protein_position","Amino_acids",  "Codons","Existing_variation", "IMPACT","DISTANCE", "STRAND","FLAGS", "SYMBOL","SYMBOL_SOURCE",
#   "HGNC_ID","CANONICAL", "MANE",  "ENSP",  "SWISSPROT", "TREMBL","UNIPARC","REFSEQ_MATCH", "SOURCE","REFSEQ_OFFSET","GIVEN_REF",
#   "USED_REF",  "BAM_EDIT",  "HGVSc", "HGVSp", "HGVS_OFFSET", "gnomAD_AF", "gnomAD_AFR_AF","gnomAD_AMR_AF","gnomAD_ASJ_AF","gnomAD_EAS_AF","gnomAD_FIN_AF",
#   "gnomAD_NFE_AF","gnomAD_OTH_AF","gnomAD_SAS_AF","CLIN_SIG","SOMATIC","PHENO", "CADD_PHRED","CADD_RAW","CHROM", "POS", "REF","ALT","OLD_ID","UKB_AF","PARTECIPANTS" 
# Keep 
col_to_keep = c("CHROM", "POS", "REF","ALT","CHROM","REF","ALT","CADD_PHRED", "Consequence","Protein_position","Amino_acids","HGVSp",
                "Existing_variation", "IMPACT", "SYMBOL",
                "UKB_AF", "EUR_ALT_FREQS", "EUR_OBS_CT","AFR_ALT_FREQS","AFR_OBS_CT","AMR_ALT_FREQS","AMR_OBS_CT",
                "EAS_ALT_FREQS","EAS_OBS_CT","SAS_ALT_FREQS","SAS_OBS_CT", "gnomAD_AF", "gnomAD_AFR_AF","gnomAD_AMR_AF","gnomAD_ASJ_AF","gnomAD_EAS_AF","gnomAD_FIN_AF",
                "gnomAD_NFE_AF","gnomAD_OTH_AF","gnomAD_SAS_AF","PARTECIPANTS" )

#____________________________________________________________________
# Select what impact to retain in the inal spreadsheet
# unique(df$IMPACT)
# "MODERATE" "HIGH"     "MODIFIER" "LOW"  
# selecting everything but LOW i have 1334 variants (1390 in the original df)

#___________________________________________________________________
df_clean = df %>% select(col_to_keep) %>%filter(IMPACT != "LOW")


#__________________________________________________________________
# Add pLI info
# import transcript
transcripts = read_tsv("/Users/luca/Desktop/VarioPath/MDT_transcript", col_names = F)
gnomad_pli = read_tsv("/Users/luca/Desktop/VarioPath/gnomad_paper_PMID32461654_supplementary_dataset_11_full_constraint_metrics.tsv") %>% 
  filter(transcript %in% transcripts$X1) %>% select("gene","pLI") 
# merge tables
df_clean = merge(df_clean,gnomad_pli, by.x = "SYMBOL", by.y="gene")
# move pLI column to where it belongs!
df_clean = relocate(df_clean, pLI, .after = CADD_PHRED)

#_________________________________________________________________
# Add MOI
# MOI source df from Karyn's paper
moi_df = readxl::read_xlsx("/Users/luca/Desktop/VarioPath/isth2020-1_genelist_supptab-4.xlsx") %>% 
  select("Inheritance","Gene symbol (HGNC)")
# Add MOI info
df_clean = merge(df_clean,moi_df, by.x = "SYMBOL", by.y="Gene symbol (HGNC)", allow.cartesian=TRUE,all.x = T)
# move MOI column to where it belongs!
df_clean = relocate(df_clean, Inheritance, .after = pLI)

#_________________________________________________________________
# created groups with genes fro every MDT
blee_coag = c("F10","F11","F12","F13A1","F13B","F2",
		"F5","F7","F8","F9","FGA","FGB","FGG",
		"GGCX","KNG1","LMAN1","MCFD2","SERPINE1",
		"SERPINF2","VKORC1","VWF")
thrombosis = c("ADAMTS13","HRG","PIGA","PLG","PROC",
		"PROS1","SERPINC1","SERPIND1","THBD")
platelet = c("ABCC4","ABCG5","ABCG8","ACTB","ACTN1","ANKRD26",
		"ANO6","AP3B1","AP3D1","ARPC1B","BLOC1S3","BLOC1S6",
		"CDC42","CYCS","DIAPH1","DTNBP1","ETV6","FERMT3","FLI1",
		"FLNA","FYB1","GATA1","GFI1B","GNE","GP1BA","GP1BB","GP6",
		"GP9","HOXA11","HPS1","HPS3","HPS4","HPS5","HPS6","IKZF5",
		"ITGA2B","ITGB3","KDSR","LYST","MECOM","MPIG6B","MPL","MYH9",
		"NBEA","NBEAL2","P2RY12","PLA2G4A","PLAU","RASGRP2","RBM8A",
		"RNU4ATAC","RUNX1","SLFN14","SRC","STIM1","STXBP2","TBXA2R",
		"TBXAS1","THPO","TUBB1","VIPAS39","VPS33B","WAS")
hered_sfero = c("ANK1","EPB41","EPB42","SLC4A1","SPTA1","SPTB")
# define the MDT names
MDTs = c("bleeding_and_coagulation",
	 "thrombosis", "platelet",
	 "hereditary_sferocytosis")
#___________________________________________________________________________________
# Create a column for each MDT and populate it with 0
for (var in MDTs) {
  df_clean[,var] = 0
}
# Sustitute 0 with 1 if the gene belong to that MDT
for (li in 1:dim(df_clean)[1]) {
            if (df_clean[li,"SYMBOL"] %in% blee_coag) {
              df_clean[li,"bleeding_and_coagulation"] = 1 
            } 
            if(df_clean[li,"SYMBOL"] %in% thrombosis) {
              df_clean[li,"thrombosis"] = 1 
            }
            if(df_clean[li,"SYMBOL"] %in% platelet) {
              df_clean[li,"platelet"] = 1 
            }
            if(df_clean[li,"SYMBOL"] %in% hered_sfero) {
              df_clean[li,"hereditary_sferocytosis"] = 1 
            }
}
dir.create("/Users/luca/Desktop/VarioPath/MDT_variatns/")
# Create MDT spefic df
for (domain in MDTs) {
  col_pos = which(colnames(df_clean) %in% domain)
  tmp_df = df_clean[which(df_clean[,..col_pos] == 1),]
  mutate_all(tmp_df, as.factor)
  filename = paste("/Users/luca/Desktop/VarioPath/MDT_variatns/","MDT_for_", domain, ".xls", sep = "")
  WriteXLS::WriteXLS(tmp_df, filename, FreezeRow = 1, col.names = T)
  message("created file: ",filename)
  message("it has ", dim(tmp_df)[1]," variatns")
}

