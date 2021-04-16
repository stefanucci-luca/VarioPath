library(tidyverse)
library(dplyr)
library(ggplot2)
library(ggforce)

# Gene list for PID
googledrive::drive_download("https://docs.google.com/spreadsheets/d/1mDGFXJ9zHWlSezQaKNfbyhdfQo7hlTTB50qXw2fSPN0/edit#gid=29422767")
genes_df = readxl::read_xlsx('VarioPath_disease_gene_list_20201012v3.xlsx',
                           trim_ws = T)

file.remove('VarioPath_disease_gene_list_20201012v3.xlsx')

pid_df = genes_df %>% filter(stringr::str_detect(pattern = "Primary\ immune\ disorders" , string = genes_df$`Chosen disease domain(s)` ))

pid_gene_list = pid_df %>% select(c(`Approved symbol (HGNC)`,`MISC_MOIs (OMIM via ENSG)`))

# Variants
var_df = data.table::fread('/Volumes/GoogleDrive/My Drive/Phd/VarioPath/variants/VarioPath_VEP_annotation_only_variopath_transcripts_variants_20200615.tab')

# PID vartiants

var_df_pid = var_df[var_df$SYMBOL %in% pid_gene_list$`Approved symbol (HGNC)`]
# remove duplicated variants
var_df_pid = var_df_pid[!duplicated(var_df_pid$Uploaded_variation),]


# Variants
patho_df = data.table::fread('/Volumes/GoogleDrive/My Drive/Phd/VarioPath/variants/VarioPath_variants_20200615.tsv')
# convert annotation of variants to the one used in VarioPath
library("stringr")
df = as.data.frame(patho_df)
df$mapply_dist = mapply(function(x,y) which.min(x==y),strsplit(df$REF,""),
                        strsplit(df$ALT,""))
df <- dplyr::mutate_all(df, as.character)
df2 = df
for (variable in 1:dim(df)[1]) {
  if (df[variable,"mapply_dist"] == 2){
    if (str_length(df2[variable,"REF"]) > str_length(df2[variable,"ALT"]) ) {
      df2[variable,"POS"] <- as.numeric(df2[variable,"POS"]) + 1
      df2[variable,"REF"] <- paste0(na.omit(strsplit(df[variable,"REF"],"")[[1]][df[variable,"mapply_dist"]:length(strsplit(df[variable,"REF"],"")[[1]])]), collapse="")
      if (str_length(df2[variable,"ALT"]) == 1) {
        df2[variable,"ALT"] <- "-"
      } else {
        df2[variable,"ALT"] <- paste0(na.omit(strsplit(df[variable,"ALT"],"")[[1]][1:length(strsplit(df[variable,"ALT"],"")[[1]])]), collapse="")
      }
    } else {
      df2[variable,"POS"] <- as.numeric(df2[variable,"POS"]) + 1
      df2[variable,"ALT"] <- paste0(na.omit(strsplit(df[variable,"ALT"],"")[[1]][df[variable,"mapply_dist"]:length(strsplit(df[variable,"ALT"],"")[[1]])]), collapse="")
      if (str_length(df2[variable,"REF"]) == 1) {
        df2[variable,"REF"] <- "-"
      } else {
        df2[variable,"REF"] <- paste0(na.omit(strsplit(df[variable,"REF"],"")[[1]][1:length(strsplit(df[variable,"REF"],"")[[1]])]), collapse="")
      }
    }
  } else if (df[variable,"mapply_dist"] > 2){
    if (str_length(df2[variable,"REF"]) > str_length(df2[variable,"ALT"]) ) {
      df2[variable,"POS"] <- as.numeric(df2[variable,"POS"]) + 1
      df2[variable,"REF"] <- paste0(na.omit(strsplit(df[variable,"REF"],"")[[1]][2:length(strsplit(df[variable,"REF"],"")[[1]])]), collapse="")
      if (str_length(df2[variable,"ALT"]) == 1) {
        df2[variable,"ALT"] <- "-"
      } else {
        df2[variable,"ALT"] <- paste0(na.omit(strsplit(df[variable,"ALT"],"")[[1]][1:length(strsplit(df[variable,"ALT"],"")[[1]])]), collapse="")
      }
    } else {
      df2[variable,"POS"] <- as.numeric(df2[variable,"POS"]) + 1
      df2[variable,"ALT"] <- paste0(na.omit(strsplit(df[variable,"ALT"],"")[[1]][2:length(strsplit(df[variable,"ALT"],"")[[1]])]), collapse="")
      if (str_length(df2[variable,"REF"]) == 1) {
        df2[variable,"REF"] <- "-"
      } else {
        df2[variable,"REF"] <- paste0(na.omit(strsplit(df[variable,"REF"],"")[[1]][1:length(strsplit(df[variable,"REF"],"")[[1]])]), collapse="")
      }
    }
  }
}
# New IDs
patho_df$NEW_ID <- paste(df2$`#CHROM`, df2$POS, df2$REF, df2$ALT , sep="_")

#df with PID genes, pathogenicity and variants
var_df_pid_patho = merge(var_df_pid, patho_df, by.x = 'Uploaded_variation', by.y = 'NEW_ID')

#Attach MOI
var_df_pid_patho$moi = ''
for (i in 1:dim(var_df_pid_patho)[1]) {
  moi = pid_gene_list$`MISC_MOIs (OMIM via ENSG)`[which(pid_gene_list$`Approved symbol (HGNC)` %in% var_df_pid_patho$SYMBOL[i])[1]]
  var_df_pid_patho$moi[i] = moi
}

df_to_plot = var_df_pid_patho %>% count(SYMBOL, moi)

df_to_plot$moi_correct = as.character(
  apply(as.data.frame(df_to_plot[,2]), 
      1, 
      function(x){
        a = str_replace_all( unique(str_split( x , " \\| ")[[1]]), pattern = 'ND|0\\.0', replacement = "") 
        return( paste( a[a!= ''], collapse = ",") )
        }
      )
  )


ggplot(data = df_to_plot,
       aes( reorder(SYMBOL, n), n, fill = moi_correct ))+
  geom_col() + 
  theme_bw() + 
  coord_flip() +
  guides(fill=guide_legend(ncol=1)) +
  scale_y_continuous(position = "right") +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 5),
        strip.text.x = element_text(angle = 0), 
        legend.position = 'top') +
  scale_fill_manual(values = ggsci::pal_igv('default')(15))
ggsave('Desktop/test_PID.pdf', device = 'pdf', height = unit(40,'cm'))
 

# To extract the variants from UKB cohort:
# bcftools view -R <( sed 's/^/chr/' PID_regions_to_extract_from_UKB.bed) ../VCF_200K_filtered_for_variopath/karyns_variant_from_200KWES_biallelic_records_only_unrelated_variant_normalised_sorted.bcf > PID_karyns_variant_from_200KWES_biallelic_records_only_unrelated_variant_normalised_sorted.bcf
# then join with VEP annotated variants
# join 
#