library(readr)
library(dplyr)
require(reshape)
library(stringr)
library(jsonlite)

# General settings
project_dir <- "~/Desktop/VarioPath/"
setwd(project_dir)

# Karyn variants pathogenyc and likely pathogenic
patho_likely_patho_df <- read_tsv("/Volumes/GoogleDrive/My Drive/Phd/VarioPath/variants/ALL_variants_P_LP_someVUS_original_karyn.txt",
                                  guess_max = 290000
)

# dim(patho_likely_patho_df)
# 299606      9

#--------------------------------------------------------------------------------------------------
### Add protine variant information on the manually curated variants
# The conversion to protein location and aminoacid substittution has been retrieved using the VariantValidator web tool.
c_to_p <- read_tsv("~/Desktop/listo_converted_curated_karyn_batch_job.txt", skip = 2, col_names = T,
                                  guess_max = 290000
)
# extract the info for the manually curated variants from Karyn's list
c_variants = as.data.frame(unlist(
                    lapply( 
                      apply(patho_likely_patho_df[which(patho_likely_patho_df$._2 == "curatedMix;"),], 1, function(x) 
                        unique(gsub(".*=","",grep("INPUT", x, value = TRUE)))),
                      '[', 1)
                      )
                   )

colnames(c_variants) <- "variants"

# remove semicolumn at the end of the variants
c_variants$variants <- str_remove_all(c_variants$variants, ";")

c_to_p_df = merge(c_variants, 
          c_to_p, 
          by.x = "variants",
          by.y="Input") %>% 
          select( "GRCh38_CHR", "GRCh38_POS", "GRCh38_REF", "GRCh38_ALT",  "variants", "HGVS_Predicted_Protein")


c_to_p_df$HGVS_Predicted_Protein <- paste0( "PROT=", c_to_p_df$HGVS_Predicted_Protein)


patho_likely_patho_df <- merge(patho_likely_patho_df, 
                               c_to_p_df, 
                               by.x= c("Chr", "Start", "REF", "ALT" ), 
                               by.y= c("GRCh38_CHR", "GRCh38_POS", "GRCh38_REF", "GRCh38_ALT"), 
                               all.x = T)

#dim(patho_likely_patho_df)
#299606     11

# There are about 70 variants missing after the merge step, try to include them. 
#--------------------------------------------------------------------------------------------------

# Try to parse INFO columns
patho_likely_patho_df = transform(patho_likely_patho_df, INFO = colsplit(INFO, split = ":|;", names = "INFO"))
patho_likely_patho_df = transform(patho_likely_patho_df, INFO_tag = colsplit(INFO_tag, split = ":|;", names = "INFO_tag"))

# extract information from the info and info_tag colum
info_to_extract = c("GENE","CLASS","CLNSIG","DNA" ,"PROT")

for (info_value in info_to_extract) {
  message("extracting info ", info_value )
  patho_likely_patho_df[,info_value] = unlist(
                                              lapply( 
                                                      apply(patho_likely_patho_df, 1, function(x) 
                                                        unique(gsub(".*=","",grep(info_value, x, value = TRUE)))),
                                                '[', 1)
                                            )
} 


keep_columns = c("Chr", 
                 "Start", 
                 "gene", 
                 "._2",
                 info_to_extract)

df_final <- patho_likely_patho_df[,which(colnames(patho_likely_patho_df) %in% keep_columns)]
colnames(df_final)[3] <- "SOURCE"

View(df_final[which(df_final$GENE == "VWF"),])

json_variant_list = toJSON(unname(split(df_final, 1:nrow(df_final))))



