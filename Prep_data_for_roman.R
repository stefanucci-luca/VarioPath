library(dplyr)
library(stringr)
# set working direcotry
setwd("~/Desktop/VarioPath/VarioPath_scripts/")
# Functions to use
# check if ref and alt are resulting in SNV. i.e. every variants composed of nucleotides > 1 and where alt != ref
check.if.SNV = function(df, col_ref, col_alt){
               require(stringr)
               apply(df, 
                     1, 
                     function(x)
                     str_length(x[eval(col_ref)]) == str_length(x[eval(col_alt)]) & str_length(x[eval(col_ref)]) == 1 )
}
# extract gene name from a column (e.g. info column)
extract.feature = function(df, col_info, feature){
  require(reshape)
  require(stringr)
  # extract information from INFO column
  unlist( 
  lapply(
  apply(df, 
        1, 
        function(x)
          unique(gsub(".*=",
                      "",
                      grep(eval(feature),
                           str_split(x[eval(col_info)], pattern = ":|;|,|\\|", simplify = T),
                                    value = TRUE)
                      )
                 )
          ),
  '[',
  1)
  )
}
# Assign variant bin. according to these guidelines:
          # 1. Not in UKB = "bin_4"
          # 2. In UKB ultra rare (<=1:10,000) = "bin_1"
          # 3. In UKB rare (1:1,000-1:10,000) => discussed in VarioPath MDTs = "bin_2"
          # 4. In UKB too frequent (=> 1:1,000) = = "bin_3"
assign.bin = function(list_AF){
  unlist(
    lapply(list_AF, 
           function(x) 
             if ( is.na(x) ) {
               return('bin_4')  
             } else if (x <= 1/10000) {
               return('bin_1')
             } else if (x > 1/10000 & x < 1/1000) {
               return('bin_2')
             }  else if (x >= 1/1000) {
               return('bin_3')  
             } 
    )
  )
}
# Import the master list with VEP annotated 300K variants from gdrive
googledrive::drive_download("https://drive.google.com/file/d/1WQgNdV-kD7Dq5Q9ovp3H4TZjVE7j3vbL/view?usp=sharing")
variant_df_for_roman = data.table::fread('VarioPath_VEP_annotation_only_variopath_transcripts_variants_20200615.tab', 
                                         sep = '\t')
# Remove downloaded file and keep inly the one stored in RAM
file.remove('VarioPath_VEP_annotation_only_variopath_transcripts_variants_20200615.tab')
# Retain just the missense variants that are the one Roman will use in its analysis
variant_df_for_roman <- variant_df_for_roman %>% 
  filter(Consequence == 'missense_variant')
#check if label with 'chr' at the beginning, if not add it
if (unique(startsWith(variant_df_for_roman$Uploaded_variation, 'chr')) == FALSE) {
  variant_df_for_roman$Uploaded_variation = paste0("chr", variant_df_for_roman$Uploaded_variation)
}
# Import the variatns used for the MDTs; these are stored locally because of pobbibly sensible data.
MDT_variants = data.table::fread('../MDT_master_spreadsheet_to_clean.tab', 
                                 sep = '\t', 
                                 select = c('NEW_ID','UKB_AF'), 
                                 stringsAsFactors = FALSE)
# extract AF calculate in UKB for the MDT variants
MDT_variants[,'AF_200K_UKB'] = as.numeric(extract.feature(MDT_variants, 'UKB_AF', 'AF'))
#check if label with 'chr' at the beginning, if not add it
if (unique(startsWith(MDT_variants$NEW_ID, 'chr')) == FALSE) {
  MDT_variants$NEW_ID = paste0("chr", MDT_variants$NEW_ID)
}


variant_df_for_roman$MDT_variant = 'No'
for (id in seq(1,length(variant_df_for_roman$MDT_variant))) {
  if (variant_df_for_roman$Uploaded_variation[id] %in% MDT_variants$NEW_ID) {
    variant_df_for_roman$MDT_variant[id] = 'Yes' 
  }
}
# Attach AF to Roman's table
variant_df_for_roman = merge( variant_df_for_roman, MDT_variants[,c('NEW_ID','AF_200K_UKB')], 
                              by.x='Uploaded_variation', 
                              by.y='NEW_ID', 
                              all.x = TRUE )
# Assign variant bins
variant_df_for_roman$bin = assign.bin(variant_df_for_roman$AF_200K_UKB)
# Extract relevant info for Roman
variant_df_for_roman_final = variant_df_for_roman[,c('SYMBOL', 'SWISSPROT', 'Uploaded_variation', 'Protein_position', 'Amino_acids', 'MDT_variant', 'bin')]
write.table(variant_df_for_roman_final,
            'variant_for_Roman_structure_effect_analysis.tsv', 
            append = F, quote = F, sep = '\t', row.names = F, col.names = T)
# Upload in gdrive
googledrive::drive_upload(media = './variant_for_Roman_structure_effect_analysis.tsv',
                          overwrite = T,
                          name = 'variant_for_Roman_structure_effect_analysis.tsv',
                          path = "https://drive.google.com/drive/folders/1vZxRNv0m4lOhVhFHP_Cq6XOy_fz1FrbO?usp=sharing")
file.remove('variant_for_Roman_structure_effect_analysis.tsv')



