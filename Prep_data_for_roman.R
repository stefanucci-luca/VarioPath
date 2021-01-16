
library(ensemblVEP)
library(dplyr)
library(stringr)



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

# Download variable list from gdrive - Variopath directory
# 
googledrive::drive_download("https://drive.google.com/file/d/1WAqePfc5qZH-mkAWtOsxrjh8ayvHEudX/view?usp=sharing")
variant_df = data.table::fread('VarioPath_variants_20200615_GENE.tsv', sep = '\t')
# set names
colnames(variant_df) <- c('CHROM',
                          'start',
                          'id',
                          'ref',
                          'alt',
                          'pathogenicity',
                          'source',
                          'INFO')
# preview df
str(variant_df)
# create a working variable
variant_df_w = variant_df
# Remove downloaded file and keep inly the one stored in RAM
file.remove('VarioPath_variants_20200615_GENE.tsv')

# For Roman only SNV variants are needed that are in coding regions.
# First I'll remove everything that is not SNV. 
variant_df_w = variant_df_w[which(check.if.SNV(variant_df_w, 'ref', 'alt')),]

# Second add the gene symbol info and column
variant_df_w$gene_symbol = extract.feature(variant_df_w, 'INFO', 'GENE')

genes_for_roman = c("F10", "F11", "F12", "F13A1", "F13B", "F2", "F5", "F7", "F8", "F9", "FGA", "FGB", 
                    "FGG", "GGCX", "KNG1", "LMAN1", "MCFD2", "SERPINE1", "SERPINF2", "VKORC1", "VWF", 
                    "ADAMTS13", "HRG", "PIGA", "PLG", "PROC", "PROS1", "SERPINC1", "SERPIND1", "THBD", 
                    "ABCC4", "ABCG5", "ABCG8", "ACTB", "ACTN1", "ANKRD26", "ANO6", "AP3B1", "AP3D1", 
                    "ARPC1B", "BLOC1S3", "BLOC1S6", "CDC42", "CYCS", "DIAPH1", "DTNBP1", "ETV6", "FERMT3", 
                    "FLI1", "FLNA", "FYB1", "GATA1", "GFI1B", "GNE", "GP1BA", "GP1BB", "GP6", "GP9", "HOXA11", 
                    "HPS1", "HPS3", "HPS4", "HPS5", "HPS6", "IKZF5", "ITGA2B", "ITGB3", "KDSR", "LYST", "MECOM", 
                    "MPIG6B", "MPL", "MYH9", "NBEA", "NBEAL2", "P2RY12", "PLA2G4A", "PLAU", "RASGRP2", "RBM8A", 
                    "RNU4ATAC", "RUNX1", "SLFN14", "SRC", "STIM1", "STXBP2", "TBXA2R", "TBXAS1", "THPO", "TUBB1", 
                    "VIPAS39", "VPS33B", "WAS", "ANK1", "EPB41", "EPB42", "SLC4A1", "SPTA1", "SPTB")

variant_df_w = dplyr::filter(variant_df_w, 
                             GENE %in% genes_for_roman)


write.table(x=variant_df_w[,c('CHROM','start','id','ref','alt','pathogenicity')],
          'tmp_var_file.vcf',
          append = F, 
          quote = F, 
          row.names = F, 
          col.names = F,
          sep = '\t')
# This bit of the script not working yet so I annotated with VEP using the online tool and selecting only coding variants 
myparam <- VEPFlags(flags=list(host="ensembldb.ensembl.org"))
ensemblVEP('tmp_var_file.vcf', myparam, verbose=TRUE)
# 
# file.remove('tmp_var_file.vcf')

googledrive::drive_download("https://drive.google.com/file/d/1kVgKbaL31WK08fAEzNxbAGpHd8LFIJRM/view?usp=sharing")
variant_df_for_roman = data.table::fread('Coding_variants_MDT_genes_for_Roman_missing_NO_GENE_variants.tsv', sep = '\t')
# Remove downloaded file and keep inly the one stored in RAM
file.remove('Coding_variants_MDT_genes_for_Roman_missing_NO_GENE_variants.tsv')

# Get variatns used for the MDTs
MDT_variants = data.table::fread('../MDT_master_spreadsheet_to_clean.tab', 
                                 sep = '\t', 
                                 select = c('NEW_ID'), 
                                 stringsAsFactors = FALSE)

MDT_variants$NEW_ID_2 = ''
for (i in 1:dim(MDT_variants)[1]) {
  MDT_variants$NEW_ID_2[i] = str_replace_all( str_remove_all(MDT_variants$NEW_ID[i], 
                                                          'chr'),
                                           "_",
                                           ":")
  }

variant_df_for_roman$MDT_variant = 'No'
for (id in seq(1,length(variant_df_for_roman$MDT_variant))) {
  if (variant_df_for_roman$V1[id] %in% MDT_variants$NEW_ID_2) {
    variant_df_for_roman$MDT_variant[id] = 'Yes' 
  }
}

feat_to_extract = c('SWISSPROT', 'SYMBOL=')

for (feat in feat_to_extract) {
  variant_df_for_roman[,feat] = extract.feature(variant_df_for_roman, 'V14', feat)
}
