#!/usr/bin/Rstudio

#import the libraries
library(readr)
library(dplyr)
library(reshape)
library(stringr)
library(tidyr)
library(pbapply)
library(googledrive)

# read in the variant list
drive_download("https://drive.google.com/file/d/1zWYevop2cHk3-8wWqCofSZTAsiZkfjsP/view?usp=sharing" )
df_vk = data.table::fread("ALL_variants_P_LP_someVUS_original_karyn_HGNC.txt", header = T)

# Removed the "Chr"
df_vk_clean <- df_vk[!which(df_vk$Chr == "Chr"),]
# Replace Chr23 with ChrX
df_vk_clean$Chr[which(df_vk_clean$Chr == 23)] = "X"

# Fix the headers to resemble VCF structure
# replace the third column with openCGA 
colnames(df_vk_clean)[1] <- "#CHROM"
colnames(df_vk_clean)[2] <- "POS"
colnames(df_vk_clean)[3] <- "openCGA_ID"
colnames(df_vk_clean)[6] <- "QUAL"
colnames(df_vk_clean)[7] <- "FILTER"

df_vk_clean[which(df_vk_clean$ALT == "<DEL>"),"ALT"] = "-"

df_vk_clean$openCGA_ID <- stringr::str_c(df_vk_clean$`#CHROM`, 
                                         df_vk_clean$POS, 
                                         df_vk_clean$REF, 
                                         df_vk_clean$ALT
                                         ,sep = ":")

#________________________________________________________________
# Change the source names and reshape it do remove duplicate entries
# split the sources in different colums
df_vk_clean_tmp <- tidyr::separate(
  data = df_vk_clean,
  col = FILTER,
  into = paste("subfilter", seq(1:4), sep = "_"), # this procudes a warning, because are some missing columns, but the last column is an empty char
  sep = ";",
  fill = "left",
  remove = F
)

df_vk_clean_tmp <- as.data.frame(df_vk_clean_tmp)

# for every column containing the sources change the entries with the new one.
for (i in paste("subfilter", seq(1:4), sep = "_")){
  position = which(colnames(df_vk_clean_tmp) == i)
  print(i)
  print(position)
  df_vk_clean_tmp[, position] <-  plyr::revalue( df_vk_clean_tmp[, position],
                                                 c(
                                                   "hgmdAll" = "curatedPub",
                                                   "hgmdVWF" = "curatedPub",
                                                   "curatedMIx" = "curatedNBR",
                                                   "ClinVarAll" = "ClinVar",
                                                   "ClinVarVWF" = "ClinVar"
                                                 )
  )
}
# reshape and attach to the main dataframe
df_vk_clean$FILTER <-  
  apply(df_vk_clean_tmp[,c("subfilter_1",
                           "subfilter_2",
                           "subfilter_3",
                           "subfilter_4")], 
        1, 
        function(x) 
          paste(str_replace_na(unique(x), "")[str_replace_na(unique(x), "") != ""], collapse ="; ")
  )

#____________________________________________________________________
# Consolidate the info columns

cols = c("INFO","INFO_tag","V10")
df_vk_clean_tmp$information_col = do.call(paste, c(df_vk_clean_tmp[,cols], sep=";"))

patterns = c( 'hgmdAll|hgmdVWF|curatedMIx|ClinVarAll|ClinVarVWF| ' )

df_vk_clean$INFO = lapply(df_vk_clean_tmp[,"information_col"], function(x) 
  sort(str_replace_all( unique(unlist(str_split(x, ":|;"))),
                        pattern = patterns ,
                        replacement = "" )[str_replace_all( unique(unlist(str_split(x, ":|;"))),
                                                            pattern = patterns ,
                                                            replacement = "" ) != ''])
)

# Clean the df from the unwanted columms
df_vk_clean = df_vk_clean[, !c("INFO_tag", "V10") ]

#____________________________________________________________________
# use the QUAL column to keep track of the variant pathogenicity

df_vk_clean$QUAL =  pbapply(df_vk_clean[,], 
                          1, 
                          function(x) 
                          unique(str_remove_all(
                                    gsub(".*=",
                                      "",
                                      grep('CLASS|CLNSIG', 
                                           x = str_split(string = x, pattern = ",|:|;", simplify = T), 
                                           value = TRUE)
                                      ), pattern = '"'
                                 )
                          )
                  )

#Remove from CLASS and CLNSIF from the info column

patterns = c( 'CLASS=DM|CLASS=DM?|CLNSIG=Pathogenic|
               CLINSIG=Pathogenic/Likely_pathogenic|
               CLINSIG=Likely_pathogenic|
               CLINSIG=not_provided|CLINSIG=Uncertain_significance|
               CLINSIG=Conflicting_interpretations_of_pathogenicity|
               CLINSIG=Likely_pathogenic(1)|\\\\n|\\n' )
# 
df_vk_clean$INFO = pbapply(df_vk_clean[,"INFO"], 
                            1, 
                            function(x) 
                            str_remove_all(sort(str_replace_all(str_split(x, ":|;"),
                                                 pattern = patterns ,
                                                 replacement = "" )
                                                  ), 
                                pattern = '\\?|"|\\\\')
                        )
#____________________________________________________________________
# write the output
data.table::fwrite(df_vk_clean, 
                   "VarioPath_variants_20200615.tsv", sep = '\t', 
                   quote = F,
                   row.names = F, 
                   col.names = T, qmethod = "escape", 
)
#____________________________________________________________________
#Upload on the variopath folder in gdrive
df_vk_upload <- drive_upload(media = "VarioPath_variants_20200615.tsv",
             path = as_id("https://drive.google.com/drive/folders/1Eu5NNuODO3uYIozSPRE3oHJSoHiUkD81?usp=sharing"),
             overwrite = T,
             name = "VarioPath_variants_20200615.tsv")
# Assign the permission
# df_vk_ipload <- df_vk_ipload %>%
#   drive_share(role = "reader", type = "group")

#___________________________________________________________________
#Remove files from local machine
file.remove("ALL_variants_P_LP_someVUS_original_karyn_HGNC.txt")
file.remove("VarioPath_variants_20200615.tsv")
