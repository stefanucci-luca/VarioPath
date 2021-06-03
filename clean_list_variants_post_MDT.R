library(tidyverse)
library(googledrive)

lis  = list(
  "MDT_for_thrombosis_2021.02.08.xlsx" = 'https://drive.google.com/file/d/1l2LEblI8VI19XnrPJm2uB5_OG_cLhfux/view?usp=sharing',
  "MDT_for_platelet_20210226_MDT_round1_completed.xlsx" = 'https://drive.google.com/file/d/15nWkN5HzdJgy8_yILdYEQc42hNCEsclK/view?usp=sharing',
  "Copy of MDT_for_bleeding_and_coagulation_20210128_MS_update_latest_JC - 28 January, 22:29.xlsx" = 'https://drive.google.com/file/d/1rm3N0U9rQqtNOI-Z8eXY-LphBiXvH4x4/view?usp=sharing'
)

for (i in 1:length(lis)) {
  googledrive::drive_download(as.character(lis[i]))
  df_tmp = xlsx::read.xlsx(names(lis[i]), sheetIndex = 1)
  file.remove(names(lis[i]))
  assign( paste0("df",i),
          value = df_tmp
          )
}

col_sel = c("MDT.decision", "GENE", "CHROM", "POS", "REF", "ALT" )

df1 = df1 %>% select(col_sel)
df2 = df2 %>% select(col_sel)
colnames(df3)[1] = "MDT.decision"
df3 = df3 %>% select(col_sel)

df_final = do.call(rbind, list(df1,df2,df3)) %>% 
  na.omit() %>% 
  filter(MDT.decision == "2. Accept but review pheno" | MDT.decision == "1. Accept" | MDT.decision == "3. Undecided" |
           MDT.decision == "Accept but review phenotype" | MDT.decision == "Accept" | MDT.decision == "Accept but review pheno" |
           MDT.decision == "Undecided")



xlsx::write.xlsx(df_final, file = paste0("/Volumes/GoogleDrive/My Drive/Phd/VarioPath/MDT/mdt_pass_variatns_", 
                                         format(Sys.Date(), "%Y%m%d"), ".xlsx" 
                                         ), 
                 sheetName = "MDT_pass", 
                 row.names = F
                                )

df_final$dot = "."
df_final_hpc = df_final %>% select('CHROM',	'POS', 'dot',	'REF',	'ALT')


write_delim(df_final_hpc, file = paste0("/Users/luca/Library/Group Containers/G69SCX94XU.duck/Library/Application Support/duck/Volumes/logincpu.hpc.cam.ac.uk – SFTP/home/ls760/UKBb/variopath/burden_test_raremetal/mdt_pass_variatns_", 
                                         format(Sys.Date(), "%Y%m%d"),".csv"), quote_escape = F, delim = "\t")
