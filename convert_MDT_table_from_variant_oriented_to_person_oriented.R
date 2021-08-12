library(dplyr)
library(stringr)
library(data.table)
library(tidyverse)

dt_no_AF = data.table::fread("/Users/ls31/Library/Group Containers/G69SCX94XU.duck/Library/Application Support/duck/Volumes/logincpu.hpc.cam.ac.uk – SFTP/rds/rds-who1000-cbrc/data/UKBb/variopath/VCF_200K_filtered_for_variopath/no_AF_filter/All_variopath_variants_in_unrelated_european_noAFfilter_VEPed.tsv", 
                  nThread = 5)

dt_platelet = xlsx::read.xlsx("/Volumes/GoogleDrive-105684671225474146017/My Drive/Desktop_Macbook_PhD/Desktop/VarioPath/MDT_variants/V20210211/MDT_for_platelet.xls", sheetIndex = 1)
dt_thrombo = xlsx::read.xlsx("/Volumes/GoogleDrive-105684671225474146017/My Drive/Desktop_Macbook_PhD/Desktop/VarioPath/MDT_variants/V20210211/MDT_for_thrombosis.xls", sheetIndex = 1)
dt_hs = xlsx::read.xlsx("/Volumes/GoogleDrive-105684671225474146017/My Drive/Desktop_Macbook_PhD/Desktop/VarioPath/MDT_variants/V20210211/MDT_for_hereditary_spherocytosis.xls", sheetIndex = 1)
dt_bleeding = xlsx::read.xlsx("/Volumes/GoogleDrive-105684671225474146017/My Drive/Desktop_Macbook_PhD/Desktop/VarioPath/MDT_variants/V20210211/MDT_for_bleeding_and_coagulation.xls", sheetIndex = 1)

dt_mdt = do.call(what = rbind,args = list(dt_platelet,dt_thrombo, dt_hs, dt_bleeding))

dt_mdt = dt_mdt %>% 
  arrange(GENE, CHROM, POS, REF, ALT, .by_group = TRUE)

dt_mdt_reshaped = summarise(dt_mdt, as_tibble(str_split(dt_mdt$PARTECIPANTS, pattern = " ", simplify = T)))
colnames(dt_mdt_reshaped) <- str_replace(colnames(dt_mdt_reshaped), pattern = "V", replacement = "partecipant_")

dt_mdt_reshaped = cbind(dt_mdt, dt_mdt_reshaped)

dt_melted = reshape2::melt(data = dt_mdt_reshaped,
                           id.vars = colnames(select(
                             dt_mdt_reshaped,!starts_with(match = "partecipant_")
                           ))) 

dt_melted_clean = dt_melted[-which(dt_melted$value == ""), ] %>% 
  arrange(GENE, CHROM, POS, REF, ALT, .by_group = TRUE)

participant_info = as_tibble(str_split(string = dt_melted_clean$value,pattern = "=|;",simplify = T,n = 3))
colnames(participant_info) = c("participantID","genotype","reads_supproting_variant")

dt_melted_clean_expanded = cbind(dt_melted_clean, participant_info)


# European list
dt_euro = read_tsv("/Users/ls31/Library/Group Containers/G69SCX94XU.duck/Library/Application Support/duck/Volumes/logincpu.hpc.cam.ac.uk – SFTP/rds/rds-who1000-cbrc/data/shared_luanluan_luca_UKB/unrelated_ethnic_groups/unrelated_UKB_200_ethnicity.txt.somalier-ancestry_not_1KGP_EUR.tsv",
                   col_names = F) %>% 
  select(1)
# Keep only unrelated european 
dt_melted_clean_expanded_unr_euro = dt_melted_clean_expanded[which(dt_melted_clean_expanded$participantID %in% dt_euro$X1),]


as.data.frame(sort(table(unique(dt_melted_clean_expanded[,c('GENE', 'CHROM', 'POS', 'REF', 'ALT',"participantID")])$participantID), decreasing = T))


multivar_people = c('1163990', '1686054', '3030231', '1256114')
View(dt_melted_clean_expanded[which(dt_melted_clean_expanded$participantID %in% multivar_people),])

write_csv(dt_melted_clean_expanded, "/Volumes/GoogleDrive-105684671225474146017/My Drive/Desktop_Macbook_PhD/Desktop/VarioPath/MDT_variants/spreadsheet_reordered_per_person.csv",
          quote_escape = "none", 
          col_names = T, 
          append = F)
