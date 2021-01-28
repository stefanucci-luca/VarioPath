
library(dplyr)
library(ggplot2)


unrelated.european = read.csv('/Users/luca/Library/Group Containers/G69SCX94XU.duck/Library/Application Support/duck/Volumes/logincpu.hpc.cam.ac.uk – SFTP/rds/rds-who1000-cbrc/data/shared_luanluan_luca_UKB/unrelated_ethnic_groups/unrelated_UKB_200_ethnicity.txt.somalier-ancestry_not_1KGP_EUR.tsv',
                              sep = '\t', header = F) %>% 
  select(V1)

dim(unrelated.european)
str(unrelated.european)

# extract the ID information from the partecipants column
get.ids = function(partecipantID_vector){
  require(stringr)
  id = unlist(
    lapply(
      lapply(partecipantID_vector, 
             function(x)
               str_extract_all(x , pattern = "[^=| |;]+") 
      ), 
      '[[', 
      1)
    ) 
  id_clean = id[id %>% str_length() > 4]
  return(id_clean)
}

MDTs = c('bleeding_and_coagulation',
         'hereditary_spherocytosis',
         'platelet',
         'thrombosis'
)

setwd('/Users/luca/Desktop/VarioPath/MDT_variants')
latest_dir = '/Users/luca/Library/Group Containers/G69SCX94XU.duck/Library/Application Support/duck/Volumes/logincpu.hpc.cam.ac.uk – SFTP/cases_MDT_V20210126'


cases_total = c()
for (MDT in MDTs) {
  
  outdir='/Users/luca/Desktop/VarioPath/MDT_variants/cases'
  file_name_out_1 = paste("cases_only_list_", MDT, ".tab", sep = '')
  file_in = file.path( latest_dir, paste( 'cases_for_', MDT, '.tsv',
                                          sep=''))
  
  MDT_cases = read.delim(file_in, header = T, sep = '\t') %>% 
    dplyr::select("GENE", "CHROM", "POS", "REF", "ALT", "PARTECIPANTS")
  
  cases = unique(get.ids(MDT_cases$PARTECIPANTS))
  cases = cases[unique(get.ids(MDT_cases$PARTECIPANTS)) %in% unrelated.european$V1]
  
  write.csv(cases, 
            file = paste0(latest_dir, '/', file_name_out_1), 
            row.names = F,
            quote = F)

  cases_total = c(cases, cases_total)
  
}

control_total = unrelated.european$V1[!unrelated.european$V1 %in% unique(cases_total)]
 
write.csv(control_total, 
          file = paste0(latest_dir, '/', 'control_partecipants_EUR_unrelated.tab'), 
          row.names = F,
          quote = F)

