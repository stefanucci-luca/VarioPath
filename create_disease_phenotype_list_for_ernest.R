
df_diseases_prenotypes = readxl::read_xlsx('/Users/luca/Library/Group Containers/G69SCX94XU.duck/Library/Application Support/duck/Volumes/logincpu.hpc.cam.ac.uk – SFTP/home/ls760/shared_luanluan_luca_UKB/Luca/definitions/dictionary_11_March_2021.xlsx', sheet = 1) 
df_conversion =  read.csv('/Users/luca/Library/Group Containers/G69SCX94XU.duck/Library/Application Support/duck/Volumes/logincpu.hpc.cam.ac.uk – SFTP/home/ls760/shared_luanluan_luca_UKB/Luca/returned_datasets/phenotype_event_names.csv') 
load('/Users/luca/Library/Group Containers/G69SCX94XU.duck/Library/Application Support/duck/Volumes/logincpu.hpc.cam.ac.uk – SFTP/home/ls760/shared_luanluan_luca_UKB/Luca/returned_datasets/ukb_first_events.rdata')

retrun.disease.code.format = function(num){
  return(sprintf("%03d", num))
}

for (dis_code in df_conversion[,3]) {
  column =  paste0( 'cal_ps_d_', retrun.disease.code.format(dis_code))
  file_name = as.character( df_conversion[dis_code,2] )
  participants_vector = paste(ukb_first_events[ which(is.na(ukb_first_events[,column]) == FALSE ) ,c("identifier",column)][,1], collapse = "\n")
  write.table(participants_vector, 
            file = file.path ( "Desktop/VarioPath/bevimed/phenotype_definition" , 
                               paste0( "participant_list_", 
                                       file_name, ".txt") ), 
            quote = FALSE, append = F, row.names = F, col.names = F  )
}


ukb_first_events[which(ukb_first_events$identifier == '1308252'),]
