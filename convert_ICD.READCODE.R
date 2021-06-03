

conver.csv.ICDtoRC = function(file_csv){
  vect_icd10_to_rc = df_conv$read_code
  names(vect_icd10_to_rc) = str_remove(df_conv$icd10_code, "\\.")
  df_tmp = suppressWarnings(
    read.csv(file_csv, header = T)
  )
  df_tmp$read_codes = suppressMessages( 
    plyr::revalue(str_remove(df_tmp$ICD10code, "\\.")
                                    , vect_icd10_to_rc)
  )
  return(df_tmp)
  }


dr_in = path.expand('/Users/luca/Library/Group Containers/G69SCX94XU.duck/Library/Application Support/duck/Volumes/logincpu.hpc.cam.ac.uk – SFTP/shared_luanluan_luca_UKB/Luca/definitions/secondary_care')
dr_out = path.expand('/Users/luca/Library/Group Containers/G69SCX94XU.duck/Library/Application Support/duck/Volumes/logincpu.hpc.cam.ac.uk – SFTP/shared_luanluan_luca_UKB/Luca/definitions/primary_care')
df_conv = openxlsx::read.xlsx( xlsxFile = "~/Desktop/VarioPath/primarycare_codings_UKB_20210520/all_lkps_maps_v2.xlsx", sheet = 12)

for (icd_file in grep(list.files(dr_in, include.dirs = F, pattern = "*.csv"), pattern='OPCS', invert=TRUE, value=TRUE))
{
  df_converted = conver.csv.ICDtoRC(file.path( dr_in, icd_file))
  write.csv(df_converted,
            file.path(dr_out, 
                      paste0(stringr::str_remove(basename(icd_file), ".csv"), "converted_to_RC.csv" )
            )
  )
}



