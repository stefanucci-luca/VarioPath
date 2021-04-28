library(tidyverse)

df_hearing = data.table::fread("/Users/luca/Library/Group Containers/G69SCX94XU.duck/Library/Application Support/duck/Volumes/logincpu.hpc.cam.ac.uk – SFTP/home/ls760/shared_luanluan_luca_UKB/ukbcvo/hearin_test_otUKBcnv.tab")
df_blobio = data.table::fread("/Users/luca/Library/Group Containers/G69SCX94XU.duck/Library/Application Support/duck/Volumes/logincpu.hpc.cam.ac.uk – SFTP/home/ls760/shared_luanluan_luca_UKB/ukbcvo/blod_biochemistry_otUKBcnv.tab")
df_uri = data.table::fread("/Users/luca/Library/Group Containers/G69SCX94XU.duck/Library/Application Support/duck/Volumes/logincpu.hpc.cam.ac.uk – SFTP/home/ls760/shared_luanluan_luca_UKB/ukbcvo/urine_test_otUKBcnv.tab")

df_complete = Reduce(function(x,y) merge(x = x, y = y, by = "f.eid"), 
                     list(df_hearing, df_blobio, df_uri))

switch_name = c( 
        "f.eid" = "ID",
        "f.30500." = "Microalbumin in urine",
        "f.30505." = "Microalbumin flag. Result below 6.7 mg/L",
        "f.30510." = "Creatinine (enzymatic) in urine micromole/L",
        "f.30520." = "Potassium in urine millimole/L",
        "f.30525." = "Potassium in urinie flag. Result below 2 millimole/L or above 200 millimole/L",
        "f.30530." = "Sodium in urine millimole/L",
        "f.30535." =  "Sodium in Urine flag. Result below 10 millimole/L or Result above 400 millimole/L",
        "f.30600." = "Albumin",
        "f.30610." = "Alkaline phosphatase",
        "f.30620." = "Alanine aminotransferase U/L",
        "f.30630." = "Apolipoprotein A g/L",
        "f.30640." = "Apolipoprotein B g/L",
        "f.30650." = "Aspartate aminotransferase U/L",
        "f.30660." = "Direct bilirubin umol/L",
        "f.30670." = "Urea mmol/L",
        "f.30680." = "Calcium mmol/L",
        "f.30690." = "Cholesterol mmol/L",
        "f.30700." = "Creatinine umol/L",
        "f.30710." = 'C-reactive protein mg/L',
        "f.30720." = "Cystatin C mg/L",
        "f.30730." = "Gamma glutamyltransferase U/L",
        "f.30740." = "Glucose mmol/L",
        "f.30750." = "HbA1c mmol/mol",
        "f.30760." = "HDL cholesterol mmol/L",
        "f.30770." = "IGF-1 nmol/L",
        "f.30780." = "LDL direct mmol/L",
        "f.30790." = "Lipoprotein A nmol/L",
        "f.30800." = "Oestradiol pmol/L",
        "f.30810." = "Phosphate mmol/L",
        "f.30820." = "Rheumatoid factor IU/ml",
        "f.30830." = "SHBG nmol/L",
        "f.30840." = "Total bilirubin umol/L",
        "f.30850." = "Testosterone nmol/L",
        "f.30860." = "Total protein g/L",
        "f.30870." = "Triglycerides mmol/L",
        "f.30880." = "Urate umol/L",
        "f.30890." = "Vitamin D nmol/L",
        "f.4233." = "Signal to noise ratio left",
        "f.4244." = "Signal to noise ratio right",
        "f.20019." = "Speech-reception treshold left",
        "f.20021." = "Speech-reception treshold right"
  )


for (repl in 1:length(switch_name)) {
  position = grep(pattern = names(switch_name)[repl], 
              x = colnames(df_complete), )
  valUK = gsub(pattern = names(switch_name)[repl], 
               replacement = switch_name[repl], 
               x = colnames(df_complete), )
  colnames(df_complete)[position] = valUK[position]
}

saveRDS(df_complete, "/Volumes/GoogleDrive/My Drive/Phd/VarioPath/MDT/MDT Platelet/Phenotype/phenotype_biochem_urine_hearing_20210428.RDS")
