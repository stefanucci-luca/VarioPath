
# load list od VTE MDT variants.
# load packages
library(tidyverse)
library(data.table)

# Functions:
#extaract cases IDs
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

# VTE MDT sheet
df_vte = xlsx::read.xlsx("~/Desktop/VarioPath/MDT_variants/V20210211/MDT_for_thrombosis.xls", sheetIndex = 1)

# UKB participants vte info
load("/Users/luca/Library/Group Containers/G69SCX94XU.duck/Library/Application Support/duck/Volumes/logincpu.hpc.cam.ac.uk – SFTP/shared_luanluan_luca_UKB/Luanluan/data/ukb_dvt_pe_combined.rdata")
df$eid = as.character(df$eid)

# create new dt with info for analysis
dt_vte = data.table( "eid" = get.ids(df_vte$PARTECIPANTS), 
                     "case" = T)

# get variant info
dt_vte$variant = lapply(dt_vte$eid, 
       function(x){
         paste(df_vte[grep(x, df_vte$PARTECIPANTS),"CHROM"],
               df_vte[grep(x, df_vte$PARTECIPANTS),"POS"],
               df_vte[grep(x, df_vte$PARTECIPANTS),"REF"],
               df_vte[grep(x, df_vte$PARTECIPANTS),"ALT"],
               sep="_")
          }
       )

dt_vte = dt_vte %>% separate(variant, sep = ',', into = c("variant_1","variant_2"))

dt_vte$variant_1 = str_remove_all(dt_vte$variant_1, "c\\(")
dt_vte$variant_2 = str_remove_all(dt_vte$variant_2, "\\)")

write.csv(dt_vte, file = "vte_id_variants.csv", quote = F, row.names = F)
dt_vte = read.csv(file = "vte_id_variants.csv", header = T)
dt_vte$eid = as.character(dt_vte$eid)
file.remove("vte_id_variants.csv")

# UKB genetic sex
df_sex = fread("/Users/luca/Library/Group Containers/G69SCX94XU.duck/Library/Application Support/duck/Volumes/logincpu.hpc.cam.ac.uk – SFTP/shared_luanluan_luca_UKB/ukbcvo/sex_otUKBcnv.tab")[,c("f.eid","f.22001.0.0")]
df_sex$f.eid = as.character(df_sex$f.eid)

# merge info
dt_all_info = merge(df_sex, dt_vte, by.x = "f.eid", by.y = "eid", all = T)
dt_all_info$case = replace_na(dt_all_info$case, FALSE)
dt_all_info = merge(dt_all_info, df, by.x = "f.eid", by.y = "eid")

# Male analysis
dt_male = as.data.table(dt_all_info[which(dt_all_info$f.22001.0.0 == 1),])

parameter = c()
logodds = c()
SE = c()
ci = c()
p = c()
gender = c()
var = c()

for (variant in unique(dt_male$variant_1) ) {
  tmp_case = dt_male[which(dt_male$variant_1 == variant),]
  if (dim(tmp_case)[1]>=10){
      tmp_control = dt_male[which(is.na(dt_male$variant_1)),]
      tmp = rbind(tmp_case, tmp_control)
      tmp$variant_1 = as.character(tmp$variant_1)
      tmp$variant_1[which(is.na(tmp$variant_1))] = "no_variant"
      model = glm(data = tmp, vte ~ variant_1 + bmi + ages, family = "binomial")
      par = parameters::parameters(model = model, df_method="wald")
      var = c(var, variant)
      parameter = c(parameter, par$Parameter[2])
      logodds = c(logodds, par$Coefficient[2])
      SE = c(SE, par$SE[2])
      ci = c(ci, paste(par$CI_low[2], par$CI_high[2], sep = "; "))
      p = c(p, par$p[2])
      gender = c(gender, "Male")
    }
  }

# Feale analysis
dt_female = as.data.table(dt_all_info[which(dt_all_info$f.22001.0.0 == 0),])


for (variant in unique(dt_female$variant_1) ) {
  tmp_case = dt_female[which(dt_female$variant_1 == variant),]
  if (dim(tmp_case)[1]>=10){
    tmp_control = dt_female[which(is.na(dt_male$variant_1)),]
    tmp = rbind(tmp_case, tmp_control)
    tmp$variant_1 = as.character(tmp$variant_1)
    tmp$variant_1[which(is.na(tmp$variant_1))] = "no_variant"
    model = glm(data = tmp, vte ~ variant_1 + bmi + ages, family = "binomial")
    par = parameters::parameters(model = model, df_method="wald")
    var = c(var, variant)
    parameter = c(parameter, par$Parameter[2])
    logodds = c(logodds, par$Coefficient[2])
    SE = c(SE, par$SE[2])
    ci = c(ci, paste(par$CI_low[2], par$CI_high[2], sep = "; "))
    p = c(p, par$p[2])
    gender = c(gender, "Female")
  }
}

df_exp = data.frame( 
  'variant' = var,
  'predictor' = parameter,
  'log(OR)' = logodds, 
  'SE' = SE,
  'ConfInt[95]' = ci,
  'p-val' = p,
  'sex' = gender )

write.csv(df_exp, "single_effect_vte.csv", quote = F, row.names = F)

# Not gender stratified 

for (variant in unique(dt_all_info$variant_1) ) {
  tmp_case = dt_all_info[which(dt_all_info$variant_1 == variant),]
  if (dim(tmp_case)[1]>=10){
    tmp_control = dt_all_info[which(is.na(dt_male$variant_1)),]
    tmp = rbind(tmp_case, tmp_control)
    tmp$variant_1 = as.character(tmp$variant_1)
    tmp$variant_1[which(is.na(tmp$variant_1))] = "no_variant"
    model = glm(data = tmp, vte ~ variant_1 + bmi + ages, family = "binomial")
    par = parameters::parameters(model = model, df_method="wald")
    var = c(var, variant)
    parameter = c(parameter, par$Parameter[2])
    logodds = c(logodds, par$Coefficient[2])
    SE = c(SE, par$SE[2])
    ci = c(ci, paste(par$CI_low[2], par$CI_high[2], sep = "; "))
    p = c(p, par$p[2])
    gender = c(gender, "no gender stratified")
  }
}

df_exp = data.frame( 
            'variant' = var,
            'predictor' = parameter,
            'log(OR)' = logodds, 
            'SE' = SE,
            'ConfInt[95]' = ci,
            'p-val' = p,
            'sex' = gender )

write.csv(df_exp, "single_effect_vte.csv", quote = F, row.names = F)
