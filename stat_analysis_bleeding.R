library(tidyverse)

load('/Users/luca/Library/Group Containers/G69SCX94XU.duck/Library/Application Support/duck/Volumes/logincpu.hpc.cam.ac.uk – SFTP/rds/project/who1000-1/rds-who1000-cbrc/data/shared_luanluan_luca_UKB/Luca/returned_datasets/ukb_first_events.rdata')
load('/Users/luca/Library/Group Containers/G69SCX94XU.duck/Library/Application Support/duck/Volumes/logincpu.hpc.cam.ac.uk – SFTP/rds/project/who1000-1/rds-who1000-cbrc/data/shared_luanluan_luca_UKB/Luanluan/data/ukb_dvt_pe_combined.rdata')
var_info = readxl::read_xls('Desktop/VarioPath/MDT_variants/V20210211/MDT_for_bleeding_and_coagulation.xls')
var_info = unique(var_info)

participants = unlist(str_split(var_info$PARTECIPANTS, pattern = " "))

position_variant_per_participant = unlist(
  lapply(participants, 
         function(x)
           grep(pattern = x, x = var_info$PARTECIPANTS )))


df_participants = as.data.frame(
  lapply(list(1:length(participants)), 
         function (x) 
           cbind( participants[x], var_info[position_variant_per_participant[x],] )
  ))

df_participants$participants = unlist(lapply(str_split(string = as.character(df_participants$participants.x.), pattern = '=' ),'[[', 1))

df_phenotypes = merge( ukb_first_events, df_participants,  by.x ='identifier', by.y = 'participants', all.x = TRUE )

df_phenotypes = merge( df_phenotypes, df,  by.x ='identifier', by.y = 'eid', all.x = TRUE )

# Vector to rename columns
#####################################
# Rename columns
v1= c('cal_p_d_001' =	'Acute_lymphoblastic_leukaemia_p_p',
'cal_p_d_002'=	'Acute_myeloid_leukaemia_p_p',
'cal_p_d_003'=	'AF_p_p',
'cal_p_d_004'=	'Aplastic_anaemia_p',
'cal_p_d_005'=	'Bleeding_after_tooth_extraction_or_other_procedure_p',
'cal_p_d_006'=	'CALIBER_DVT_p',
'cal_p_d_007'=	'CALIBER_PE_p',
'cal_p_d_008'=	'Cataract_p',
'cal_p_d_009'=	'Central_nervous_system_bleeding_p',
'cal_p_d_010'=	'Cutaneous_bruising_or_bleeding_p',
'cal_p_d_011'=	'Epistaxis_p',
'cal_p_d_012'=	'Gastrointestinal_bleeding_p',
'cal_p_d_013'=	'Haemangioma_p',
'cal_p_d_014'=	'Haemarthrosis_p',
'cal_p_d_015'=	'Hearing_loss_p',
'cal_p_d_016'=	'Hereditary_deficiency_of_other_clotting_factors_p',
'cal_p_d_017'=	'Hereditary_factor_IX_deficiency_p',
'cal_p_d_018'=	'Hereditary_factor_VIII_deficiency_p',
'cal_p_d_019'=	'Hereditary_factor_XI_deficiency_p',
'cal_p_d_020'=	'IHD_p',
'cal_p_d_021'=	'ISSTR_p',
'cal_p_d_022'=	'Menorrhagia_p',
'cal_p_d_023'=	'Myelodysplastic_syndrome_p',
'cal_p_d_024'=	'Non-haematological_malignancy_p',
'cal_p_d_025'=	'Oral_cavity_bleeding_p',
'cal_p_d_026'=	'Other_cardiac_bleeding_p',
'cal_p_d_027'=	'Other_eye_bleeding_p',
'cal_p_d_028'=	'Other_gynaecological_bleeding_p',
'cal_p_d_029'=	'Other_haematological_malignancy_p',
'cal_p_d_030'=	'Other_haematuria_p',
'cal_p_d_031'=	'Other_obstetric_bleeding_p',
'cal_p_d_032'=	'Other_respiratory_system_bleeding_p',
'cal_p_d_033'=	'Other_traumatic_bleeding_p',
'cal_p_d_034'=	'PAD_p',
'cal_p_d_035'=	'Postpartum_haemorrhage_p',
'cal_p_d_036'=	'primary_thrombocytopenia_p',
'cal_p_d_037'=	'Qualitative_platelet_defects_p',
'cal_p_d_038'=	'Renal_disease_p',
'cal_p_d_039'=	'secondary_thrombocytopenia_p',
'cal_p_d_040'=	'Stroke_NOS_p',
'cal_p_d_041'=	'T2D_p',
'cal_p_d_042'=	'thrombocytopenia_unspecified_p',
'cal_p_d_043'=	'thyroid_p',
'cal_p_d_044'=	'TIA_p',
'cal_p_d_045'=	'von_Willebrand_disease_p',
'cal_p_d_046'=	'PID_allergic_purpura_p',
'cal_p_d_047'=	'PID_alopecia_p',
'cal_p_d_048'=	'PID_ankylosing_spondilytis_p',
'cal_p_d_049'=	'PID_autoimmune_haemolytic_anemia_p',
'cal_p_d_050'=	'PID_autoimmune_hepatatis_p',
'cal_p_d_051'=	'PID_autoimmune_polyglandular_failure_p',
'cal_p_d_052'=	'PID_blistering_disorder_p',
'cal_p_d_053'=	'PID_Coeliac_disease_p',
'cal_p_d_054'=	'PID_Idiopathic_thrombocytopenic_purpura_p',
'cal_p_d_055'=	'PID_juvenile_RA_p',
'cal_p_d_056'=	'PID_lupus_erythematosus_p',
'cal_p_d_057'= 'PID_MS_p',
'cal_p_d_058'= 'PID_other_connective_tissue_p',
'cal_p_d_059'=	'PID_psoriasis_p',
'cal_p_d_060'=	'PID_RA_p',
'cal_p_d_061'=	'PID_reactive_arthropaties_p',
'cal_p_d_062'=	'PID_Sarcoidosis_p',
'cal_p_d_063'=	'PID_SLE_p',
'cal_p_d_064'= 'PID_systemic_sclerosis_p',
'cal_p_d_065'=	'PID_T1D_p',
'cal_p_d_066'=	'PID_thyroiditis_p',
'cal_p_d_067'=	'PID_thyrotoxicosis_p',
'cal_p_d_068'=	'PID_vasculitis_p',
'cal_p_d_069'=	'PID_vitiligo_p',
'cal_p_d_070'=	'PID_IBD_p',
'cal_p_d_071'=	'PID_chron_p',
'cal_p_d_072'=	'PID_ulcerative_colitis_p',
'cal_ps_d_001' =	'Acute_lymphoblastic_leukaemia',
'cal_ps_d_002'=	'Acute_myeloid_leukaemia',
'cal_ps_d_003'=	'AF',
'cal_ps_d_004'=	'Aplastic_anaemia',
'cal_ps_d_005'=	'Bleeding_after_tooth_extraction_or_other_procedure',
'cal_ps_d_006'=	'CALIBER_DVT',
'cal_ps_d_007'=	'CALIBER_PE',
'cal_ps_d_008'=	'Cataract',
'cal_ps_d_009'=	'Central_nervous_system_bleeding',
'cal_ps_d_010'=	'Cutaneous_bruising_or_bleeding',
'cal_ps_d_011'=	'Epistaxis',
'cal_ps_d_012'=	'Gastrointestinal_bleeding',
'cal_ps_d_013'=	'Haemangioma',
'cal_ps_d_014'=	'Haemarthrosis',
'cal_ps_d_015'=	'Hearing_loss',
'cal_ps_d_016'=	'Hereditary_deficiency_of_other_clotting_factors',
'cal_ps_d_017'=	'Hereditary_factor_IX_deficiency',
'cal_ps_d_018'=	'Hereditary_factor_VIII_deficiency',
'cal_ps_d_019'=	'Hereditary_factor_XI_deficiency',
'cal_ps_d_020'=	'IHD',
'cal_ps_d_021'=	'ISSTR',
'cal_ps_d_022'=	'Menorrhagia',
'cal_ps_d_023'=	'Myelodysplastic_syndrome',
'cal_ps_d_024'=	'Non-haematological_malignancy',
'cal_ps_d_025'=	'Oral_cavity_bleeding',
'cal_ps_d_026'=	'Other_cardiac_bleeding',
'cal_ps_d_027'=	'Other_eye_bleeding',
'cal_ps_d_028'=	'Other_gynaecological_bleeding',
'cal_ps_d_029'=	'Other_haematological_malignancy',
'cal_ps_d_030'=	'Other_haematuria',
'cal_ps_d_031'=	'Other_obstetric_bleeding',
'cal_ps_d_032'=	'Other_respiratory_system_bleeding',
'cal_ps_d_033'=	'Other_traumatic_bleeding',
'cal_ps_d_034'=	'PAD',
'cal_ps_d_035'=	'Postpartum_haemorrhage',
'cal_ps_d_036'=	'primary_thrombocytopenia',
'cal_ps_d_037'=	'Qualitative_platelet_defects',
'cal_ps_d_038'=	'Renal_disease',
'cal_ps_d_039'=	'secondary_thrombocytopenia',
'cal_ps_d_040'=	'Stroke_NOS',
'cal_ps_d_041'=	'T2D',
'cal_ps_d_042'=	'thrombocytopenia_unspecified',
'cal_ps_d_043'=	'thyroid',
'cal_ps_d_044'=	'TIA',
'cal_ps_d_045'=	'von_Willebrand_disease',
'cal_ps_d_046'=	'PID_allergic_purpura',
'cal_ps_d_047'=	'PID_alopecia',
'cal_ps_d_048'=	'PID_ankylosing_spondilytis',
'cal_ps_d_049'=	'PID_autoimmune_haemolytic_anemia',
'cal_ps_d_050'=	'PID_autoimmune_hepatatis',
'cal_ps_d_051'=	'PID_autoimmune_polyglandular_failure',
'cal_ps_d_052'=	'PID_blistering_disorder',
'cal_ps_d_053'=	'PID_Coeliac_disease',
'cal_ps_d_054'=	'PID_Idiopathic_thrombocytopenic_purpura',
'cal_ps_d_055'=	'PID_juvenile_RA',
'cal_ps_d_056'=	'PID_lupus_erythematosus',
'cal_ps_d_057'= 'PID_MS',
'cal_ps_d_058'= 'PID_other_connective_tissue',
'cal_ps_d_059'=	'PID_psoriasis',
'cal_ps_d_060'=	'PID_RA',
'cal_ps_d_061'=	'PID_reactive_arthropaties',
'cal_ps_d_062'=	'PID_Sarcoidosis',
'cal_ps_d_063'=	'PID_SLE',
'cal_ps_d_064'= 'PID_systemic_sclerosis',
'cal_ps_d_065'=	'PID_T1D',
'cal_ps_d_066'=	'PID_thyroiditis',
'cal_ps_d_067'=	'PID_thyrotoxicosis',
'cal_ps_d_068'=	'PID_vasculitis',
'cal_ps_d_069'=	'PID_vitiligo',
'cal_ps_d_070'=	'PID_IBD',
'cal_ps_d_071'=	'PID_chron',
'cal_ps_d_072'=	'PID_ulcerative_colitis')
#####################################
# groups

miscellaneous = c(
'Acute_lymphoblastic_leukaemia'
'Acute_myeloid_leukaemia'
'AF'
'Aplastic_anaemia',
'CALIBER_DVT'
'CALIBER_PE'
'Stroke_NOS'
'T2D'
'IHD'
'ISSTR'
'Myelodysplastic_syndrome'
'Non-haematological_malignancy'
'thyroid'
'TIA'
'PAD'

'Cataract'
'Hearing_loss'
'primary_thrombocytopenia'
'Qualitative_platelet_defects'
'Renal_disease'
'secondary_thrombocytopenia',
'thrombocytopenia_unspecified'

bleeding=c(
'Bleeding_after_tooth_extraction_or_other_procedure',
'Central_nervous_system_bleeding',
'Cutaneous_bruising_or_bleeding',
'Epistaxis',
'Gastrointestinal_bleeding',
'Haemangioma',
'Haemarthrosis',
'Hereditary_deficiency_of_other_clotting_factors',
'Hereditary_factor_IX_deficiency',
'Hereditary_factor_VIII_deficiency',
'Hereditary_factor_XI_deficiency',
'Menorrhagia',
'Oral_cavity_bleeding',
'Other_cardiac_bleeding',
'Other_eye_bleeding',
'Other_gynaecological_bleeding',
'Other_haematological_malignancy',
'Other_haematuria',
'Other_obstetric_bleeding',
'Other_respiratory_system_bleeding',
'Other_traumatic_bleeding',
'Postpartum_haemorrhage',
'von_Willebrand_disease')

pid = c(
'PID_allergic_purpura',
'PID_alopecia',
'PID_ankylosing_spondilytis',
'PID_autoimmune_haemolytic_anemia',
'PID_autoimmune_hepatatis',
'PID_autoimmune_polyglandular_failure',
'PID_blistering_disorder',
'PID_Coeliac_disease',
'PID_Idiopathic_thrombocytopenic_purpura',
'PID_juvenile_RA',
'PID_lupus_erythematosus',
'PID_MS',
'PID_other_connective_tissue',
'PID_psoriasis',
'PID_RA',
'PID_reactive_arthropaties',
'PID_Sarcoidosis',
'PID_SLE',
'PID_systemic_sclerosis',
'PID_T1D',
'PID_thyroiditis',
'PID_thyrotoxicosis',
'PID_vasculitis',
'PID_vitiligo',
'PID_IBD',
'PID_chron',
'PID_ulcerative_colitis')
#############################################

colnames(df_phenotypes) = str_replace_all(colnames(df_phenotypes), v1)

gene_list = unique(na.omit(df_phenotypes$GENE))

lapply(X = gene_list,
       function(x){
         df_tmp = df_phenotypes %>% filter( is.na(GENE) == TRUE | GENE == x );
         df_tmp[which(is.na(df_tmp[,"GENE"])),"GENE"] = "no_variants"
         lapply( colnames(df_tmp)[5:6], 
               function(y){ 
                 if (sum(is.na(df_tmp[,y])) < length(df_tmp[,y])) {
                   pm = parameters::parameters(glm( as.factor(df_tmp[,eval(y)]) ~ as.factor(df_tmp$GENE) + df_tmp$sex, family = 'binomial'))
                   #return(paste(x, c(pm$Coefficient, pm$p) ))         
                   y
                 }
               }
            )
          }
         )

