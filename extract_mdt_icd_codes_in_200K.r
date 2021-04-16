#!/usr/bin/R

# using R version 3.6
library(ukbtools)
library(dplyr)
library(stringr)

group_icd=c("R040", "R233", "R041", "K920", "K921", "K922", "K625", "K25", "K26", 
            "K27", "K28", "K290", "K64" , "I84", "K766", "K552", "R31", "T81", "T792", 
            "T81", "N920", "N921", "N922", "N924", "O72", "O720", "O71", "O722", 
            "M2500", "M2501", "M2502", 'M2503', "M2504", "M2505", "M2506", "M2507", 
            "M2508", "M2509","I610", "I611", "I612", "I613", "I614", "I615", "I616", 
            "I617", "I618", "I619", 'I620', "I621", "I629", "N923", "R58", "D699") 

gourp_of_columns_to_look_at=c("diagnoses_main_icd10",
                              "diagnoses_secondary_icd10",
                              "diagnoses_main_icd9",
                              "diagnoses_secondary_icd9")

# function to extract the icd codes from every line (i.e. partecipants)
# takes 3 arguments: 
#       dataframe, 
#       columns to use to subset the dataframe. they can be ICD9 or 10, or main and secondary diagnosis.
#       icd_codes. which are the selected codes relevant for the MDTs
# the output can be appended as new columns to the UK Biobank df
extract_individual_icd_codes = function(data_frame, column_subset, icd_codes) {
                                    unlist(    apply(                                                              
                                                                data_frame[,column_subset],                 # subset the df to the relevant columns (this allows to look in a different way to primary and secondary ICD codes)
                                                                1,                                          # apply the function by row - per individual
                                                                function(y)                                 # to every partecipant
                                                                unlist(
                                                                    lapply(y,                               # look across all the selected ICD columns
                                                                        function(x) 
                                                                        if ( x %in% icd_codes ) {           # if the ICD code belong to the subset then return it
                                                                            return(x)
                                                                        }
                                                                    )     
                                                                )
                                                            ) )
}
#
select_icd_record_columns = function(data_frame, icd_column) {
                                    colnames(data_frame %>% 
                                                    dplyr::select(matches(eval(icd_column)))
                                    )
}
#
ukbdf = ukb_df('ukb44092', n_threads= 'max')
# Subset for the 200K release
ukb200k = read.csv("sample_list_200k.txt", header = F)		# read the list of samples that belong to the 200K subset

ukbdf200k = ukbdf %>% 						# use the participants id to subset the entire UK BioBank cohort
		filter(eid %in% ukb200k$V1)

for (i in gourp_of_columns_to_look_at){
    col_sub=select_icd_record_columns(ukbdf200k, i)
    ukbdf200k[,i] = ""
    ukbdf200k[,i] = extract_individual_icd_codes(ukbdf200k, col_sub, group_icd)
}
#