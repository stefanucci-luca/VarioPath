#!/usr/bin/Rscript
# In R version 3.6
# Load the libraries
library(ukbtools)
library(dplyr)
library(stringr)

#Set working directory
setwd("/home/ls760/UKBb/info_200K")
dir.create("hpo_info")

# create info df
ukbdf = ukb_df('ukb44092', n_threads= 'max')

# Subset for the 200K release
ukb200k = read.csv("sample_list_200k.txt", header = F)		# read the list of samples that belong to the 200K subset
ukbdf200k = ukbdf %>% 						# use the participants id to subset the entire UK BioBank cohort
		filter(eid %in% ukb200k$V1)
message("The subset has ", dim(ukbdf200k)[1], " rows and ", dim(ukbdf200k)[2], " columns in the datasheet") # print to stdout the info on the subset df

# Extract columns with ICD information for icd9, icd10:
for (icd_version in 9:10){
	ukbdf_icd = ukbdf200k %>%  
	            dplyr::select(matches(paste("^diagnoses.*icd", icd_version, "_f4127*", sep = ""))) # extract the icd columns
			# NOTE: this script consider only the ICD codes used in diagnosis,
	    		# either as first or secondary diagnosis (see fields 41270, 41202 and 41204 on UKBB research website)
			# there are also other ICD codes used to describe cause of death, type of cancer, date of disease onset (which is very interesting)
	message("there are ", dim(ukbdf_icd)[2], " ICD columns version ", icd_version, " in the datasheet") # print the number of columns with ICD values
	# message("First ten entries are:\n", head(colnames(ukbdf_icd)))

	icd_vector = as.vector(t(ukbdf_icd)) # Transform the df to a vector of ICD
	icd_vector <-icd_vector[!is.na(icd_vector)] # Remove NA from the vector above
	message("There are ", length(icd_vector), " ICD codes recorded in version ", icd_version)

	icd_frequency = as.data.frame(table(icd_vector)) 				# Get the ICD code frequencies 
	icd_frequency = icd_frequency[order(icd_frequency$Freq, decreasing = T),]	# Sort the frequencies df in a decreasing order (i.e. most common first) 
	
# Convert ICD codes from UKB ones [UKB reconvert ICD codes, most of the time it is just removing the dot] to original ICDs and add code description
	icd_frequency$icd_meaning =  apply(						# append result of the command to a new df column
					   as.data.frame(icd_frequency[,"icd_vector"]), # select columns with ICD and transform it to a temporary df
					   1,						# apply function by row
					   function(x)
					   ukb_icd_code_meaning(icd.code=eval(x), icd.version=icd_version)[[2]]		#extract code meaning (it has real ICD code and its meaning)		
					   )
# Extract ICD code from the meaning column (first entry of the str_split vector)
	icd_frequency$code = as.character(						# Save to "code" column of the df 
					  lapply(					# apply function to column 'icd_meaning'
						 strsplit(				# split the column on the space (i.e. ' ') character)
							  as.character(icd_frequency$icd_meaning), 		# convert to character to be compatible with the str_split function
							  split=" ", 
						          fixed=2), 			# Ask to return just 2 vectors. The first contains the ICD code, the second contains the description
						 '[', 					# extract from the list
						 1)					# extract just the first entry (I.e. ICD code)
					 )
# Append meaning to another column (second entry of the str_split vector)
	icd_frequency$meaning = unlist(							# Save to "meaning" column of the df 
				       lapply(						# apply function to column 'icd_meaning'. The nested apply is used to extract and convert to string the results
					      lapply(					# of the strsplit function. They would normally be vectors, but the second lapply convert to string 
						     strsplit(
							      as.character(icd_frequency$icd_meaning), 
							      split=" ", 		# split the column on the space (i.e. ' ') character)
							      fixed=2), 		# Return just 2 vectors.
						     '[', 				# extract from the list
						     -1), 				# extract everything but the first entry (all the meanings)
					     toString)					# toString function convert the vector of the first lapply result to a string
					)
	
# Output the table in a tsv file	
	write.table(x= icd_frequency[,c("code", "meaning", "Freq")],
	            file= paste("hpo_info/frequency_ICD", icd_version,"in_the_200k_ukb_cohort.tsv", sep="_"),
	            append=F,
	            quote=F,
	            sep="\t",
		    row.names=F,
	 	    col.names=c(paste("ICD-", icd_version, sep=""), 
				paste("ICD-", icd_version, "_meaning", sep=""),
				"Frequency_in_200K_cohort"
				)
		    )  
}





