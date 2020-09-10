#!/usr/bin/env Rscript

require(GenomicRanges)
require(biomaRt)
require(dplyr)
require(stringr)

# # Load humand GRCh38 database info from ensembl
# # Ensembl Vesion 86
ensembl_data <- useEnsembl(biomart = "ENSEMBL_MART_ENSEMBL", 
                           version = "101", # this command stick the research to version 101 of ENSEMBL. GRCh38.p13
                           dataset = "hsapiens_gene_ensembl",)

##################################
# Mart Setup
##################################

# The command in this section are helpful to set up the research ensembl mart.
# The command are in comments, to get them to work remove the hash. 

# # Show available
# listMarts()

# # Show available dataset
# listDatasets(ensembl_data)

# # Show available filters
# e.g. ens_hs_transcript or ensembl_transcript_id to use the trascripts as filter
#listFilters(ensembl_data) 

# # show attributes (i.e. the values we want to retrieve)
# # e.g. chromosome_name start_position end_position strand ensembl_exon_id
# listAttributes(ensembl_data) 
# searchAttributes(ensembl_data, "exon") # Select 

#################################
# Define the research attributes
#################################

# these attrivutes are the ones returned in the db after the filtering. 
# NB biomaRt can't research gene and trasncript information in the same research.
# For this reason gene name are attached below after a second research.

# Transcript attributes to have in the final dataframe

attributes_tx = c( "chromosome_name", "exon_chrom_start", "exon_chrom_end", "strand", 
                   "ensembl_transcript_id", 
                   "ensembl_exon_id"
                   ) 

# ########################
# ### Import Transcritps file
# ########################
 
# Read csv from google drive - shared folder
# NB. the file is readin the csv copy I made on my local machine

VarioPath_df <- read.csv("/Users/luca/Desktop/VarioPath/disease_gene_list_20200601_transcripts_local_copy_made_on_20200909.csv",
                         header = TRUE, sep = ",", stringsAsFactors = FALSE)

# return a few info to check if the file is ok
message( "number of VarioPath genes = ", dim(VarioPath_df)[1])
print(head(VarioPath_df))

#############################
#   Quick QCs on the data
#############################

# Check how many entries have double occurrences

# Double genes
gene_occurrence <- as.data.frame(sort(table(VarioPath_df$Approved.symbol), decreasing = TRUE))
ratio_gene_occurrence <- sum(gene_occurrence$Freq) / length(gene_occurrence$Freq)

if (ratio_gene_occurrence > 1) {
        cat("Some genes occurre more than once")
        print(head(gene_occurrence))
} else if (ration_gene_occurrence == 1) {
        cat("All genes occurre just once")
}

# Double trasnscript IDs
tx_occurrence <- as.data.frame(sort(table(VarioPath_df$LRG...MANE...UP.ENST.no.versioning), decreasing = TRUE))
ratio_txID_occurrence <- sum(tx_occurrence$Freq) / length(tx_occurrence$Freq)

if (ratio_txID_occurrence > 1) {
        print("Some genes occurre more than once:")
        print(head(tx_occurrence))
} else if (ratio_txID_occurrence == 1) {
        cat("All genes occurre just once")
}

# Remove lines with wrong ENSEMBL transcript IDs (i.e. NA and NCBIs)

vect_missing_values <- which(! grepl(x = VarioPath_df[,dim(VarioPath_df)[2]], patter= "ENST"))

if (length(vect_missing_values) > 0 ) {

        # Return message stating which lines have a problem. 
        message("There are lines that are missing a valid ENSEMBL transcript ID in the last column. See stdout fro more info")
        cat("The following lines are missing a valid ENSEMBL transcript ID in the last column:")
        print(VarioPath_df[vect_missing_values,])
        
        # Print a mesage that says which lines will be removed
        message("These lines will be removed")
        cat("These lines ("); cat(length(vect_missing_values)); cat(") will be removed:")
        cat(vect_missing_values)
        
        # Remove the missing values lines
        VarioPath_df_clean <- VarioPath_df[ - vect_missing_values,]

        # extract list trascript no versioning
        VarioPath_trx <- VarioPath_df_clean[,dim(VarioPath_df_clean)[2]]
} else if (length(vect_missing_values) == 0 ) {
        message("All entries have a ENSEMBL transcript ID")
        VarioPath_trx <- VarioPath_df[,dim(VarioPath_df)[2]]
        }

#######################################
# Query the Mart ensembl database
#######################################

Tx_info <- getBM(attributes = attributes_tx,
                 filters = "ensembl_transcript_id",
                 values = VarioPath_trx,  
                 uniqueRows = TRUE,
                 mart = ensembl_data)

# Get the gene symbols
# Biomart webserver/database set up can't do gene-based queries and transcript-based ones
# For this reason gene names (i.e. symbols) have to be retrieved in a second query
Tx_info_symbol <- getBM(attributes = c("ensembl_transcript_id", "hgnc_symbol"),
                        filters = c("ensembl_transcript_id"),
                        values = VarioPath_trx, 
                        uniqueRows = TRUE,
                        mart = ensembl_data)

# Merge Transcript dataframe and gene name one
transcript_df <- merge(Tx_info,
                      Tx_info_symbol,
                      by= "ensembl_transcript_id")


#############################
#   Quick QCs
#############################

# are all the genes retrieved?

if (length(unique(transcript_df$hgnc_symbol)) - length(unique(VarioPath_df_clean$Approved.symbol)) < 0) {
        message("the query is missing some data, try to recover them by using uniprot trasncript")
        missing_transcript <- VarioPath_df_clean[which(! VarioPath_df_clean[,dim(VarioPath_df_clean)[2]] %in% 
                                                               transcript_df$ensembl_transcript_id ),]
        missing_symbol <- VarioPath_df_clean[which(! VarioPath_df_clean$Approved.symbol %in% 
                                                           transcript_df$hgnc_symbol ),]
        cat("the query is missing these genes:")
        print(missing_symbol)
        cat("the query is missing these transcripts:")
        print(missing_transcript)
} else if (length(unique(transcript_df$hgnc_symbol)) - length(unique(VarioPath_df_clean$Approved.symbol)) == 0) {
        print("the query has found all the genes.")     
        }



