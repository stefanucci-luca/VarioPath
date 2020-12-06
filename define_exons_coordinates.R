#!/usr/bin/env Rscript

require(GenomicRanges)
require(biomaRt)
require(dplyr)
require(plotly)
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
# Define the Mart research attributes
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

VarioPath_df <- read.csv("/Users/luca/Desktop/VarioPath/disease_gene_list_20200915.csv",
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
tx_occurrence <- as.data.frame(sort(table(VarioPath_df$chosen.transcript), decreasing = TRUE))
ratio_txID_occurrence <- sum(tx_occurrence$Freq) / length(tx_occurrence$Freq)

if (ratio_txID_occurrence > 1) {
        print("Some genes occurre more than once:")
        print(head(tx_occurrence))
        wrong_ensID <- tx_occurrence[which(tx_occurrence$Freq > 1),]
        wrong_ensID_original_table <- VarioPath_df[which(VarioPath_df$chosen.transcript %in% wrong_ensID[,1]),]
        print(wrong_ensID_original_table)
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


proteins_info = getSequence(id=Tx_info_symbol$ensembl_transcript_id,
                  type="ensembl_transcript_id",
                  seqType="peptide", 
                  mart=ensembl_data)
proteins_info$length_protein = str_length(proteins_info$peptide)

# Merge Transcript dataframe and gene name one
transcript_df <- merge(Tx_info,
                      Tx_info_symbol,
                      by= "ensembl_transcript_id")

transcript_df_protein = merge(transcript_df, proteins_info, by="ensembl_transcript_id")

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

###################################
# Exons length
###################################

transcript_df$exon_length <- transcript_df$exon_chrom_end - transcript_df$exon_chrom_start
orf_length = aggregate(transcript_df$exon_length, by=list(transc=transcript_df$ensembl_transcript_id), FUN=sum)

transcript_df2 <- merge(transcript_df,
                       orf_length,
                       by.x= "ensembl_transcript_id",
                       by.y= "transc")
colnames(transcript_df2)[9] <- "orf_length"

transcript_df2[which(transcript_df2$hgnc_symbol == "VWF"),]

write.table(transcript_df2, "~/Desktop/VarioPath/Exon_information.tsv")

####################################
# plot gene size
####################################

library(ggplot2)
library(dplyr)
df_geno_plot <- arrange(unique(transcript_df2[,c("hgnc_symbol","orf_length")]), hgnc_symbol) 
ggplot(df_geno_plot, aes(x = orf_length, y = ..count..)) + geom_histogram( bins = 1000) + xlim(0,100000)

library(readr)
library(dplyr)
require(reshape)
library(stringr)
library(jsonlite)

# General settings
project_dir <- "~/Desktop/VarioPath/"
setwd(project_dir)

# Karyn variants pathogenyc and likely pathogenic
patho_likely_patho_df <- read_tsv("/Volumes/GoogleDrive/My Drive/Phd/VarioPath/variants/ALL_variants_P_LP_someVUS_original_karyn.txt",
                                  guess_max = 290000
)


# Try to parse INFO columns
patho_likely_patho_df = transform(patho_likely_patho_df, INFO = colsplit(INFO, split = ":|;", names = "INFO"))
patho_likely_patho_df = transform(patho_likely_patho_df, INFO_tag = colsplit(INFO_tag, split = ":|;", names = "INFO_tag"))

# extract information from the info and info_tag colum
info_to_extract = c("GENE")

for (info_value in info_to_extract) {
  message("extracting info ", info_value )
  patho_likely_patho_df[,info_value] = unlist(
    lapply( 
      apply(patho_likely_patho_df, 1, function(x) 
        unique(gsub(".*=","",grep(info_value, x, value = TRUE)))),
      '[', 1)
  )
} 


df1 <- as.data.frame(table(patho_likely_patho_df$GENE))

df3 = merge(transcript_df_protein, df1, by.x="hgnc_symbol", by.y="Var1") 
df3 = merge(transcript_df2,  df3, by="hgnc_symbol")

var_per_len_prot = df3 %>% 
                   select(c("hgnc_symbol","orf_length","length_protein","Freq")) %>% 
                   unique()
var_per_len_prot$per_100 <- var_per_len_prot$length_protein / 100
var_per_len_prot$var_norm <- var_per_len_prot$Freq/var_per_len_prot$per_100

# Plot the number of variants per gene
plot1 = ggplot(data = var_per_len_prot, mapping= aes(x = orf_length, y = var_norm, 
                                                     label = hgnc_symbol, 
                                                     text = length_protein)) + 
  geom_point()+
  #geom_label(position = "jitter", label.size = 0.01)
  theme_minimal() + 
  theme()
plotly_build(plot1)

# Showing only the first 50
ggplot(df3[order(df3$Freq, decreasing = T ),][1:50,], 
       aes(x = reorder(hgnc_symbol, Freq), y = Freq)) + 
  geom_col() + 
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 90, size = 7)) 
# Normalise by the length of the the protein
df3$Freq_norm_orf = df3$Freq / (df3$orf_length / 100 )
# show the frequency
ggplot(df3, aes(x = reorder(hgnc_symbol, Freq_norm_orf), y = Freq_norm_orf)) + 
  geom_col() + 
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 90, size = 2))
# Showing only the first 50
ggplot(df3[order(df3$Freq_norm_orf, decreasing = T ),][1:50,], 
       aes(x = reorder(hgnc_symbol, Freq_norm_orf), y = Freq_norm_orf)) + 
  geom_col() + 
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 90, size = 7)) 


# Plot variants information 
# only MNV
patho_likely_patho_df_MNV = patho_likely_patho_df[which(str_length(patho_likely_patho_df$REF) > 1 
                                  & str_length(patho_likely_patho_df$ALT) > 1),]
# remove MNV
patho_likely_patho_df_NO_MNV = patho_likely_patho_df[- which(str_length(patho_likely_patho_df$REF) > 1 
                                                        & str_length(patho_likely_patho_df$ALT) > 1),]
# Large deletion
patho_likely_patho_df_DEL = patho_likely_patho_df[which(patho_likely_patho_df$ALT == "<DEL>"),]
# Remove large deletion
patho_likely_patho_df_NO_MNV_DEL = patho_likely_patho_df_NO_MNV[- which(patho_likely_patho_df_NO_MNV$ALT == "<DEL>"),]
# Select only small INDEL
patho_likely_patho_df_SMALL_INDEL = patho_likely_patho_df_NO_MNV_DEL[ 
  which(
    str_length(patho_likely_patho_df_NO_MNV_DEL$REF) != 1 |
      str_length(patho_likely_patho_df_NO_MNV_DEL$ALT) != 1),]
# Select only SNV
patho_likely_patho_df_ONLY_SNV = patho_likely_patho_df_NO_MNV_DEL[ 
    which(
      str_length(patho_likely_patho_df_NO_MNV_DEL$REF) == 1 &
      str_length(patho_likely_patho_df_NO_MNV_DEL$ALT) == 1),]

dim(patho_likely_patho_df_MNV);
dim(patho_likely_patho_df_DEL);
dim(patho_likely_patho_df_SMALL_INDEL);
dim(patho_likely_patho_df_ONLY_SNV)

######################################
# Consequences
######################################

variant_VEP = read_delim("~/Desktop/VarioPath/annotation_karyn_list_withVEP_only_correct_trasncripts.tab",delim = "\t" )

variant_VEP$aa_change = unlist(lapply(str_split(string = variant_VEP$HGVS_OFFSET,
                 pattern = "\\."), `[`, 3))

variant_VEP$aa_position_changed = gsub("([0-9]+).*$", "\\1", variant_VEP$aa_change )

table_aa_change = as.data.frame(table(paste(variant_VEP$SYMBOL, variant_VEP$aa_position_changed, sep = "_" )))
table_aa_change = table_aa_change[-grep(pattern = "*_NA", x = table_aa_change$Var1),]
table_aa_change$name <- lapply( str_split(table_aa_change$Var1, pattern = "_"), '[', 1)
table_aa_change <- table_aa_change[order(table_aa_change$Freq, decreasing = T),]  
plot_1 = ggplot(table_aa_change[1:1000,], 
       aes(x = reorder(Var1, Freq), y = Freq)) + 
  geom_point() + 
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 90, size = 4)) 

plotly_build(plot_1)

##########################################
# Slide #2
# ##########################################

# Slide #2: table with numbers of SNVs and indels  (so 2 columns, 3 rows)
# o   In total
# o   HIGH impact (LOFs, frameshifts – respectively) – he prefers these to be called LOFs / frameshift as opposed to HIGH as more meaningful
# o   MODERATE impact (missenses, inframe indels - respectively) – he prefers these to be called nsSNPs / inframe indels as opposed to MODERATE as more meaningful to ‘proteinologists’  (I’ve just made that word up)

message('number of variants in the list: ', dim(variant_VEP))
message('number of SNV: ', length(which(variant_VEP$VARIANT_CLASS == "SNV")))
message('number of INDEL / frameshifts: ', length(which(variant_VEP$VARIANT_CLASS != "SNV")))

message('number of SNV and LOFs / frameshifts: ', length(which(variant_VEP$VARIANT_CLASS == "SNV" & variant_VEP$IMPACT == "HIGH")))
message('number of INDEL and LOFs / frameshifts: ', length(which(variant_VEP$VARIANT_CLASS != "SNV" & variant_VEP$IMPACT == "HIGH")))

message('number of SNV and nsSNPs / inframe indels: ', length(which(variant_VEP$VARIANT_CLASS == "SNV" & variant_VEP$IMPACT == "MODERATE")))
message('number of INDEL and nsSNPs / inframe indels: ', length(which(variant_VEP$VARIANT_CLASS != "SNV" & variant_VEP$IMPACT == "MODERATE")))

message('number of SNV and nsSNPs / inframe indels: ', length(which(variant_VEP$VARIANT_CLASS == "SNV" & variant_VEP$IMPACT == "LOW")))

