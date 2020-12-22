#!/usr/bin/Rscript
library(readr)
library(dplyr)
require(reshape)
library(stringr)
library(jsonlite)
require(biomaRt)
library(purrr)
library(ggplot2)

project_dir <- "~/Desktop/VarioPath/"
setwd(project_dir)

# get pLI info from gnoma paper PMID32461654 ( Nat. 2020 )
df_gnomad = read_delim('gnomad_paper_PMID32461654_supplementary_dataset_11_full_constraint_metrics.tsv', delim = '\t')
# create a suisbet of the entire dataframe
pli_df=df_gnomad[,c("transcript", "pLI")]
# renoame transcript columns for consistency with the ensembl output
colnames(pli_df)[1] = "ensembl_transcript_id"

# Get transcript information from ENSEMBL

# Load humand GRCh38 database info from ensembl
# Ensembl Vesion 86
ensembl_data <- useEnsembl(biomart = "ENSEMBL_MART_ENSEMBL", 
                           version = "101", # this command stick the research to version 101 of ENSEMBL. GRCh38.p13
                           dataset = "hsapiens_gene_ensembl",)

# Transcript attributes to have in the final dataframe
attributes_tx = c( "chromosome_name", "exon_chrom_start", "exon_chrom_end", "strand", 
                   "ensembl_transcript_id", 
                   "ensembl_exon_id"
) 
# Query the Mart ensembl database
Tx_info <- getBM(attributes = attributes_tx,
                 filters = "ensembl_transcript_id",
                 values = df_gnomad$transcript,  
                 uniqueRows = TRUE,
                 mart = ensembl_data)

# calculate exon lengths
Tx_info$exon_length = Tx_info$exon_chrom_end - Tx_info$exon_chrom_start
# Calculate ORF length
orf_length = aggregate(Tx_info$exon_length, 
                       by=list(ensembl_transcript_id=Tx_info$ensembl_transcript_id), 
                       FUN=sum)
colnames(orf_length)[2] = c("length_ORF")

# Get the gene symbols
# Biomart webserver/database set up can't do gene-based queries and transcript-based ones
# For this reason gene names (i.e. symbols) have to be retrieved in a second query
Tx_info_symbol <- getBM(attributes = c("ensembl_transcript_id", "hgnc_symbol"),
                        filters = c("ensembl_transcript_id"),
                        values = df_gnomad$transcript,   
                        uniqueRows = TRUE,
                        mart = ensembl_data)

proteins_info = getSequence(id=df_gnomad$transcript,
                            type="ensembl_transcript_id",
                            seqType="peptide", 
                            mart=ensembl_data)

proteins_info$length_protein = str_length(proteins_info$peptide)
# merge pLI info to ORF length and protein length
complete_df = list(orf_length, 
                   proteins_info, 
                   pli_df) %>% 
              reduce(inner_join, 
                     by = "ensembl_transcript_id")
# group list according to ORF lengths
complete_df$group_orf <- cut(complete_df$length_ORF, 
                         breaks= c( seq(quantile(complete_df$length_ORF)[1],
                                        quantile(complete_df$length_ORF)[3],
                                        100), 
                                    seq(quantile(complete_df$length_ORF)[3],
                                        quantile(complete_df$length_ORF)[5],
                                        500) )
                         )

complete_df$group_protein <- cut(complete_df$length_protein, 
                                 breaks= seq(min(complete_df$length_protein),
                                             max(complete_df$length_protein),
                                             500)
                          )
# remove eventual NA
complete_df = na.omit(complete_df)
# plot correlation ORF to pLI
ggplot(complete_df, aes( complete_df$pLI, y = complete_df$length_ORF)) + 
  geom_point() 
# trascriptome pLI
ggplot(complete_df, aes( x= pLI, ..count..)) + 
  geom_density() 
# trascriptome pLI grouped by ORF length
ggplot(complete_df, aes(pLI, ..count.., color=group_orf )) + 
  geom_density() + 
  facet_wrap('group_orf', scales = 'free_y')
# trascriptome pLI grouped by protein length
ggplot(complete_df, aes(pLI, ..count.., color=group_protein )) + 
  geom_density() + 
  facet_wrap('group_protein', scales = 'free_y')
