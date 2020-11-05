
# Import libraries
library(Gviz)
library(rtracklayer)
library(trackViewer)
library(readr)
library(dplyr)
library(stringr)
require(reshape)
require(GenomicRanges)
require(biomaRt)

# Import TG gene list
tg_genes <- readxl::read_xlsx("~/Desktop/TG_genes_by_diseases.xlsx")
tg_genes_list <- tg_genes$Gene_Symbol_HGNC

# Karyn variants pathogenyc and likely pathogenic
patho_likely_patho_df <- read_tsv("/Volumes/GoogleDrive/My Drive/Phd/VarioPath/variants/ALL_variants_P_LP_someVUS_original_karyn.txt",
                                  guess_max = 290000
                                  )

# Try to parse INFO columns
patho_likely_patho_df = transform(patho_likely_patho_df, INFO = colsplit(INFO, split = ":|;", names = "INFO"))
patho_likely_patho_df = transform(patho_likely_patho_df, INFO_tag = colsplit(INFO_tag, split = ":|;", names = "INFO_tag"))
# assign genes to variants
patho_likely_patho_df$gene = unlist(
                             lapply( 
                             apply(patho_likely_patho_df, 1, function(x) 
                             unique(gsub(".*=","",grep("GENE", x, value = TRUE)))),
                             '[', 1)
                             )

a <- table(patho_likely_patho_df$gene)

# Subset the spreadsheet to a smaller one, just with the information needed
p_lp_df = tibble( "chrom" = patho_likely_patho_df$Chr,
                      "start" = patho_likely_patho_df$Start, 
                      "symbol" = patho_likely_patho_df$gene,
                      "ref" = patho_likely_patho_df$REF,
                      "alt" = patho_likely_patho_df$ALT
                    ) %>% filter(!is.na(as.numeric(start)))


############################################################################################################################
############################################################################################################################
# # Load humand GRCh38 database info from ensembl
# # Ensembl Vesion 86
ensembl_data <- useEnsembl(biomart = "ENSEMBL_MART_ENSEMBL", 
                           version = "101", # this command stick the research to version 101 of ENSEMBL. GRCh38.p13
                           dataset = "hsapiens_gene_ensembl",)
attributes_tx = c( "chromosome_name", "exon_chrom_start", "exon_chrom_end", "strand", 
                   "ensembl_transcript_id", 
                   "ensembl_exon_id"
) 
VarioPath_df <- read.csv("/Volumes/GoogleDrive/My Drive/Phd/VarioPath/genes/disease_gene_list_20201012.csv",
                         header = TRUE, sep = ",", stringsAsFactors = FALSE)
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
} else if (ratio_gene_occurrence == 1) {
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

# Merge Transcript dataframe and gene name one
transcript_df <- merge(Tx_info,
                       Tx_info_symbol,
                       by= "ensembl_transcript_id")



##########################################################################################################################################
##########################################################################################################################################



for (gene in tg_genes_list) {
  message("producing data for ", gene)
  # variants
  gene_variant_df <- p_lp_df[which(p_lp_df$symbol == gene),]
  gene_variant_gr <- GRanges(gene_variant_df$chrom, 
                        IRanges(start = as.numeric(gene_variant_df$start), 
                                width = 1))
  # gene space
  gene_transcript_df <- transcript_df[which(transcript_df$hgnc_symbol == gene),]
  gene_transcript_gr <- GRanges(gene_transcript_df$chromosome_name, 
                     IRanges(start = gene_transcript_df$exon_chrom_start, 
                             end = gene_transcript_df$exon_chrom_end))
  # Variants founded in resources
  
  resource.gr <- GRanges(gene_transcript_df$chromosome_name, 
                      IRanges(start = gene_transcript_df$exon_chrom_start, 
                              end = gene_transcript_df$exon_chrom_end))
  resource.gr$SNPsideID <- "bottom"

  # Producing Plot
  # visual on P and LP variants
  gene_variant_gr$color <- rep(c("blue","Orange"), length.out = length(gene_variant_gr))
                # random 10 color => sample.int(10, length(gene_variant_gr), replace=TRUE) # Color according to Pathogenic or likely pathogenic
  gene_variant_gr$score <- sample.int(25, length(gene_variant_gr), replace = TRUE) # Height of the SNP
  gene_variant_gr$border <- sample(c("gray20"), length(gene_variant_gr), replace=TRUE)
  #gene_variant_gr$label <- as.character(1:length(gene_variant_gr))
  gene_variant_gr$alpha <- 1
  gene_variant_gr$SNPsideID <- "top"
  # Visual on variants founded in resources
  idx <- sample.int(length(gene_variant_gr), length(gene_variant_gr)/2)
  gene_variant_gr$color[idx] <- "grey"
  gene_variant_gr$SNPsideID[idx] <- "bottom"
  gene_variant_gr$border[idx] <- "black"
  gene_variant_gr$alpha[idx] <- 0.3
  gene_variant_gr$score[idx] <- 1
  # Visual on genes bodies/transcript
  gene_transcript_gr$fill <- c("#51C6E6")
  gene_transcript_gr$height <- c(0.03)
  # Plot
  lolliplot(gene_variant_gr, 
            gene_transcript_gr, 
            type="circle", 
            cex=.8,
            xaxis = FALSE, 
            yaxis=FALSE,
            ylab="")
  grid.text(paste("variants for", gene), x=.5, y=.5, just="top", 
            gp=gpar(cex=1.5, fontface="bold"))
}



# create the GRanges for the P_LP variants
p_lp_df.gr <- GRanges(p_lp_df$chrom, 
                     IRanges(start = p_lp_df$start, 
                             width = 1))
p_lp_df.gr$gene_symbol = p_lp_df$symbol



# Function to plot
plot_lollipop_variopath <- function(gene_length, p_lp, all_var) {
  
  
}







## shape must be "circle", "square", "diamond", "triangle_point_up", or "triangle_point_down"
available.shapes <- c("circle", "square", "diamond", 
                      "triangle_point_up", "triangle_point_down")
sample.gr$shape <- sample(available.shapes, size = length(sample.gr), replace = TRUE)
sample.gr$legend <- paste0("legend", as.numeric(factor(sample.gr$shape)))
lolliplot(sample.gr, features, type="circle", legend = "legend", xaxis = FALSE)

sample.gr$SNPsideID <- sample(c("top", "bottom"), 
                              length(sample.gr),
                              replace=TRUE)


