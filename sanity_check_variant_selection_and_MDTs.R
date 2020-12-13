

# Allele frequency calculated in UK BioBank script and my calculation (with plink2)
# Example 3 variants don't correlate well between my maths and UKB's one.

plot(df_clean$AF_ukb_calc, df_clean$EUR_ALT_FREQS)
cor(df_clean$AF_ukb_calc, df_clean$EUR_ALT_FREQS, method = 'pearson')
not_corr_df = df_clean[which(df_clean$AF_ukb_calc > 0.0008 & df_clean$EUR_ALT_FREQS < 0.0002),]

not_corr_df$recalc <- ( (dum_df$EUR_ALT_FREQS * dum_df$EUR_OBS_CT) + 
                          (dum_df$AFR_ALT_FREQS * dum_df$AFR_OBS_CT) + 
                          (dum_df$SAS_ALT_FREQS * dum_df$SAS_OBS_CT) +
                          (dum_df$EAS_ALT_FREQS * dum_df$EAS_OBS_CT) +
                          (dum_df$AMR_ALT_FREQS * dum_df$AMR_OBS_CT) ) / 
  (dum_df$EUR_OBS_CT + dum_df$AFR_OBS_CT + dum_df$SAS_OBS_CT + dum_df$EAS_OBS_CT + dum_df$AMR_OBS_CT)
plot(not_corr_df$AF_ukb_calc, not_corr_df$recalc)


# Dragana's email 
# Here is the plot I was talking about:
# Slide 3 in the presentation linked below.
# Maybe interesting for Luca as well for your sanity checks. 
# https://docs.google.com/presentation/d/1vIOTmltmXSlS8doJ2pvpSwc3fhmbuSDAziY8LsBODKY/edit#slide=id.p [docs.google.com]



#  From Karyn email
# Following up on the discussion this morning, re: the number of variants for the top BPD genes, I would do, for each gene:
#   Number of variants in curatedPub, curatedNBR and ClinVar BEFORE any filtering based on pathogenicity (so in the source files)
# Number of variants in curatedPub, curatedNBR and ClinVar AFTER filtering based on pathogenicity (so in the official variant file)
# Number of vairants in the UKB population after your filtering
# 
# I will include the source files for curatedNBR and ClinVar in    VarioPath > variant   .
# curatedPub cannot be moved there and I’ll tell you where you can access it.