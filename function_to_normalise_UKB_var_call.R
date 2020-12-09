library("stringr")
df = data.table::fread("karyns_variant_from_200KWES_biallelic_records_only_unrelated_variant_normalised_filtered_reshaped_reid.tab")
df$mapply_dist = mapply(function(x,y) which.min(x==y),strsplit(df$REF,""),
                        strsplit(df$ALT,""))

df <- dplyr::mutate_all(df, as.character)

df2 = df

for (variable in 1:dim(df)[1]) {
  if (df[variable,"mapply_dist"] == 2){
    if (str_length(df2[variable,"REF"]) > str_length(df2[variable,"ALT"]) ) {
      df2[variable,"POS"] <- as.numeric(df2[variable,"POS"]) + 1
      df2[variable,"REF"] <- paste0(na.omit(strsplit(df[variable,"REF"],"")[[1]][df[variable,"mapply_dist"]:length(strsplit(df[variable,"REF"],"")[[1]])]), collapse="")
      if (str_length(df2[variable,"ALT"]) == 1) {
        df2[variable,"ALT"] <- "-"
      } else {
      df2[variable,"ALT"] <- paste0(na.omit(strsplit(df[variable,"ALT"],"")[[1]][1:length(strsplit(df[variable,"ALT"],"")[[1]])]), collapse="")
      }
    } else {
      df2[variable,"POS"] <- as.numeric(df2[variable,"POS"]) + 1
      df2[variable,"ALT"] <- paste0(na.omit(strsplit(df[variable,"ALT"],"")[[1]][df[variable,"mapply_dist"]:length(strsplit(df[variable,"ALT"],"")[[1]])]), collapse="")
      if (str_length(df2[variable,"REF"]) == 1) {
        df2[variable,"REF"] <- "-"
      } else {
        df2[variable,"REF"] <- paste0(na.omit(strsplit(df[variable,"REF"],"")[[1]][1:length(strsplit(df[variable,"REF"],"")[[1]])]), collapse="")
      }
    }
  } else if (df[variable,"mapply_dist"] > 2){
    if (str_length(df2[variable,"REF"]) > str_length(df2[variable,"ALT"]) ) {
      df2[variable,"POS"] <- as.numeric(df2[variable,"POS"]) + 1
      df2[variable,"REF"] <- paste0(na.omit(strsplit(df[variable,"REF"],"")[[1]][2:length(strsplit(df[variable,"REF"],"")[[1]])]), collapse="")
      if (str_length(df2[variable,"ALT"]) == 1) {
        df2[variable,"ALT"] <- "-"
      } else {
        df2[variable,"ALT"] <- paste0(na.omit(strsplit(df[variable,"ALT"],"")[[1]][1:length(strsplit(df[variable,"ALT"],"")[[1]])]), collapse="")
      }
    } else {
      df2[variable,"POS"] <- as.numeric(df2[variable,"POS"]) + 1
      df2[variable,"ALT"] <- paste0(na.omit(strsplit(df[variable,"ALT"],"")[[1]][2:length(strsplit(df[variable,"ALT"],"")[[1]])]), collapse="")
      if (str_length(df2[variable,"REF"]) == 1) {
        df2[variable,"REF"] <- "-"
      } else {
        df2[variable,"REF"] <- paste0(na.omit(strsplit(df[variable,"REF"],"")[[1]][1:length(strsplit(df[variable,"REF"],"")[[1]])]), collapse="")
      }
    }
  }
}

df2$NEW_ID <- paste(df2$CHROM, df2$POS, df2$REF, df2$ALT , sep="_")

write.table(df2, "karyns_variant_from_200KWES_biallelic_records_only_unrelated_variant_normalised_filtered_reshaped_reid_my_norm.tab", sep="\t", quote = F, row.names=F, col.names=T)
