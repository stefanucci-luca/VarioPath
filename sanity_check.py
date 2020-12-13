#!/usr/bin/env python3
import os
import shutil
import tempfile
from time import sleep
import pysam
import pandas as pd
import numpy as np
import matplotlib as plt
import re
import seaborn as sns

# write in the path to the VCF or VCF like files to read
vcf_to_read=['/home/ls760/hgmd_pro_2019.4_hg38.vcf',
              '/home/ls760/clinvar_20200615.vcf'] 
# write in the path to the file(s) that resemble a bed file but can't be parsed from any reader
bed_to_read=['/home/ls760/curatedNBR_CSVD-TG-PAH-EAHAD-GRID-GPS_P-LP-VUS_GRCh38_v2_AssemConvert_sorted.bed']


# Creates the list wjere I want to store the infromation from all the VCFs
CHROM = []
POS = []
REF = []
ALT = []
sourx = []                                        # create a list with the sources of the variants.
INFO=[]                                           # INFO column should contain all the information that are not chrom,pos,ref,alt
GENE=[]
PATHO=[]

actual_sir=os.getcwd()
with tempfile.TemporaryDirectory() as tmpdirname: # Create a tmp dir that store the files and indexes
  dir=tmpdirname                                  # it will be removed as soon as the for loop finishes     
  os.chdir(dir)
  for file in vcf_to_read:
    shutil.copy(file, dir)                       # Copy temprarily the file in a tmp working dir 
    filename = os.path.basename(file)
    print(filename)
    tbx = pysam.tabix_index(filename, preset='vcf')   # Creates the indexes for the files
    vcf = pysam.VariantFile(tbx)                  # read the files in
    vcf.close                                     # close the connections to the files
    for rec in vcf.fetch():
      info=[]
      CHROM.append(rec.chrom) 
      POS.append(rec.pos) 
      REF.append(rec.ref)
      if  isinstance(rec.alts, tuple):
        ALT.append(rec.alts[0])
      else:
        ALT.append(rec.alts)
      sourx.append(filename.split("_")[0])
      if filename.split("_")[0] == 'hgmd':
        PATHO.append(','.join(rec.info['CLASS']))
        GENE.append(rec.info['GENE'])
      else:
        PATHO.append(','.join(rec.info['CLNSIG']))
        if 'GENEINFO' in rec.info:
          GENE.append(rec.info['GENEINFO'].split(':')[0])
        else:
          GENE.append("gene information not present")
      for i in range(len(list(rec.info))):
        #info.append(rec.info[list(rec.info)[i]])
        strings = [list(rec.info)[i], rec.info[list(rec.info)[i]]] # it works on a single record. It creates a list of INFO_ID and INFO value.
        info.append('='.join(map(str, strings)))                   # concatenate values in a string e.g 'RS"=rs123 
      INFO.append(','.join(map(str, info)))                        # concatenate the string in an INFO list for all teh entries 
  shutil.rmtree(dir, ignore_errors=False)
  os.chdir(actual_sir)    

# It is not possible to parse this file with existing parser. I created this simple one that extract the info I want.
# Everything else goes in the INFO column
for file in bed_to_read:
  with open(file) as bed_file:
    for line in bed_file:
      CHROM.append(line.split("\t")[0])
      POS.append(line.split("\t")[1])
      splitted = line.split("_")
      gene_sor = re.search('(_(.*)-NM|_(.*)-ENST|_(.*)-NR|_(.*)-c.)', line)
      GENE.append(gene_sor.group(1).split('_')[-1].split('-')[0])
      POS_37 = re.findall(r'\d{5,20}', line)[2]
      index_pos = [splitted.index(i) for i in splitted if POS_37 in i]
      REF.append(splitted[index_pos[0]+1])
      ALT.append(splitted[index_pos[0]+2])
      INFO.append(line.split("\t")[4])
      sourx.append(os.path.basename(file).split("_")[0])


# From the lists create a dataframe
df = pd.DataFrame({'CHROM':CHROM, 
                   'POS':POS, 
                   'REF':REF, 
                   'ALT':ALT,
                   'GENE':GENE,
                   'source':sourx,
                   'INFO':INFO})


# Remove duplicates variant entries and
# Pivot the dataframe to have source of information on the same line. 
df2 = df.groupby(['CHROM', 'POS', 'REF', 'ALT','GENE']).agg({'source' : lambda x: ','.join(x), 
                                                      'INFO': lambda x: ';' .join(x)})
df2 = df2.sort_values(by=['CHROM', 'POS'])

df2.to_csv('/home/ls760/aggregated_karyn_resources.csv')

# Number of variants in the 
n_variatn_df = pd.read_csv("~/Desktop/VarioPath/number_variants_in_MDT_genes_clean.tsv",  delimiter='\t', skiprows=1, dtype={'n_karyn_list':np.int32, 'n_variant_in_UKB':np.int32}) 
n_variatn_df.head
n_variatn_df.astype({'n_karyn_list': 'int32', 'n_variant_in_UKB': 'int32' })

sns.relplot(x='n_karyn_list', y='n_variant_in_UKB', data = n_variatn_df, )
plt.show()

# list pre pathogenicity filter

sns.relplot(x='n_karyn_list', y='n_variant_in_UKB', data = df2 )
plt.show()
