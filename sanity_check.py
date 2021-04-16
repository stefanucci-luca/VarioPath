#!/usr/bin/env python3
import os
import shutil
import tempfile
import pysam
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import re
import seaborn as sns
from random import random

# write in the path to the VCF or VCF like files to read
vcf_to_read=['/home/km369/.random_stuff/2019-4/hgmd_pro_2019.4_hg38.vcf',
              '/home/km369/project-cbrc/VarioPath/clinvar_20200615.vcf'] 
# write in the path to the file(s) that resemble a bed file but can't be parsed from any reader
bed_to_read=['/home/km369/project-cbrc/VarioPath/curatedNBR_CSVD-TG-PAH-EAHAD-GRID-GPS_P-LP-VUS_GRCh38_v2_AssemConvert.bed']


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
      if  isinstance(rec.alts, tuple):          # in hgmd the ALT come as tuple. isinsatnce check if tuple
        ALT.append(rec.alts[0])                 # if tuple take only the first allele. It never has more than one allele
      else:
        ALT.append(rec.alts)                    # if not tuple just add to the list
      sourx.append(filename.split("_")[0])      
      if filename.split("_")[0] == 'hgmd':
        PATHO.append(rec.info['CLASS'])
        GENE.append(rec.info['GENE'])
      else:
        if 'CLNSIG' in rec.info:
          PATHO.append(','.join(rec.info['CLNSIG']))
        else:
          PATHO.append(','.join(rec.info['CLNSIGINCL'][0].split(":")[1]))
        if 'GENEINFO' in rec.info:
          GENE.append(rec.info['GENEINFO'].split(':')[0])
        else:
          GENE.append("gene information not present")
      for i in range(len(list(rec.info))):
        #info.append(rec.info[list(rec.info)[i]])
        strings = [list(rec.info)[i], rec.info[list(rec.info)[i]]] # it works on a single record. It creates a list of INFO_ID and INFO value.
        info.append('='.join(map(str, strings)))                   # concatenate values in a string e.g 'RS"=rs123 
      INFO.append(','.join(map(str, info)))                        # concatenate the string in an INFO list for all teh entries 
  os.chdir(actual_sir)    
  shutil.rmtree(dir, ignore_errors=False)

  
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
      PATHO.append(splitted[-2])


# From the lists create a dataframe
df = pd.DataFrame({'CHROM':CHROM, 
                   'POS':POS, 
                   'REF':REF, 
                   'ALT':ALT,
                   'GENE':GENE,
                   'PATHOGENICITY':PATHO,
                   'SOURCE':sourx,
                   'INFO':INFO})



# Remove duplicates variant entries and
# Pivot the dataframe to have source of information on the same line. 
df2 = df.groupby(['CHROM', 'POS', 'REF', 'ALT','GENE']).agg({'SOURCE' : lambda x: ','.join(map(str, x)), 
                                                              'PATHOGENICITY': lambda y: ',' .join(map(str, y)),
                                                              'INFO': lambda z: ';' .join(map(str, z))})
df2 = df2.sort_values(by=['CHROM', 'POS'])

df2.to_csv('/home/ls760/aggregated_karyn_resources.csv', sep='\t')

# Number of variants in the 
n_variatn_df = pd.read_csv("~/Desktop/VarioPath/number_variants_in_MDT_genes_clean.tsv",  delimiter='\t', skiprows=1, dtype={'n_karyn_list':np.int32, 'n_variant_in_UKB':np.int32}) 
n_variatn_df.head
n_variatn_df.astype({'n_karyn_list': 'int32', 'n_variant_in_UKB': 'int32' })

var_ori=pd.read_csv("~/Desktop/VarioPath/TO_REMOVE", delimiter=',', skiprows=0) 
var_ori.columns = ('gene', 'n_variant_unfiltered')

n_variatn_df.gene=n_variatn_df.gene.astype(str)
var_ori.gene=var_ori.gene.astype(str)
df4=n_variatn_df.merge(var_ori, how='left')

moi=pd.read_csv("/Users/luca/Desktop/VarioPath/MDT_variants/aggregated_spreadsheet.tsv", delimiter='\t', skiprows=0) 
moi = moi[['GENE', 'MOI_original_column']].drop_duplicates(keep='first')

df4=df4.merge(moi, how='left', left_on='gene', right_on='GENE')

# created groups with genes fro every MDT
blee_coag=["F10","F11","F12","F13A1","F13B",
           "F2","F5","F7","F8","F9","FGA","FGB",
           "FGG","GGCX","KNG1","LMAN1","MCFD2",
           "SERPINE1","SERPINF2","VKORC1","VWF"]
thrombosis=["ADAMTS13","HRG","PIGA","PLG","PROC",
		          "PROS1","SERPINC1","SERPIND1","THBD"]
platelet=["ABCC4","ABCG5","ABCG8","ACTB","ACTN1","ANKRD26",
          "ANO6","AP3B1","AP3D1","ARPC1B","BLOC1S3","BLOC1S6",
          "CDC42","CYCS","DIAPH1","DTNBP1","ETV6","FERMT3","FLI1",
          "FLNA","FYB1","GATA1","GFI1B","GNE","GP1BA","GP1BB","GP6",
          "GP9","HOXA11","HPS1","HPS3","HPS4","HPS5","HPS6","IKZF5",
          "ITGA2B","ITGB3","KDSR","LYST","MECOM","MPIG6B","MPL","MYH9",
          "NBEA","NBEAL2","P2RY12","PLA2G4A","PLAU","RASGRP2","RBM8A",
          "RNU4ATAC","RUNX1","SLFN14","SRC","STIM1","STXBP2","TBXA2R",
          "TBXAS1","THPO","TUBB1","VIPAS39","VPS33B","WAS"]
hered_sfero = ["ANK1","EPB41","EPB42","SLC4A1","SPTA1","SPTB"]

mdt_list=['blee_coag', 'thrombosis', 'platelet', 'hered_sfero']

df4['MDT']='MDT_value'
for i in mdt_list:
  for j in range(0,len(df4['gene'])):
    if df4.gene[j] in eval(i):
      df4.MDT[j]=str(i)

sns.set_style("whitegrid")  #set polt style
plot1=sns.relplot(x="n_karyn_list", y="n_variant_in_UKB", hue="MOI_original_column", col='MDT', data = df4, legend=True, size="MOI_original_column", sizes=(20,20), facet_kws={'sharey': False, 'sharex': False})
plot1.set_axis_labels('pathogenic and likely pathogenic variants', 'pathogenic and likely pathogenic variants in UK BioBank cohort')
for i in range(len(plot1.axes[0])):
  a = pd.concat({'x': df4.n_karyn_list, 'y': df4.n_variant_in_UKB, 'val': df4.gene}, axis=1)
  for j, point in a.iterrows():
    if point.val in eval(list(plot1.axes_dict.keys())[i]):
      if point.x > np.quantile(df4.n_karyn_list,0.7) or point.y > np.quantile(df4.n_variant_unfiltered,0.7):
        plot1.axes[0][i].text(point['x']+random(), point['y']+random(), str(point['val']), size=7)
plt.savefig("/Users/luca/Desktop/VarioPath/MDT_variants/sanity_check/plot1_filtered_and_UKB_cohort_20201214.png", 
            dpi="figure",
            orientation='landscape', 
            format="png",
            transparent=True, 
            bbox_inches='tight')

plot2=sns.relplot(x="n_karyn_list", y="n_variant_unfiltered",  hue="MOI_original_column", col='MDT', data = df4, legend=True,size="MOI_original_column", sizes=(20,20), facet_kws={'sharey': False, 'sharex': False})
plot2.set_axis_labels('pathogenic and likely pathogenic variants', 'all variants in list')
for i in range(len(plot2.axes[0])):
  a = pd.concat({'x': df4.n_karyn_list, 'y': df4.n_variant_unfiltered, 'val': df4.gene}, axis=1)
  for j, point in a.iterrows():
    if point.val in eval(list(plot2.axes_dict.keys())[i]):
      if point.x > np.quantile(df4.n_karyn_list,0.7) or point.y > np.quantile(df4.n_variant_unfiltered,0.7):
        plot2.axes[0][i].text(point['x']+random(), point['y']+random(), str(point['val']), size=7)
plt.savefig("/Users/luca/Desktop/VarioPath/MDT_variants/sanity_check/plot2_filtered_not_filtered_20201214.png", 
            dpi="figure",
            orientation='landscape', 
            format="png",
            transparent=True, 
            bbox_inches='tight')

plot3=sns.relplot(x="n_variant_in_UKB", y="n_variant_unfiltered",  hue="MOI_original_column", col='MDT', data = df4, legend=True,size="MOI_original_column", sizes=(20,20), facet_kws={'sharey': False, 'sharex': False})
plot3.set_axis_labels('pathogenic and likely pathogenic variants in UK BioBank cohort', 'all variants in list')
for i in range(len(plot3.axes[0])):
  a = pd.concat({'x': df4.n_variant_in_UKB, 'y': df4.n_variant_unfiltered, 'val': df4.gene}, axis=1)
  for j, point in a.iterrows():
    if point.val in eval(list(plot3.axes_dict.keys())[i]):
      if point.x > np.quantile(df4.n_variant_in_UKB,0.7) or point.y > np.quantile(df4.n_variant_unfiltered,0.7):
        plot3.axes[0][i].text(point['x']+random(), point['y']+random(), str(point['val']), size=7)
plt.savefig("/Users/luca/Desktop/VarioPath/MDT_variants/sanity_check/plot3_not_filtered_and_UKB_cohort_20201214.png", 
            dpi="figure",
            orientation='landscape', 
            format="png",
            transparent=True, 
            bbox_inches='tight')