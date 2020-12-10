import pysam
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pydrive.auth import GoogleAuth
from pydrive.drive import GoogleDrive

gauth = GoogleAuth()
gauth.LocalWebserverAuth() # client_secrets.json need to be in the same directory as the script
drive = GoogleDrive(gauth)

ID_variants_dir='1Eu5NNuODO3uYIozSPRE3oHJSoHiUkD81'
fileList = drive.ListFile({'q': "'1Eu5NNuODO3uYIozSPRE3oHJSoHiUkD81' in parents and trashed=false"}).GetList()
for file in fileList:
  print('Title: %s, ID: %s' % (file['title'], file['id']))
  if(file['title'] == 'VarioPath_variants_20200615.tsv'):
    file.GetContentFile('VarioPath_variants_20200615')


df_curatedPub = pysam.TabixFile('VarioPath_variants_20200615',  )
# Read in curatedPub file
CHROM = []
POS = []
REF = []
ALT = []
sourx = []

for rec in df_curatedPub.fetch():
    CHROM.append(rec.contig) 
    POS.append(rec.pos) 
    REF.append(rec.REF)
    ALT.append(rec.ALT)
    source = 'curatedPub'

df = pd.DataFrame(data, columns=['CHROM', 'POS', 'REF', 'ALT'])

# Number of variants in the 
n_variatn_df = pd.read_csv("~/Desktop/VarioPath/number_variants_in_MDT_genes_clean.tsv",  delimiter='\t', skiprows=1, dtype={'n_karyn_list':np.int32, 'n_variant_in_UKB':np.int32}) 
n_variatn_df.head
n_variatn_df.astype({'n_karyn_list': 'int32', 'n_variant_in_UKB': 'int32' })

sns.relplot(x='n_karyn_list', y='n_variant_in_UKB', data = n_variatn_df, )
plt.show()
