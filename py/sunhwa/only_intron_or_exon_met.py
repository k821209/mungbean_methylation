#!/usr/bin/python3 
import numpy as np
import pandas as pd
file_CG_exon = 'S_total.CG.q.dict.exonmethylv.txt'
file_CHG_exon = 'S_total.CHG.q.dict.exonmethylv.txt'
file_CHH_exon = 'S_total.CHH.q.dict.exonmethylv.txt'

file_CG_intron = 'S_total.CG.q.dict.intronmethylv.txt'
file_CHG_intron = 'S_total.CHG.q.dict.intronmethylv.txt'
file_CHH_intron = 'S_total.CHH.q.dict.intronmethylv.txt'

df_CG_exon = pd.read_csv(file_CG_exon,sep='\t')
df_CHG_exon = pd.read_csv(file_CHG_exon,sep='\t')
df_CHH_exon = pd.read_csv(file_CHH_exon,sep='\t')

df_CG_intron = pd.read_csv(file_CG_intron,sep='\t')
df_CHG_intron = pd.read_csv(file_CHG_intron,sep='\t')
df_CHH_intron = pd.read_csv(file_CHH_intron,sep='\t')


df_CG = pd.merge(df_CG_exon,df_CG_intron,on='genename',how='outer',suffixes=('_exon','_intron'))
df_CHG = pd.merge(df_CHG_exon,df_CHG_intron,on='genename',how='outer',suffixes=('_exon','_intron'))
df_CHH = pd.merge(df_CHH_exon,df_CHH_intron,on='genename',how='outer',suffixes=('_exon','_intron'))

df_CG = df_CG.dropna(axis=0,how='any')
df_CHG = df_CHG.dropna(axis=0,how='any')
df_CHH = df_CHH.dropna(axis=0,how='any')

df_CG_t = df_CG[(df_CG['Methylation Level_intron'] == 0) & (df_CG['Methylation Level_exon'] > 0)]
print(df_CG_t)
