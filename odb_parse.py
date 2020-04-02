#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May 18 17:20:39 2019

Script to parse ODB files. Can either find one to one orthologues or one to many.

@author: shanedenecke
"""

########## Import modules
import os
import pandas as pd
import argparse
import warnings
from Bio import SeqIO
import gc 
import shutil 
from pathlib import Path

sys_home = str(Path.home())


########## Read in ARGS
CLI=argparse.ArgumentParser()
CLI.add_argument("-node",type=str,default='Metazoa',help='Choose node to cluster species at. Can be either "Metazoa or Arthropod"')
CLI.add_argument("-taxid",nargs="*",type=str,default=['7227_0','9606_0','29058_0'],help='put in space separate list of taxid. Should be in the format NUMBER_0. e.g. "7227_0" for Drosophila')
CLI.add_argument("-output",type=str,default='id',help='Choose the output type. For a table of IDs put "id". For sequences put "seq"')
CLI.add_argument("--algorithm",type=str,default='oto',help='Choose what parsing algorithm you want to use. oto=one to one; otm=one to many. mtm: many to many"')
CLI.add_argument("--home",type=str,default=sys_home+'/Applications/Custom_Applications/OrthoDB_source/',help='Choose what parsing algorithm you want to use. oto=one to one; otm=one to many (only can be used in pairwise comparisons"')

args = CLI.parse_args()


#args.algorithm='mtm'

########## Define Functions
warnings.filterwarnings("ignore")
def unicount(series):
    return len(list(set(series)))

def database_find(taxid,data):
    x=set(list(data[data['taxid']==taxid].loc[:,'database']))
    if 'ENSEMBL' in x:
        return('ENSEMBL')
    elif'NCBIgid' in x:
        return('NCBIgid')
    elif'UniProt' in x:
        return('UniProt')
    elif'NCBIproteinGI' in x:
        return('NCBIproteinGI')
    elif'InterPro' in x:
        return('InterPro')
    elif'Trinity' in x:
        return('Trinity')
    else:
        print('NO MATCH')




########### Import arguments
home=args.home
out_type=args.output
algo=args.algorithm

if(args.node=='Arthropod' or args.node=='arthropod'):
    odb_node=pd.read_csv(home+'odb10v0_OG2genes.6656.tab',sep='\t',header=None)    
if(args.node=='Metazoa' or args.node=='metazoa'):
    odb_node=pd.read_csv(home+'odb10v0_OG2genes.33208.tab',sep='\t',header=None)
if len(args.taxid)==1:
    taxid_list=[line.rstrip('\n') for line in open(args.taxid[0])]
else:
    taxid_list=args.taxid




########### Import and clean Key data
odb_key=pd.read_csv(home+'/odb10v1_gene_xrefs_extended_arth_sub.tab',sep='\t',header=None)

if len(taxid_list)>3:
    taxid_names=''
else: 
    taxid_names="_".join(taxid_list)

odb_key.columns=['odb','xref','database']
odb_key[['taxid','gene']]=odb_key['odb'].str.split(':',expand=True) ## split taxid and gene columns 
odb_key_sp=odb_key[odb_key['taxid'].isin(taxid_list)] ### subset for individual species


########## Get databases for each species and add to odb_key
dbs=[database_find(x,odb_key_sp) for x in taxid_list]
dbs_dict=dict(zip(taxid_list,dbs))

dataframe_list_db_add=[] ### loop to add database column to dataframes
for k, v in dbs_dict.items():
   dataframe_db_add=odb_key_sp[((odb_key['database']==v) & (odb_key['taxid']==k))]
   dataframe_list_db_add.append(dataframe_db_add)

odb_key_sub=pd.concat(dataframe_list_db_add) ## merge dataframes with database added
odb_key_sub=odb_key_sub.drop(['database','gene','taxid'],axis=1) ## get rid of useless columns
odb_key_sub.drop_duplicates(subset='odb',keep='first',inplace=True)

del odb_key_sp


##########  Import noe data
odb_node[['taxid','gene']]=odb_node[1].str.split(':',expand=True) ### split odb column into taxid and ortho group
odb_node.columns=['OG','junk','taxid','gene']
node_tax_subset=odb_node[odb_node['taxid'].isin(taxid_list)]
ortho_counts=node_tax_subset.groupby('OG')['taxid'].agg({'total_genes':'count','unique_taxids':unicount}) ## aggregate function
del odb_node 


gc.collect()


########### Perform Filtering Algorithms

if algo=='oto': ### Filter for one to one orthologues
    all_species=ortho_counts[ortho_counts['total_genes']==len(taxid_list)] ## filter for all genes that are present in all species listed in taxid file
    one_to_one=all_species[all_species['unique_taxids']==len(taxid_list)] ## filter for OGs which have exactly same number of total genes as species they are present in
    og_groups=[str(x) for x in list(one_to_one.index)]
    og_groups=list(one_to_one.index)
    og_genes=node_tax_subset[node_tax_subset['OG'].isin(og_groups)]
    og_genes.rename(columns={'junk':'odb'},inplace=True)
    del node_tax_subset
elif algo=='otm':
    prim_species=taxid_list[0]
    prim_ogs=list(node_tax_subset[node_tax_subset.taxid==prim_species].OG)
    prim_ogs2=[x for x in prim_ogs if prim_ogs.count(x)==1]
    one_to_many=ortho_counts.loc[prim_ogs2,:]
    og_groups=list(one_to_many.index)
    og_genes=node_tax_subset[node_tax_subset['OG'].isin(og_groups)]
    og_genes.rename(columns={'junk':'odb'},inplace=True)
    del node_tax_subset
elif algo=='mtm':
    prim_species=taxid_list[0]
    many_to_many=ortho_counts
    og_groups=list(many_to_many.index)
    og_genes=node_tax_subset[node_tax_subset['OG'].isin(og_groups)]
    og_genes.rename(columns={'junk':'odb'},inplace=True)
    del node_tax_subset




########### Give Output

if out_type=='id':
    og_final=og_genes.drop(['gene'],axis=1)
    m=pd.merge(og_final,odb_key_sub,on='odb',how='left')
    m['xref']=m.apply(lambda row: row['xref'] if type(row['xref'])==str else row['odb'],axis=1)
    m=m.drop(['odb'],axis=1)
    
    if algo=='oto':
        final=pd.pivot_table(m,index='OG',values='xref',columns='taxid',aggfunc=lambda x: x) ## pivot table to simply get all otos
    if algo=='otm' or algo=='mtm':
        merge_dict={}
        for i in taxid_list:
            if i==prim_species:
                final=pd.pivot_table(m[m.taxid==i],index='OG',values='xref',columns='taxid',aggfunc=lambda x: x)
            else:
                temp=m[m.taxid==i]
                temp.rename(columns={'xref':i},inplace=True)
                temp.drop('taxid',axis=1,inplace=True)
                merge_dict.update({i:temp})
        for k,v in merge_dict.items():
            if algo=='otm':
                final=pd.merge(final,v,on='OG',how='left')
            elif algo=='mtm':
                final=pd.merge(final,v,on='OG',how='outer')
                #final[prim_species]=[float('NaN') if type(x)==numpy.ndarray else x for x in list(final[prim_species])]
          
    ### output final table
    final['OG']=final.index
    final.to_csv('temp.tsv',index=False,sep='\t')
    os.system('cat temp.tsv')  
        
            
if out_type=='seq':
        
    
    fa_sub={} 
    og_list=set(og_genes['odb'])
    for recs in SeqIO.parse(home+'/odb_arthropoda_augment.faa', 'fasta'):
        if recs.id in og_list:
            fa_sub.update({recs.description:recs.seq})

### write ortho group fastas
    try:
        shutil.rmtree('./og_sequences'+taxid_names) ### remove directory if already exists
    except:
        pass
    
    os.mkdir('./og_sequences'+taxid_names)
    for i in og_groups:
        
        og_sub=list(og_genes[og_genes['OG']==i]['odb'])
        og_fasta={k:v for (k,v) in fa_sub.items() if k in og_sub}
        
        with open('./og_sequences'+taxid_names+'/'+str(i)+'.faa', 'w') as fp:
            for k,v in og_fasta.items():
                fp.write('>')
                fp.write(k)
                fp.write('\n')
                fp.write(str(v))
                fp.write('\n')