#!/usr/bin/env Rscript

########### Import packages
shhh <- suppressPackageStartupMessages
shhh(library(data.table))
shhh(library(tidyverse))


########### Set up Arguments
args = commandArgs(trailingOnly=TRUE)
#setwd('/mnt/disk/shane/Transporter_ID/Arthropod_ABC_pipeline')
#args[1]='/mnt/disk/shane/Transporter_ID/Arthropod_ABC_pipeline/ABC_search/TesSpe/"total_ABC_recip_blast.tsv'



########### Setup functions
abc.number=function(x){
  l=list()
  for(i in unique(x$family)){
    sub=x[family==i]
    sub$family=paste0(sub$family,'_',1:nrow(sub))
    l[[i]]=sub
  }
  return(rbindlist(l)$family)
}


########### Parse reciprical BLAST table
recip_blast=fread(args[1]) %>% select(V1,V2,V4,V5) %>% rename(query=V1,subject=V2,evalue=V4,qlen=V5)
recip_blast$subject=gsub('^.+__(.+)__.+$','\\1',recip_blast$subject)
species=unlist(strsplit(args[1],'\\/'))[3]
target.fam=as.character(str_match(unlist(strsplit(args[1],'\\/'))[4],'ABC[A-Z]+'))

########### Initialize empty lists for "for" loop
abc.total=list()
filter.list=list()
unsorted.list=list()


########### For loop for sorting into families
for(i in unique(recip_blast$query)){
  
  ## Extract relevant info from BLAST table
  sub=recip_blast[query==i]
  top_hit_fam=sub[query==i]$subject[1]
  fams=table(sub$subject)
  evalues=sub$evalue
  qlen=sub$qlen[1]
  
    
  if(!grepl('ABC[A-Z]+',sub$subject[1])){ ## filter out any that have top blast hit not an ABC
    filter.list[[i]]=data.table(geneid=i,family=top_hit_fam)
  }else if(min(evalues)>1e-5){ #### remove any where top hit is above 1e-5
    filter.list[[i]]=data.table(geneid=i,family=top_hit_fam)
  }else if(max(fams)/sum(fams)>.7){ ### take 4/5 cases are from the same family
    abc.total[[i]]=data.table(geneid=i,family=names(fams)[which(fams==max(fams))])
  }else if(length(unique(na.omit(sub$subject[1:3])))==1){ ### where top 3 are from the same family
    abc.total[[i]]=data.table(geneid=i,family=top_hit_fam)
  } else if(evalues[1]==0 & grepl('ABC',top_hit_fam)){ ## Keep if evalue ==0
    abc.total[[i]]=data.table(geneid=i,family=top_hit_fam)
  }else if((evalues[2]/evalues[1]> 1e5) & grepl('ABC',top_hit_fam)){ ## Keep where top hit is ABC and overwhelmingly significant
    abc.total[[i]]=data.table(geneid=i,family=top_hit_fam)
  }else if(('ABCBF' %in% names(fams) & ('ABCBH' %in% names(fams)))){
    abc.total[[i]]=data.table(geneid=i,family='ABCBAMG') ### if we can't tell among ABC-B families then we classify them as AMG (ambiguous)
  }else if(qlen>300){ ## keep where 4 out of 5 are ABCs but can be any family
    abc.total[[i]]=data.table(geneid=i,family='ABC_Unsorted')
    unsorted.list[[i]]=sub
  }
  else{
    filter.list[[i]]=data.table(geneid=i,family=top_hit_fam)
  }
}
  
########### Clean Family sort output
filter.table=rbindlist(filter.list) 
abc.table=rbindlist(abc.total) %>% arrange(family) %>% data.table()
abc.table$family=paste0(species,'__',abc.table$family,'__',abc.table$geneid)


########### Create a table and file for ABCs which were filtered out
if(nrow(filter.table)>0){
  colnames(filter.table)=c('query','name') 
  filter.output=filter.table %>% merge(recip_blast,by='query')
  write.csv(filter.output,paste0(gsub('/total_ABC_recip_blast.tsv','',args[1]),'/ABC_filtered_BLAST_out.csv'))
}

########### Print Good ABC table to terminal
if(nrow(abc.table)>0){
  colnames(abc.table)=c('code','name')
  cat(format_csv(abc.table))
}
