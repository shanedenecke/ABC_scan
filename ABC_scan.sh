#!/usr/bin/env bash

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

############ Describe help setting
if [ "$1" == "-h" ]; then
  echo "
  Welcome ! This shell script is designed to search ABC transporters in non-model arthropods
  
  Arguments:
  
  -proteome: Path to folder containing one or more Arthropod proteomes that you wish to search. 
  -threads: Self explanatory
  -outdir: Output diretory where all of your outputs will be located. Note: Directory will be created and files will be put into relevant subdiretories automatically
  
  Optional arguments
  
  -prefix: What do you want the prefix of your sequence names to be? Defaults to TesSpe (Test Species)
  -print: Do you want to print the final output to the terminal
  -hmm_profile: A HMM profile to search with. Defaults to PF00005
  -database:  protein database for reciprocal blast. Defaults to Drosophila_Tribolium_Human_Tetranychus database. 
  -minlen: The minimum length of proteins for filtering. Defaults to 250
  -motif: Motif to search for in MEME format. Defaults to 3 common ABC formats
  -domain_filter: Do you want to filter based on DOMAIN? default NA
  "
  exit 0
fi

############ Add arguments
while [ "$#" -gt 0 ]; do
  case "$1" in
    -proteome) PROTEOME="$2"; shift 2;;
    -hmm_profile) HMM_PROFILE="$2"; shift 2;;
    -database) DATABASE="$2"; shift 2;;
    -prefix) PREFIX="$2"; shift 2;;
    -threads) THREADS="$2"; shift 2;;
    -outdir) OUTDIR="$2"; shift 2;;
    -print) PRINT="$2"; shift 2;;
    -minlen) MINLEN="$2"; shift 2;;
    -motif) MOTIF="$2"; shift 2;;
    -domain_filter) DOMAIN_FILTER="$2"; shift 2;;
  esac  
done


############ Debugging for developer
#cd /mnt/disk/shane/Transporter_ID/Arthropod_ABC_pipeline
#PROTEOME=./proteomes/AedAeg_unigene.faa
#SCRIPT_DIR='/home/sdenecke/Applications/ABC_scan'

#cd /home/shanedenecke/Dropbox/quick_temp
#PROTEOME=./abc_compare/BemTab_unigene.faa
#SCRIPT_DIR='/home/shanedenecke/Applications/Custom_Applications/ABC_scan'
#PREFIX='TEST1'
#THREADS=2


############ Add Defaults
if [ -z $HMM_PROFILE ];then HMM_PROFILE=$SCRIPT_DIR/ABC_tran.hmm; fi
if [ -z $MOTIF ];then MOTIF=$SCRIPT_DIR/ABC_meme_clean.txt; fi
if [ -z $DATABASE ];then DATABASE=$SCRIPT_DIR/combined_marked_proteome.faa; fi
if [ -z $OUTDIR ];then OUTDIR="./ABC_search"; fi
if [ -z $PRINT ];then PRINT="NO"; fi
if [ -z $THREADS ];then THREADS=4; fi
if [ -z $MINLEN ];then MINLEN=250; fi
if [ -z $DOMAIN_FILTER ];then DOMAIN_FILTER=NA; fi
if [ -z $PREFIX ];then PREFIX=TesSpe; fi


############ Make working directory
mkdir -p $OUTDIR
mkdir $OUTDIR/$PREFIX
echo ' Now searching the '$PREFIX' proteome'


############ perform HMM search
hmmsearch --cpu $THREADS --notextw -E 100 $HMM_PROFILE $PROTEOME | sed -n '/Scores for complete sequences/,/------ inclusion threshold/p' | sed '$ d' | awk 'NR > 4 { print }' | awk '/^$/{exit} {print $0}' | sed -e "s/\s\{1,\}/\t/g" | cut -f 2- > $OUTDIR/$PREFIX/total_ABC_HMM.table


############ clean HMM search and output fasta
cut -f 9 $OUTDIR/$PREFIX/total_ABC_HMM.table | sed 's/\s+//g' | $SCRIPT_DIR/unigene_fa_sub.sh $PROTEOME - > $OUTDIR/$PREFIX/total_ABC_pre_blast.faa


############ run reciprocal blast
echo 'BLAST AWAY'
blastp -query $OUTDIR/$PREFIX/total_ABC_pre_blast.faa -db $DATABASE -outfmt "6 qseqid sseqid pident evalue qlen sstart send" -evalue 1e-3 -max_target_seqs 5 -max_hsps 1 -num_threads $THREADS > $OUTDIR/$PREFIX/total_ABC_recip_blast.tsv


############ sort into families 
Rscript $SCRIPT_DIR/ABC_Family_Sort.R $OUTDIR/$PREFIX/total_ABC_recip_blast.tsv > $OUTDIR/$PREFIX/total_ABC_dict.csv


############ get preliminary fasta files 
$SCRIPT_DIR/fasta_rename.py $PROTEOME $OUTDIR/$PREFIX/total_ABC_dict.csv > $OUTDIR/$PREFIX/total_ABC_marked_proteome.faa
grep -A 1 "__ABC" $OUTDIR/$PREFIX/total_ABC_marked_proteome.faa | sed '/--/d' > $OUTDIR/$PREFIX/Preliminary_ABCs.faa
rm $OUTDIR/$PREFIX/total_ABC_marked_proteome.faa


############ Search for presence of nucleotide binding domains (PF00005)
hmmsearch --domtblout $OUTDIR/$PREFIX/HMM_PF00005_output.tsv $HMM_PROFILE $OUTDIR/$PREFIX/Preliminary_ABCs.faa > $OUTDIR/$PREFIX/hmm_junk.txt
cat $OUTDIR/$PREFIX/HMM_PF00005_output.tsv | tail -n +4 | head -n -10 > $OUTDIR/$PREFIX/HMM_PF00005_output_clean.tsv

############ Search sequences for Motif with mast
mast -hit_list $MOTIF $OUTDIR/$PREFIX/Preliminary_ABCs.faa -notext -nohtml -nostatus | tail -n +3 | head -n -1 | sed -E 's/\s+/___/g' | sed 's/___/\t/g' | cut -f 1,3 > $OUTDIR/$PREFIX/Mem_motif.tsv


############ Make summary table
Rscript $SCRIPT_DIR/ABC_table_summarize.R --minlen 250 --motif $MOTIF --domain_filter $DOMAIN_FILTER --indir $(readlink -f $OUTDIR/$PREFIX)

############ Organize final output directory
mkdir $OUTDIR/$PREFIX/junk
mv $OUTDIR/$PREFIX/total_ABC_HMM.table $OUTDIR/$PREFIX/total_ABC_pre_blast.faa $OUTDIR/$PREFIX/ABC_filtered_out.csv $OUTDIR/$PREFIX/HMM_PF00005_output.tsv $OUTDIR/$PREFIX/HMM_PF00005_output_clean.tsv $OUTDIR/$PREFIX/Mem_motif.tsv $OUTDIR/$PREFIX/junk
mv $OUTDIR/$PREFIX/hmm_junk.txt $OUTDIR/$PREFIX/Preliminary_ABCs.faa $OUTDIR/$PREFIX/total_ABC_dict.csv $OUTDIR/$PREFIX/total_ABC_recip_blast.tsv $OUTDIR/$PREFIX/junk

############ Print Table to console if "Print" option selected
if [ $PRINT == "YES" ]; then cat $OUTDIR/$PREFIX/Final_ABC_table.tsv; fi


