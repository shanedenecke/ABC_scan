
# ABC_scan
Tool to scan ABC transporter species. Only tested in arthropods. See Denecke et. al 2021 (under review) for more infomation.
Designed to work on linux like systems. A web based application also performs the same function http://chrysalida.imbb.forth.gr:3838/ABC_scan/.
For any questions or comments please email me at shanedenecke@gmail.com.


## Dependencies
* R
	* tidyverse
	* data.table
	* ggtree
	* ggsci
	* agricolae
	* gplots
	* ape
	* treeio
	* argparser
	* seqinr
* python3
	* Biopython
	* re 
	* os 
	* sys
	* pandas
	
* HMMER package (http://hmmer.org/)
* MAFFT (https://mafft.cbrc.jp/alignment/software/)
* Custom_Applications (https://github.com/shanedenecke/Custom_Applications)
* BLAST+ (https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download)
* MEME suite (https://meme-suite.org/meme/)



## Instructions

First make sure all dependencies are installed in R and python and all programs are in your PATH. Also make sure that the "ABC_scan" folder is made an executable. 

To run the program simply call the "abc_scan.sh" file from your terminal environment. There are 3 necessary arguments


The full help options are displayed below. 

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
