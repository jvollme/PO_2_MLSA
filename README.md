## PO_2_MLSA.py v1.2
Automizes the creation of concatenated alignments of __single copy__ genes with ___orthologues___ in comparison organisms for the construction of MLSA-based phylogenetic trees. Can optionally automize pylogenetic tree calculation using [RaXML][]

Poorly aligned positions may be automatically filtered from the alignments using [Gblocks][] 

Optionally the resulting concatenated alignments can be specified to contain all gapped alignment positions or to be stripped either of _all_ gapped positions or of only such gapped positions in the flanking regions of each composite ortholog

This script is supposed to be part of a pipeline consisting of:

1. Conversion of Genbank/Embl-Files to __annotated__(!) fastas using [CDS_extractor.pl][] by Andreas Leimbach

2. Calculation of orthologs and paralogs using [proteinortho5][] (with the _'-single'_ and _'-self'_ arguments!)

3. The creation of concatenated MLSA-sequences based on:
    - the fasta sequences of step 1.
    - the [proteinortho5][]-results from step 2.

The output-file will be in _fasta_ format (but including gapped positions, so remember to use _'fasta_wgap'_ if and when loading into [Arb][].

####Dependancies:
 - [proteinortho5][]
 - [CDS_extractor.pl][]
 - [muscle][] | [clustalw][] | [clustalo][]
 - [Gblocks][]
 - [RaXML][]

####Usage:  
````
PO_2_MLSA.py [-h] -po PO_RESULTFILE [-f FASTA_PATH] [-tp TEMP_PATH]
             [-kt] [-am {muscle,clustalw,clustalw2,clustalo}]
             [-ap ALIGNER_PATH] [-dg {none,all,flanking}] [-s]
             [-t NTHREADS] [-gb {n,no,f,false,y,yes,t,true}]
             [-gbp GBLOCKS_PATH] [-op OUT_PATH]
             [-mt {raxml,raxml_bs,raxml_rapidbs,none}]
             [-tbp TREEBUILDER_PATH] [-sd SEED_NR] [-bs NR_BOOTSTRAPS]
             [-v]
````

####Optional arguments:
````
  -h, --help            show this help message and exit
  -po PO_RESULTFILE, --proteinortho PO_RESULTFILE
                        (String) file with proteinortho5 results
  -f FASTA_PATH, --fastas FASTA_PATH
                        (String) path to fasta files (produced by CDS_extractor) 
                        Default=current working directory
  -tp TEMP_PATH, --temp_path TEMP_PATH
                        Path for temporary files (will be created if does not exist)
                        Default= current working directory
  -kt, --keep_temp      Keep all temporary files and intermediate results
                        (Multifasta_files containing the singe genes involved in the MLSA-alignments are stored in any case)
                        Default: delete all temporary files
  -am {muscle,clustalw,clustalw2,clustalo}, --alignmethod {muscle,clustalw,clustalw2,clustalo}
                        The Alignmentmethod to use.
                        Default='muscle'
  -ap ALIGNER_PATH, --aligner_path ALIGNER_PATH
                        (OPTIONAL: set path to the aligner of choice IF not included in global PATH variable
  -dg {none,all,flanking}, --degap {none,all,flanking}
                        (Only meant for use, if Gblocks is not installed
                        Specify if and which gaps to remove:
                        	'none' keep all gapped positions in the final MLSA-alignment
                        	'all' remove ALL gapped postitions in the alignments
                        	'flanking' remove flanking gapped positions from all individual alignments
                        Default='all'
  -s, --silent          non-verbose mode
  -t NTHREADS, --threads NTHREADS
                        Maximum number of threads to use for alignment steps
                        Default=1
  -gb {n,no,f,false,y,yes,t,true}, --gblocks {n,no,f,false,y,yes,t,true}
                        calls gblocks (if installed) to remove gapped positions and poorly aligned regions
                        (Overrides '-dg'|'--degap'
                        choices:
                        	[n|no|f|false]: will NOT use Gblocks
                        	[y|yes|t|true]: WILL use Gblocks
                        Default=true (WILL use Gblocks)
  -gbp GBLOCKS_PATH, --gblocks_path GBLOCKS_PATH
                        (OPTIONAL: set path to Gblocks IF not included in global PATH variable)
  -op OUT_PATH, --out_path OUT_PATH
                        Path to output (will be created if it does not exist)
                        Default=current working directory
  -mt {raxml,raxml_bs,raxml_rapidbs,none}, --make_tree {raxml,raxml_bs,raxml_rapidbs,none}
                        Generate ML phylogenetic trees using RAxML with the substitution model "PROTGAMMAAUTO"
                        	choices:	"raxml": single tree without bootstraps (using new rapid hill climbing)
                        		raxml_bs: thorough bootstrap analyses and search for best ML tree
                        		raxml_rapidbs: rapid bootstrap analyses and search for best ML tree in one run
                        		none
                        Default=none
  -tbp TREEBUILDER_PATH, --tree_builder_path TREEBUILDER_PATH
                        Path to treebuilder (currently only raxml supported) if not listed in $PATH
  -sd SEED_NR, --seed SEED_NR
                        Integer to provide as seed for RAxML
                        0=seed generated randomly
                        Default=random seed
  -bs NR_BOOTSTRAPS, --bootstraps NR_BOOTSTRAPS
                        Number of bootstraps(if any)
                        default=1000
  -v, --version         show version information and then quit (don't run complete script)
````
[proteinortho5]: https://www.bioinf.uni-leipzig.de/Software/proteinortho/
[CDS_extractor.pl]: https://github.com/aleimba/bac-genomics-scripts.git
[muscle]: http://www.drive5.com/muscle/
[clustalw]: http://www.clustal.org/clustal2/
[clustalw2]: http://www.clustal.org/clustal2/
[clustalo]: http://www.clustal.org/omega/
[gblocks]: http://molevol.cmima.csic.es/castresana/Gblocks.html
[raxml]: http://sco.h-its.org/exelixis/web/software/raxml/index.html
[arb]: http://www.arb-home.de/
