usage: PO_2_MLSA.py [-h] -po PO_RESULTFILE [-f FASTA_PATH] [-tp TEMP_PATH]
                    [-kt] [-am {muscle,clustalw,clustalw2,clustalo}]
                    [-ap ALIGNER_PATH] [-dg {none,all,flanking}] [-s]
                    [-t NTHREADS] [-gb {n,no,f,false,y,yes,t,true}]
                    [-gbp GBLOCKS_PATH] [-op OUT_PATH]
                    [-mt {raxml,raxml_bs,raxml_rapidbs,none}]
                    [-tbp TREEBUILDER_PATH] [-sd SEED_NR] [-bs NR_BOOTSTRAPS]
                    [-v]

==PO_2_MLSA.py v1.2 by John Vollmers==
Creates concatenated alignments of UNIQUE Genes with orthologues in comparison organisms for the construction of MLSA-based phylogenetic trees.
Optionally the resulting concatenated alignments may contain all gapped alignmentpositions or may be stripped either of ALL gapped positions or of all gapped positions in the flanking regions of each composite ortholog
This script is supposed to be part of a pipeline consisting of:
	A.)Conversion of Genbank/Embl-Files to ANNOTATED(!) Fastas using CDS_extractor.pl by Andreas Leimbach
	B.)Calculation of orthologs and paralogs using proteinortho5 (WITH the '-single' and '-self' arguments!)
	C.)The creation of concatenated MLSA-sequences based on:
		-the fasta sequences of step A
		-the proteinortho5-results from step B

The output-file will be in fasta format (but including gapped positions, so remember to use 'fasta_wgap' when loading into Arb!). However it's absolutely no problem to include other common output-alignmentformats on request!

optional arguments:
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
