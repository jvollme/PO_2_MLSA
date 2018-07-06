## PO_2_MLSA.py v1.5.2
Automizes the creation of concatenated alignments of __single copy__ genes with ___orthologues___ in comparison organisms for the construction of MLSA-based phylogenetic trees. Can optionally automize pylogenetic tree calculation using [RaXML][]

Single-copy markers present in EVERY comparison organism EXACTLY ONCE are identified from a pre-calculated [proteinortho5][] resultfile, and extracted from corresponding supplied protein-fastas

Unalignable C- and N-terminal fregments will be clipped from each detected markergene and remaining poorly aligned positions will be automatically filtered from each single marker alignment using [Gblocks][] 

Optionally the resulting concatenated alignments can be specified to contain all gapped alignment positions or to be stripped either of _all_ gapped positions or of only such gapped positions in the flanking regions of each composite ortholog

This script is supposed to be part of a pipeline consisting of:

1. Extraction of protein seqeunces (with __unique__ identifiers!) fom Genbank/Embl-Files to one multi-fasta per comparison-organism using, e.g. [CDS_extractor.pl][] by Andreas Leimbach

2. Calculation of orthologs and paralogs using [proteinortho5][] (with the _'-single'_ and _'-self'_ arguments!)

3. The creation of concatenated MLSA-sequences based on:
    - the fasta sequences of step 1.
    - the [proteinortho5][]-results from step 2.

The output-file will be in _fasta_ format (but including gapped positions, so remember to use _'fasta_wgap'_ if and when loading into [Arb][].

#### Dependancies:
 - [Python v2.7+][]
 - [Biopython][]
 - [proteinortho5][]
 - [CDS_extractor.pl][]
 - [muscle][]
 - [Gblocks][]
 - [RaXML][] (optional, only when creating trees with RaXML)

#### Usage:  
````
usage: PO_2_MLSA.py [-h] -po PO_FILE [-f FASTA_PATH] [-t NTHREADS] [--filter]
                    [-op OUT_PATH] [-tp TEMP_PATH] [-kt] [-mb MUSCLE_BINARY]
                    [-gbb GBLOCKS_BINARY] [-rmlb RAXML_BINARY]
                    [-mt {raxml,raxml_rapidbs,nj,nj_bs,none}] [-sd SEED_NR]
                    [-bs NR_BOOTSTRAPS]
                    [--nj_substmodel]
                    [-v] [-s] [--debug]

````

#### Arguments:

##### Basic arguments:
````
  -h, --help            show a help message and exit
  -po PO_FILE, --proteinortho PO_FILE
                        (String) file with proteinortho5 results
  -f FASTA_PATH, --fastas FASTA_PATH
                        (String) path to fasta files (produced by CDS_extractor) 
                        Default = current working directory
  -t NTHREADS, --threads NTHREADS
                        Maximum number of threads to use for alignment steps
                        Default = 1
  --filter              OPTIONAL: Filter outliers based on single marker phylogenies in order to remove HGT events (Default: False)
````

##### Output-related arguments:
````
  -op OUT_PATH, --out_path OUT_PATH
                        Path to output (will be created if it does not exist)
                        Default = current working directory
  -tp TEMP_PATH, --temp_path TEMP_PATH
                        Path for temporary files (will be created if does not exist)
                        Default =  current working directory
  -kt, --keep_temp      Keep all temporary files and intermediate results
                        (Multifasta_files containing the singe genes involved in the MLSA-alignments are stored in any case)
                        Default: delete all temporary files
````

##### Dependancies/external binaries-related arguments:
````
  -mb MUSCLE_BINARY, --muscle_binary MUSCLE_BINARY
                        Muscle binary to use (including path if not in $PATH). default: assume "muscle" in $PATH
  -gbb GBLOCKS_BINARY, --gblocks GBLOCKS_BINARY
                        Gblocks binaries to use (include path to binaries if not listed in $PATH). Default: "Gblocks" (assumed to be in PATH)
  -rmlb RAXML_BINARY, --raxml RAXML_BINARY
                        RaxML excecutable biaries to use (Include path to binary if not listed in $PATH). Default: check for common naming of raxml binaries in $PATH, prioritizing binaries with "PTHREADS" in the name
````

##### Tree-calculation-related arguments:
````
  -mt {raxml,raxml_rapidbs,nj,nj_bs,none}, --make_tree {raxml,raxml_rapidbs,nj,nj_bs,none}
                        Generate ML phylogenetic trees using RAxML with the substitution model "PROTGAMMAAUTO"
                        choices: "raxml": single tree without bootstraps (using new rapid hill climbing)
                        	 "raxml_rapidbs": rapid bootstrap analyses and search for best ML tree in one run
                        	 "none"
                        Default = none
  -sd SEED_NR, --seed SEED_NR
                        Integer to provide as seed for RAxML
                        0 = seed generated randomly
                        Default = random seed
  -bs NR_BOOTSTRAPS, --bootstraps NR_BOOTSTRAPS
                        Number of bootstraps(if any)
                        default = 1000
  --nj_substmodel {identity,ident,blosum,pam,blosum30,blosum35,blosum40,blosum45,blosum50,blosum55,blosum60,blosum62,blosum65,blosum70,blosum75,blosum80,blosum85,blosum90,blosum95,blosum100,pam30,pam60,pam90,pam120,pam180,pam250,pam300}
                        Substitution model for distance matrix calculation. Only "identity" and variances of "blosum" and "pam" are offered. 
                        abbreviations:
                        	"identity"
                        	"blosum" = blosum62
                        	"pam" = pam120, default = "identity"
````

##### Other arguments:
````
  -v, --version         show version information and then quit (don't run complete script)
  -s, --silent          non-verbose mode
  --debug               Log extra info for debugging
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
[biopython]: https://biopython.org/
[Python v2.7+]: https://www.python.org/
