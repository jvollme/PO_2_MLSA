#!/usr/bin/python
#created 21.06.2014 by John Vollmers
import os, sys, logging, argparse, time, multiprocessing, random, traceback
from Bio.Align.Applications import ClustalwCommandline, MuscleCommandline, ClustalOmegaCommandline
from Bio.Phylo.Applications import RaxmlCommandline
from Bio import AlignIO, SeqIO
from Bio import Phylo
from Bio.Alphabet import IUPAC
from subprocess import call
from Bio.Phylo.Consensus import *
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor

myparser=argparse.ArgumentParser(description="\n==PO_2_MLSA.py v1.3 by John Vollmers==\nCreates concatenated alignments of UNIQUE Genes with orthologues in comparison organisms for the construction of MLSA-based phylogenetic trees.\nOptionally the resulting concatenated alignments may contain all gapped alignmentpositions or may be stripped either of ALL gapped positions or of all gapped positions in the flanking regions of each composite ortholog\nThis script is supposed to be part of a pipeline consisting of:\n\tA.)Conversion of Genbank/Embl-Files to ANNOTATED(!) Fastas using CDS_extractor.pl by Andreas Leimbach\n\tB.)Calculation of orthologs and paralogs using proteinortho5 (WITH the '-single' and '-self' arguments!)\n\tC.)The creation of concatenated MLSA-sequences based on:\n\t\t-the fasta sequences of step A\n\t\t-the proteinortho5-results from step B\n\nThe output-file will be in fasta format (but including gapped positions, so remember to use 'fasta_wgap' when loading into Arb!). However it's absolutely no problem to include other common output-alignmentformats on request!", formatter_class=argparse.RawTextHelpFormatter)
myparser.add_argument("-po", "--proteinortho", action = "store", dest = "PO_file", help = "(String) file with proteinortho5 results", required = True)
myparser.add_argument("-f", "--fastas", action = "store", dest = "fasta_path", help = "(String) path to fasta files (produced by CDS_extractor) \nDefault = current working directory", default = "")
myparser.add_argument("-tp", "--temp_path", action = "store", dest = "temp_path", default = "", help = "Path for temporary files (will be created if does not exist)\nDefault =  current working directory")
myparser.add_argument("-kt", "--keep_temp", action = "store_true", dest = "keep_temp", default = False, help = "Keep all temporary files and intermediate results\n(Multifasta_files containing the singe genes involved in the MLSA-alignments are stored in any case)\nDefault: delete all temporary files")
myparser.add_argument("-am", "--align_method", action = "store", dest = "alignmeth", choices = ["muscle", "clustalw", "clustalw2", "clustalo"], default = "muscle", help = "The Alignmentmethod to use.\nDefault = 'muscle'")
myparser.add_argument("-ap", "--aligner_path", action = "store", dest = "aligner_path", default = "", help = "(OPTIONAL: set path to the aligner of choice IF not included in global PATH variable")
#myparser.add_argument("-dg", "--degap", action = "store", dest = "degap", choices = ["none", "all", "flanking"], default = "all", help = "(Only meant for use, if Gblocks is not installed\nSpecify if and which gaps to remove:\n\t'none' keep all gapped positions in the final MLSA-alignment\n\t'all' remove ALL gapped postitions in the alignments\n\t'flanking' remove flanking gapped positions from all individual alignments\nDefault = 'all'")
myparser.add_argument("-F", "--filter", action = "store", dest = "afilter", choices = ["none", "degap_all", "degap_flanking", "gblocks"], default = "gblocks", help = "Specify if and how alignments should be filtered:\n\t'none' do not filter alignments (keep all gapped positions)\n\t'degap_all' remove ALL gapped postitions in the alignments (use only if gblock is not available)\n\t'degap_flanking' remove flanking gapped positions from all individual alignments (use only if gblocks is not available)\n'gblocks' use gblocks with default settings to degap and filter alignments (Default)")
myparser.add_argument("-s", "--silent", action = "store_true", dest = "no_verbose", help = "non-verbose mode")
myparser.add_argument("-t", "--threads", action = "store", dest = "nthreads", type = int, default = 1, help = "Maximum number of threads to use for alignment steps\nDefault = 1")
#myparser.add_argument("-gb", "--gblocks", action = "store", dest = "gblocks", choices = ["n", "no", "f", "false", "y", "yes", "t", "true"], default = "true", help = "calls gblocks (if installed) to remove gapped positions and poorly aligned regions\n(Overrides '-dg'|'--degap'\nchoices:\n\t[n|no|f|false]: will NOT use Gblocks\n\t[y|yes|t|true]: WILL use Gblocks\nDefault = true (WILL use Gblocks)")
myparser.add_argument("-gbp", "--gblocks_path", action = "store", dest = "gblocks_path", default = "", help = "(OPTIONAL: set path to Gblocks IF not included in global PATH variable)")
myparser.add_argument("-op", "--out_path", action = "store", dest = "out_path", default = ".", help = "Path to output (will be created if it does not exist)\nDefault = current working directory")
myparser.add_argument("-mt", "--make_tree", action = "store", dest = "tree_method", choices = ["raxml", "raxml_bs", "raxml_rapidbs", "nj", "nj_bs", "none"], default = "none", help = "Generate ML phylogenetic trees using RAxML with the substitution model \"PROTGAMMAAUTO\"\n\tchoices:\t\"raxml\": single tree without bootstraps (using new rapid hill climbing)\n\t\traxml_bs: thorough bootstrap analyses and search for best ML tree\n\t\traxml_rapidbs: rapid bootstrap analyses and search for best ML tree in one run\n\t\tnone\nDefault = none")
myparser.add_argument("-tbp", "--tree_builder_path", action = "store", dest = "treebuilder_path", default = "", help = "Path to treebuilder (currently only raxml supported) if not listed in $PATH")
myparser.add_argument("-sd", "--seed", action = "store", dest = "seed_nr", type = int, default = 0, help = "Integer to provide as seed for RAxML\n0 = seed generated randomly\nDefault = random seed")
myparser.add_argument("-bs", "--bootstraps", action = "store", dest = "nr_bootstraps", type = int, default = 1000, help = "Number of bootstraps(if any)\ndefault = 1000")
#myparser.add_argument("-ctba", "--custom_tree_builder_args", action = "store", dest = "custom_tree_builder_args", default = None, help = "custom arguments for raxml. CAUTION: will overide Only use if you know ExACTLY what you are doing!")
myparser.add_argument("-v", "--version", action = "store_true", dest = "showversion", default = False, help = "show version information and then quit (don't run complete script)")
myparser.add_argument("--debug", action = "store_true", dest = "debug", default = False, help = "Log extra info for debugging")
myparser.add_argument("--nj_substmodel", action = "store", dest = "subst_model",\
 choices = ['identity', 'blosum', 'pam','benner6', 'benner22', 'benner74', 'blosum100',\
 'blosum30', 'blosum35', 'blosum40', 'blosum45', 'blosum50', 'blosum55', 'blosum60',\
 'blosum62', 'blosum65', 'blosum70', 'blosum75', 'blosum80', 'blosum85', 'blosum90',\
 'blosum95', 'feng', 'fitch', 'genetic', 'gonnet', 'grant', 'ident', 'johnson', 'levin',\
 'mclach', 'miyata', 'nwsgappep', 'pam120', 'pam180', 'pam250', 'pam30', 'pam300',\
 'pam60', 'pam90', 'rao', 'risler', 'structure'],\
  default = "ident", help = "Substitution model for distance matrix calculation. All models listed in Bio.Phylo.TreeConstruction.DistanceCalculator.protein_models are available. default settings:\n\t\"identity\"\n\t\"blosum\" = blosum62\n\t\"pam\" = pam120")
#myparser.add_argument("--existing_align", action = "store", dest = existing_align, default = None, help = "contignue calculation from existing concatenated (and filtered!) alignment generated with PO_2_MLSA.py
args = myparser.parse_args()

#TOdo: add option "return_selection" to store selection of MLSA genes as unagligned multifastas or only lists of fasta-headers

version = "v1.3"
available_cores = multiprocessing.cpu_count() #counts how many cores are available, to check if the user-argument for threads can be fulfilled
aln_length, proc_aln_length, OG_number = 0, 0, 0
wstrings, estrings, lstrings = [], [], [] #warning, error and log messages respectively

if not os.path.exists(args.out_path):
	os.makedirs(args.out_path)
outputfilename = os.path.join(args.out_path, "concatenated_orthologs_%s_%s.fasta" %(args.alignmeth, os.path.basename(args.PO_file)))
logfilename = os.path.join(args.out_path, "PO_2_MLSA_%s.log" % time.strftime("%Y%m%d%H%M%S"))
raxml_prog = "raxmlHPC"
verbose = True
docontinue = True
gblocks = True
if args.afilter != "gblocks":
	gblocks = False
gblocks_path = args.gblocks_path
erroroccured, warningoccured = False, False
hline = "-" * 50
if args.subst_model in ["identity", "blosum", "pam"]:
	if args.subst_model == "identity":
		args.subst_model = "ident"
	elif args.subst_model == "blosum":
		args.subst_model = "blosum62"
	elif args.subst_model == "pam":
		args.subst_model = "pam120"



def setlogger(logfilename): #shall replace datsanerror, datsanlogmessage  datswarning
	#format logger
	logFormatter = logging.Formatter("[%(levelname)s]  %(message)s")
	mylogger = logging.getLogger()
	if args.debug:
		mylogger.setLevel(logging.DEBUG)
	else:
		mylogger.setLevel(logging.INFO)
		
	#set logger to write to logfile
	fileHandler = logging.FileHandler(logfilename)
	fileHandler.setFormatter(logFormatter)
	mylogger.addHandler(fileHandler)
	
	#set loger to write to stderr also
	consoleHandler = logging.StreamHandler()
	consoleHandler.setFormatter(logFormatter)
	if args.no_verbose:
		consoleHandler.setlevel(logging.WARNING)
	mylogger.addHandler(consoleHandler)
	
	return mylogger

def checkargs():
	#print "CHECKING ARGS"
	global raxml_prog
	global args, verbose
	mylogger.debug("checkargs()")
		
	if args.no_verbose:
		verbose = False
		
	if available_cores < args.nthreads:
		mylogger.warning("less than %d cores/threads available! Setting nthreads to %d" %(args.nthreads, available_cores))
		args.nthreads = available_cores
		
	if args.fasta_path != "" and (not os.path.exists(args.fasta_path) or not os.path.isdir(args.fasta_path)):
		raise OSError("fasta_path: '%s' does not exist or is no directory!" % args.fasta_path)
		
	if args.aligner_path != "":
		if os.path.exists(args.aligner_path) and os.path.isdir(args.aligner_path):
			if not os.path.exists(os.path.join(args.aligner_path, args.alignmeth)):
				raise OSError("can't find '%s' at %s" %(args.alignmeth, args.aligner_path))
		else:
			raise OSError("aligner_path: '%s' does not exist or is not a directory!" % args.aligner_path)
	elif args.aligner_path == "":
		test_aligner = which(args.alignmeth)
		if test_aligner == None:
			raise OSError("Can't find %s in any directory in the PATH variable. Please provide a path" % args.alignmeth)
			
	if not os.path.exists(args.PO_file) or not os.path.isfile(args.PO_file):
		raise OSError("cannot find proteinortho-resultfile: %s" % args.PO_file)
		
	if args.temp_path != "" and docontinue:
		if not os.path.exists(args.temp_path) or not os.path.isdir(args.temp_path):
			if verbose:
				print "Creating directory for temporary and intermediate result files: %s" % args.temp_path
			os.mkdir(args.temp_path)
		mylogger.info("-Will store temporary and intermediate result files in %s" % os.path.abspath(args.temp_path))
		
	if args.out_path != "" and docontinue:
		if not os.path.exists(args.out_path) or not os.path.isdir(args.out_path):
			if verbose:
				print "Creating directory for final result files: %s" % args.out_path
			os.mkdir(args.out_path)
		mylogger.info("-Will store final result files in %s" % os.path.abspath(args.out_path))
		
	if args.alignmeth == "clustalo" and docontinue: #check clustalo version
		clustalomega_cline = ClustalOmegaCommandline("clustalo", version = True)
		clustalo_version = clustalomega_cline()[0].rstrip().split(".")
		if int(clustalo_version[0]) < 1 or (int(clustalo_version[0]) == 1 and int(clustalo_version[1]) < 2) : #only accept versions 1.2 and newer
			raise OSError("found clustalo version is v%s ! Version 1.2 or higher is required!" % ".".join(clustalo_version))
			
	if gblocks:
		if gblocks_path == "":
			test_gblocks = which("Gblocks")
			if test_gblocks == None:
				raise OSError("can't locate Gblocks in any path in PATH variable. please provide a Path to Gblocks using the '-gbp' agrument, or choose a different filtering option ('-F')")
			elif verbose:
				print "Located Gblocks executable: %s" % test_gblocks
		else:
			if os.path.exists(gblocks_path) and os.path.isdir(gblocks_path):
				if os.path.exists(os.path.join(gblocks_path, "Gblocks")) and os.path.isfile(os.path.join(gblocks_path, "Gblocks")):
					if verbose:
						print "Located Gblocks executable: %s" % os.path.join(gblocks_path, "Gblocks")
			elif gblocks_path.endswith("Gblocks") and os.path.isfile(gblocks_path):
				if verbose:
					print "Located Gblocks executable: %s" % gblocks_path
			else:
				raise OSError("Gblocks executable could not be found in the specified path: %s" % gblocks_path)
				
	if args.tree_method != "none":
		#print "CHECKING raxml-binaries"
		if args.treebuilder_path == "":
			if which("raxmlHPC") == None and which("raxmlHPC-PTHREADS") == None and which("raxmlHPC-PTHREADS-SSE3") == None and which("raxmlHPC-SSE3") == None:
				raise OSError("ERROR: No raxmlHPC binaries found in any directory within $PATH! please provide a PATH to raxml binaries!")
				args.tree_method = "none"
			else:
				if which("raxmlHPC-PTHREADS-SSE3") != None:
					raxml_prog = which("raxmlHPC-PTHREADS-SSE3")
					mylogger.info("found "+raxml_prog+ "! will use this tool")
				elif which("raxmlHPC-PTHREADS") != None:
					raxml_prog = which("raxmlHPC-PTHREADS")
					mylogger.info("found "+raxml_prog+ "! will use this tool")
				else:
					mylogger.warning("multithreading is only supported with 'PTHREADS'-versions of raxml. Not sure if your raxml binaries support this.\t\nif raxml calculations fail, recombile raxml with 'PTHREADS'-option")
					if which("raxmlHPC-SSE3") != None:
						raxml_prog = which("raxmlHPC-SSE3")
						mylogger.info("found "+raxml_prog+ "! will use this tool")
					else:
						raxml_prog = which("raxmlHPC")
						mylogger.info("found "+raxml_prog+ "! will use this tool")
				try:
					checkraxml_cline = RaxmlCommandline(raxml_prog, version = True)
					versiontext = checkraxml_cline()[0]
					startv = versiontext.find("version ") + len("version ")
					endv = versiontext[startv:].find(" ")
					version = versiontext[startv:startv+endv].split(".")
					if int(version[0]) < 8 or (int(version[0]) == 8 and int(version[1]) == 0 and int(version[2]) < 20):
						mylogger.warning("This script was devised for and tested with RAxML v8.0.20. Your version is v%s !\n\tThis may very well still work, but if it doesn't it's YOUR fault!" % ".".join(version))
				except:
					mylogger.warning("This script was devised for and tested with RAxML v8.0.20.\n\tNot sure which version of RAxML you're using, but it sure as hell isn't v7 or v8!\n\tThis may very well still work, but if it doesn't it's YOUR fault!")
		elif os.path.exists(args.treebuilder_path) and os.path.isfile(args.treebuilder_path):
			if args.nthreads > 1 and not "PTHREADS" in args.treebuilder_path:
				mylogger.warning("multithreading is only supported with 'PTHREADS'-versions of raxml. Not sure if your choosen binaries support this.\t\nif raxml calculations fail, recombile raxml with 'PTHREADS'-option")
			try:
				checkraxml_cline = RaxmlCommandline(version = True)
				versiontext = checkraxml_cline()[0]
				startv = versiontext.find("version ")+len("version ")
				endv = versiontext[startv:].find(" ")
				version = versiontext[startv:startv+endv].split(".")
				if int(version[0]) < 8 or (int(version[0]) == 8 and int(version[1]) == 0 and int(version[2]) < 20):
						mylogger.warning("This script was devised for and tested with RAxML v8.0.20. Your version is v%s !\n\tThis may very well still work, but if it doesn't it's YOUR fault!" % ".".join(version))
				raxml_prog = args.treebuilder_path
				
			except:
				mylogger.warning("Correct raxML-version not found under %s !\nWill NOT calculate ML trees!" % args.treebuilder_path)
				args.tree_method = "none"

def which(thisfile):
	mylogger.debug("which(%s)", thisfile)
	for path in os.environ["PATH"].split(":"):
		if os.path.exists(path + "/" + thisfile):
			return path + "/" + thisfile
	return None

def read_PO_file(filename):
	mylogger.debug("read_PO_file(%s)" % filename)
	open_PO_file = open(filename, 'r')
	firstline = True
	org_index = 3 #default index of beginning of organism-result-columns in proteinortho4+
	headers = []
	MLSA_gene_counter = 0
	MLSA_list = []
	
	for line in open_PO_file:
		if firstline:
			if line.startswith("#species"):
				mylogger.warning("It seems you are using results from Proteinortho4. It's recommended to use Proteinortho5 or later!\n\t(However, the derivation of MLSA genes should still work)")
			elif line.startswith("# Species"):
				if verbose:
					print "\nProteinortho-results seem to be based on Proteinortho5 or later. Good."
			else:
				mylogger.warning("Cannot clearly recognize Format of Proteinortho-results! This may produce erroneous results!\n\t For best results use Proteinortho5!")
			headers = line.split("\t")[org_index:]
			#remove line_end_symbols:
			for h in range(len(headers)):
				headers[h] = headers[h].rstrip()
			MLSA_list = [[header, []] for header in headers]
			firstline = False
			
		elif not line.startswith("#"):
			if not "," in line:
				if not "*" in line:
					zeilentokens = line.split("\t")[org_index:]
					#remove line end symbols:
					for zt in range(len(zeilentokens)):
						zeilentokens[zt] = zeilentokens[zt].rstrip()
					if len(zeilentokens) != len(headers):
						mylogger.error("ERROR: Different numbers of columns in header_line and Locus_tag_lines in "+filename)
						break
					for h in range(len(headers)):
						MLSA_list[h][1].append(zeilentokens[h].rstrip())
						
	open_PO_file.close()
	#check if everything is ok before continuing
	for org in MLSA_list:
		if len(org[1]) != len(MLSA_list[0][1]):
			raise OSError("something went wrong while reading proteinortho-results. Not the same number of comparison genes for all Organisms!")
			
	mylogger.info("\nrecognized %d unambigious ( = unique) 'Orthologeous Groups' (OGs) shared by all %d comparison organisms" %(len(MLSA_list[0][1]), len(headers)))
	OG_number = len(MLSA_list[0][1])
	
	return headers, MLSA_list, OG_number

def read_fasta_seqs(headers, MLSA_list):
	mylogger.debug("read_fasta_seqs(headers, MLSA_list)")
	record_dict = dict.fromkeys(headers, {})
	
	for org in record_dict:
		record_dict[org] = {}
		
	for h in range(len(headers)):
		mylogger.debug("reading file: %s" % headers[h])
		if MLSA_list[h][0] != headers[h]:
			mylogger.error("ERROR: MLSA_list and headers are not synchronized! WTF!?!?!?!?")
			break
		if not os.path.exists(os.path.join(args.fasta_path, headers[h])) or not (os.path.isfile(os.path.join(args.fasta_path, headers[h])) or os.path.islink(os.path.join(args.fasta_path, headers[h]))):
			mylogger.error("ERROR: Can't find fasta file: %s !\n\tPlease provide a path to a directory, that contains all fasta-files listed in %s" %(os.path.join(args.fasta_path, headers[h]), args.PO_file))
			return
		open_fastafile = open(os.path.join(args.fasta_path, headers[h]), 'r')
		temprecord_dict = dict.fromkeys(MLSA_list[h][1], None)
		for record in SeqIO.parse(open_fastafile, "fasta", alphabet = IUPAC.protein):
			if record.id in temprecord_dict:
				record_dict[headers[h]][record.id] = record
		#add another check_function here!!
		open_fastafile.close()
		
	return record_dict

def write_seq_files(record_dict, MLSA_list):
	mylogger.debug("write_seq_files(record_dict, MLSA_list)")
	selection_filelist = []

	for orgindex in range(len(MLSA_list)):
		tempfilename = os.path.join(args.out_path, "OGselection_MLSA_single_unaligned_multifasta_%s" % MLSA_list[orgindex][0])
		selection_filelist.append(tempfilename)
		open_tempfile = open(tempfilename, 'w')
		for locus in MLSA_list[orgindex][1]:
			SeqIO.write([record_dict[MLSA_list[orgindex][0]][locus]], open_tempfile, "fasta")
		open_tempfile.close()
		
	return selection_filelist

def write_temp_files(record_dict, MLSA_list, prefix):
	mylogger.debug("write_temo_files(record_dict, MLSA_list, %s)" % prefix)
	unaligned_filelist = []
	OG_counter = 0 #numbering of the used 'Orthologeous Groups'(OGs), to differentiate between the different temp_files
	
	for locusindex in range(len(MLSA_list[0][1])):
		OG_counter += 1
		tempfilename = os.path.join(args.temp_path, prefix + str(OG_counter).zfill(5) + ".fasta")
		unaligned_filelist.append(tempfilename)
		open_tempfile = open(tempfilename, 'w')
		for orgindex in range(len(MLSA_list)):
			SeqIO.write([record_dict[MLSA_list[orgindex][0]][MLSA_list[orgindex][1][locusindex]]], open_tempfile, "fasta")
		open_tempfile.close()
		
	return unaligned_filelist

def make_alignments(unaligned_filelist): 
	mylogger.debug("run_multiprocess_alignment(%s, unaligned_filelist)" % args.alignmeth)
	mylogger.info("\n%s\nAligning Orthologeous Groups (OGs) using %s and %d cpus" %(hline, args.alignmeth, args.nthreads))
	full_thread_mp_groups = len(unaligned_filelist) // args.nthreads
	remaining_mp_group_threads = len(unaligned_filelist) % args.nthreads
	aligned_filelist = []
	startindex = 0
	
	for mp_group in range(full_thread_mp_groups):
		endindex = startindex + args.nthreads
		mylogger.debug("aligned_filelist.extend(run_multiprocess_alignment(%s, unaligned_filelist[%d:%d]))" %(args.alignmeth, startindex, endindex))
		aligned_filelist.extend(run_multiprocess_alignment(args.alignmeth, unaligned_filelist[startindex:endindex]))
		sys.stdout.write("\raligned %d of %d OGs using %s" %(endindex, len(unaligned_filelist), args.alignmeth))
		sys.stdout.flush()
		startindex = endindex
	
	#finish off any remaining alignment jobs (in case the total number of alignment-jobs was not evenly divisible by the number of cpus)
	if remaining_mp_group_threads > 0:
		endindex = len(unaligned_filelist)
		aligned_filelist.extend(run_multiprocess_alignment(args.alignmeth, unaligned_filelist[startindex:endindex]))
		sys.stdout.write("\raligned %d of %d OGs using %s" %(endindex, len(unaligned_filelist), args.alignmeth))
		sys.stdout.flush()
	
#	print "\n".join(aligned_filelist)
#	print hline
	aligned_filelist.sort() #filist gets all jumbled up by multiprocessing --> sort back into original order according to OG number!
#	print "\n".join(aligned_filelist)
	return aligned_filelist

def clustalw(inputfile, mp_output):
		outputfile = inputfile.replace("unaligned_temp_fasta_", "SINGLEalignment_CLUSTALW_temp_fasta_", 1)
		clustalw_cline = ClustalwCommandline(os.path.join(args.aligner_path, args.alignmeth), INFILE = inputfile, outfile = outputfile, type = "PROTEIN", align = True, quiet = True, OUTORDER = "INPUT")
		try:
			clustalw_cline()
			mp_output.put(outputfile)
		except Exception:
			raise RuntimeError("Your clustalw version is older than v2 (probably v1.83). You should use version 2 or newer (called clustalw2 on many systems)")

def clustalw2(inputfile, mp_output):
		outputfile = inputfile.replace("unaligned_temp_fasta_", "SINGLEalignment_CLUSTALW2_temp_fasta_", 1)
		clustalw_cline = ClustalwCommandline(os.path.join(args.aligner_path, args.alignmeth), INFILE = inputfile, outfile = outputfile, type = "PROTEIN", align = True, quiet = True, OUTORDER = "INPUT")
		clustalw_cline()
		mp_output.put(outputfile)

def clustalo(inputfile, mp_output):#Todo: find out a way to check clustalo version
		outputfile = inputfile.replace("unaligned_temp_fasta_", "SINGLEalignment_CLUSTALO_temp_fasta_", 1)
		clustalomega_cline = ClustalOmegaCommandline(os.path.join(args.aligner_path, args.alignmeth), infile = inputfile, outfile = outputfile, seqtype = "Protein", threads = args.nthreads, verbose = False, force = True, outputorder = "input-order")
		clustalomega_cline()
		mp_output.put(outputfile)

def muscle(inputfile, mp_output):
		mylogger.debug("muscle(%s)" % inputfile)
		outputfile = inputfile.replace("unaligned_temp_fasta_", "SINGLEalignment_MUSCLE_temp_fasta_", 1)
		muscle_cline = MuscleCommandline(os.path.join(args.aligner_path, args.alignmeth), input = inputfile, out = outputfile, quiet = True) #add 'stable = True' to the end of this list, if the stable-bug in muscle is fixed (remove the correct_for_muscle_bug() method in that case)
		muscle_cline()
		mp_output.put(outputfile)

def correct_for_muscle_bug(aligned_filelist, seq_filelist):
	#this def is a necessary workaround for the missing "stable" function of "muscle", which has been deactivated due to the discovery of a bug in this function
	#without that function, the sequence order in the alignment might change, and that would be catastrophic for MLSA-analyses!
	#This workaround method could become obsolete, when future muscle versions reimplement a correct '-stable' option (hopefully in version 3.9)
	#In that case, add the option", stable = True" to the MuscleCommandline-call in 'call_muscle()'.
	mylogger.debug("correct_for_muscle_bug(aligned_filelist, seq_filelist)")
	mylogger.info("\n%s\ncorrecting sequence order in muscle alignments" % ("-" * 50))
	
	for f in range(len(aligned_filelist)):
		sys.stdout.write("\rcorrecting alignmentfile %d of %d" %(f + 1, len(aligned_filelist)))
		alignmentfile = aligned_filelist[f]
		alignmenthandle = open(alignmentfile, 'r')
		alignment = AlignIO.read(alignmenthandle, "fasta")
		alignmenthandle.close()
		seqs = []
		seq_file = seq_filelist[f]
		seq_handle = open(seq_file, 'r')
		
		for record in SeqIO.parse(seq_handle, "fasta"):
			seqs.append(record)
		seq_handle.close()
		corrected = None
		
		for x in range(len(seqs)):
			for y in range(len(alignment)):
				#print str(x) + " " + seqs[x].id + " == " + alignment[y].id
				if seqs[x].id == alignment[y].id:
					#print "True"
					if x == 0:
						corrected = alignment[y:y + 1]
					else:
						corrected.append(alignment[y])
				#else:
				#	print "False"
				
		alignmenthandle = open(alignmentfile, 'w')
		AlignIO.write([corrected], alignmenthandle, "fasta")
		alignmenthandle.close()
		
	return corrected

def run_multiprocess_alignment(targetfunction, unaligned_files_portion):
	mylogger.debug("run_multiprocess_alignment(%s, %s)" %(targetfunction, unaligned_files_portion))
	if __name__ == '__main__': #just making sure function is only called within its intended context
		group_results = []
		if len(unaligned_files_portion) > args.nthreads: #just an additional safety catch
			raise RuntimeError("More alignment jobs sent to queue than specified cpus! there must be something wrong with this script!")
			
		mp_output = multiprocessing.Queue()
		processes = [multiprocessing.Process(target=eval(targetfunction), args=(x, mp_output)) for x in unaligned_files_portion]
		
		for p in processes:
			p.start()
			
		for p in processes:
			p.join()
			
		group_results = [mp_output.get() for p in processes]
		return group_results
		
	else:
		raise RuntimeError("FORBIDDEN TO CALL THIS (MULTIPROCESSING) FUNCTION FROM AN EXTERNAL MODULE\n-->ABORTING")

def write_temp_alignments(alignmentlist, prefix):
	mylogger.debug("write_temp_alignments(alignmentlist, %s)" % prefix)
	proc_aligned_filelist = []
	OG_counter = 0 #numbering of the used 'Orthologeous Groups'(OGs), to differentiate between the different temp_files

	for index in range(len(alignmentlist)):
		OG_counter += 1
		alignfilename = os.path.join(args.temp_path, prefix + str(OG_counter).zfill(5) + ".fasta")
		proc_aligned_filelist.append(alignfilename)
		open_tempfile = open(alignfilename, 'w')
		AlignIO.write([alignmentlist[index]], open_tempfile, "fasta")
		open_tempfile.close()
		
	return proc_aligned_filelist

def read_alignments(input_filelist):
	mylogger.info("read_alignments(input_filelist)")
	alignmentlist = []
	if args.alignmeth == "muscle" or args.alignmeth == "clustalo" or len(input_filelist) == 1: #Last condition assumes that ONLY final result files would be passed as a filelist of only one file
		aformat = "fasta"# just a workaround for this method
	else:
		aformat = "clustal"
		
	for alignmentfile in input_filelist:
		currentalignfile = open(alignmentfile, 'r')
		alignmentlist.append(AlignIO.read(currentalignfile, aformat))
		currentalignfile.close()
		
	if len(alignmentlist) != len(input_filelist):
		raise OSError("could not read in all temporary alignmentfiles!")
		
	return alignmentlist

def remove_gaps_from_alignment_borders(alignmentlist): #optional. Better to use Gblocks if available
	mylogger.debug("remove_gaps_from_alignment_borders(alignmentlist)")
	global outputfilename
	global proc_aln_length
	outputfilename = outputfilename.replace("concatenated_orthologs_", "concatenated_orthologs_FlankingBordersDegapped_")
	processed_alignmentlist = []
	
	for currentalignment in alignmentlist:
		start, stop = 0, len(currentalignment[0])
		for c in range(stop): #identify start of alignable region
			if not "-" in currentalignment[:, c]:
				start = c
				break
		for c in range(1, stop): #identify end of alignable region
			nc = c*(-1)
			if not "-" in currentalignment[:, nc]:
				stop += nc # = maxlength minus endposition
				break
		processed_alignmentlist.append(currentalignment[:, start:stop + 1])
		
	#print "keep_temp: " + str(args.keep_temp)
	if args.keep_temp == True:
		#print "keep_temp is 'True' --> will produce prossessed single alignment files"
		write_temp_alignments(processed_alignmentlist, "Processed_removedFLANKINGgaps_SINGLEalignment_MUSCLE_temp_fasta_OG")
		
	return processed_alignmentlist

def remove_gaps_from_complete_alignments(alignmentlist): #optional. Better to use Gblocks if available
	mylogger.debug("remove_gaps_from_complete_alignments(alignmentlist)")
	global outputfilename
	global proc_aln_length
	outputfilename = outputfilename.replace("concatenated_orthologs_", "concatenated_orthologs_AllDegapped_OG")
	processed_alignmentlist = []
	counter = 0
	maxnumber = len(alignmentlist)
	
	for currentalignment in alignmentlist:
		counter += 1
		sys.stdout.write("\rprocessing alignment %d of %d" %(counter, maxnumber))
		sys.stdout.flush()
		processed_alignment = currentalignment[:, 0:0] #initialize alignment without sequence
		for c in range(len(currentalignment[0])):
			if not "-" in currentalignment[:, c]:
				processed_alignment += currentalignment[:, c:(c + 1)]
		processed_alignmentlist.append(processed_alignment)
	#print "keep_temp: " + str(args.keep_temp)
	
	if args.keep_temp == True:
		#print "keep_temp is 'True' --> will produce prossessed single alignment files"
		write_temp_alignments(processed_alignmentlist, "Processed_removedALLgaps_SINGLEalignment_MUSCLE_temp_fasta_OG")
		
	return processed_alignmentlist

def rename_for_gblocks(align_file):
	mylogger.debug("rename_for_gblocks(%s)" % align_file)
	#temporarily rename sequences, so that gblocks doesn't freak out
	mylogger.info("renaming sequences before gblocks")
	alignment = read_alignments([align_file])[0]
	index = 1
	temp_name_dict = {}
	tempfile_name = "%s_deltemp_indexed_alignments_%s" %(align_file, time.strftime("%Y%m%d%H%M%S"))
	
	for al in alignment:
		mylogger.info("renaming %s to %d" %(al.id, index))
		temp_name_dict[str(index)] = al.id
		al.id = str(index)
		al.description = ""
		index += 1
		
	write_final_alignment(tempfile_name, alignment)
	
	return tempfile_name, temp_name_dict

def call_Gblocks(file_name): #this calls Gblocks with standard settings. I tried not to overload the argument list for this python script
	mylogger.debug("call_Gblocks(%s)" % file_name)
	tempfile_name, temp_name_dict = rename_for_gblocks(file_name)
	gblocks_args = ['-t=p', '-e=-gb', '-d=n']
	
	gblocks_command = [os.path.join(gblocks_path, "Gblocks"), tempfile_name] + gblocks_args
	call(gblocks_command)
	rename_after_gblocks(tempfile_name + "-gb", temp_name_dict, file_name + "-gb")
	
	os.remove(tempfile_name)
	os.remove(tempfile_name + "-gb")

	return file_name + "-gb"

def rename_after_gblocks(align_file, temp_name_dict, final_filename):
	mylogger.debug("rename_after_gblocks(%s, temp_name_dict, %s)" %(align_file, final_filename))
	mylogger.info("renaming sequences back after gblocks")
	alignment = read_alignments([align_file])[0]
	
	for al in alignment:
		#print "renaming " + al.id + " to " + temp_name_dict[al.id]
		al.id = temp_name_dict[al.id]
		al.description = "gblocks-filtered concatenated coregenome"
		
	write_final_alignment(final_filename, alignment)

def concatenate_alignments(alignment_list, headers):
	mylogger.debug("concatenate_alignments(alignment_list, headers)")
	global aln_length
	suffix_list = ["_cds_aa.fasta", ".fasta", ".fas", ".faa", ".fa"]#list of most probable input-sequence suffixes for removal from sequence identifier
	counter = 0
	maxnumber = len(alignment_list)
	concatenated_alignment = alignment_list[0][:, 0:0]
	
	for alignment in alignment_list:
		counter += 1
		sys.stdout.write("\rconcatenating alignment %d of %d" %(counter, maxnumber))
		sys.stdout.flush()
		concatenated_alignment += alignment
	del alignment_list #not needed anymore --> free up memory
	
	if len(concatenated_alignment) == len(headers): #len(concatenated_alignment) counts the NUMBER of comparison sequences, not their length
		for index in range(len(headers)): #renaming the concatenated sequences (otherwise they will just be named after the first ortholog)
			newid = headers[index]
			for suffix in suffix_list: #this removes the most probable suffixes from the sequence names in the final alignments
					if headers[index].endswith(suffix):
						newid = headers[index].rstrip(suffix)
						break
			concatenated_alignment[index].id = newid
			concatenated_alignment[index].description = "concatenated_unique_coregenome"
	else:
		raise RuntimeError("Records in concatenated_alignment and headers are not synchronized!")
		
	if len(concatenated_alignment) != 0:
		 aln_length = len(concatenated_alignment[0])
		
	return concatenated_alignment

def write_final_alignment(outputfilename, concatenated_alignment):
	mylogger.debug("write_final_alignment(%s, concatenated_alignment)" % outputfilename)
	outputfile = open(outputfilename, 'w')
	AlignIO.write(concatenated_alignment, outputfilename, "fasta")
	outputfile.close()

def raxml_rapidbs(alignmentfile): #parameters should be a dictionary (This dictionary thing was introduced, so that the script can be more easily adapted to accept custom commandline-parameters for raxml by the user)
	mylogger.debug("raxml_rapidbs(%s)" % alignmentfile)
	mylogger.info("Calculating phylogenies: 'rapid bootstrap analyses and search for best-scoring ML Tree in one run' using raxmlHPC and " + str(args.nthreads) + " threads")
	
	outname = "MLSA_rapidBS%d_%s_final_tree" %(args.nr_bootstraps, time.strftime("%Y%m%d%H%M%S"))
	raxml_cline = RaxmlCommandline(raxml_prog, sequences = alignmentfile, algorithm = "a", model = "PROTGAMMAAUTO", name = outname, parsimony_seed = args.seed_nr, rapid_bootstrap_seed = args.seed_nr, num_replicates = args.nr_bootstraps, threads = args.nthreads)
	mylogger.info("-->" + str(raxml_cline))
	raxml_cline()
	mylogger.info("-->SUCCESS")
	
	#the resultfiles will be: "RAxML_bipartitions.rapidBS_final_tree" and "RAxML_bipartitionsBranchLabels.rapidBS_final_tree"
	#Labels on nodes or branches, respectively
	outputfiles = ["RAxML_bipartitions." + outname, "RAxML_bipartitionsBranchLabels." + outname]

def raxml_bs(alignmentfile):
	mylogger.debug("raxml_bs(%s)" % alignmentfile)
	mylogger.info("Calculating phylogenies: Thorough bootstrap analyses with raxml")
	mylogger.info("\tDetermining best ML tree of 20 raxmlHPC runs using %d threads" % args.nthreads)
	
	raxml_cline = RaxmlCommandline(raxml_prog, model = "PROTGAMMAAUTO", name = "best_delme_tempfile", parsimony_seed = args.seed_nr, num_replicates = 20, sequences = alignmentfile, threads = args.nthreads)
	mylogger.info("\t-->" + str(raxml_cline))
	raxml_cline()
	#the resultfile will be :"RAxML_bestTree.best_delme_tempfile"
	mylogger.info("\t-->SUCCESS")
	
	mylogger.info("\tDoing bootstrap analyses with %d runs using raxmlHPC using %d threads" %(args.nr_bootstraps, args.nthreads))
	
	raxml_cline = RaxmlCommandline(raxml_prog, model = "PROTGAMMAAUTO", sequences = alignmentfile, name = "boot_delme_tempfile", parsimony_seed = args.seed_nr, bootstrap_seed = args.seed_nr, num_replicates = args.nr_bootstraps, threads = args.nthreads)
	mylogger.info("\t-->" + str(raxml_cline))
	raxml_cline()
	#the resultfile will be: "RAxML_bootstrap.boot_delme_tempfile"
	mylogger.info("\t-->SUCCESS")
	
	mylogger.info("\tDrawing bipartitions of bootstrap trees onto best ML tree using raxmlHPC using " + str(args.nthreads) + " threads")
	
	outname = "MLSA_raxmlBS" + str(args.nr_bootstraps) + "_" + time.strftime("%Y%m%d%H%M%S") + "_" + "final_tree"
	raxml_cline = RaxmlCommandline(raxml_prog, model = "PROTGAMMAAUTO", parsimony_seed = args.seed_nr, algorithm = "b", starting_tree = "RAxML_bestTree.best_delme_tempfile", bipartition_filename = "RAxML_bootstrap.boot_delme_tempfile", name = outname)
	mylogger.info("\t-->" + str(raxml_cline))
	raxml_cline()
	
	#The resultfiles will be: RAxML_bipartitions.final_tree" and "RAxML_bipartitionsBranchLabels.final_tree"
	outputfiles = ["RAxML_bipartitions." + outname, "RAxML_bipartitionsBranchLabels." + outname]
	mylogger.info("\t-->SUCCESS")
	return outputfiles

def raxml(alignmentfile):
	mylogger.debug("raxml(%s)" % alignmentfile)
	outname = "MLSA_raxml_%s_final_tree" % time.strftime("%Y%m%d%H%M%S")
	mylogger.info("Calculating phylogeny: Determining best ML tree of 20 raxmlHPC runs using %d threads" % args.nthreads)
#	raxml_cline = RaxmlCommandline(raxml_prog, sequences = alignmentfile, model = "PROTGAMMAAUTO", name = outname,  parsimony_seed = args.seed_nr, num_replicates = 20, threads = args.nthreads)
	raxml_cline = RaxmlCommandline(raxml_prog, sequences = alignmentfile, model = "PROTGAMMAAUTO", name = outname,  parsimony_seed = args.seed_nr, num_replicates = 2, threads = args.nthreads)
	mylogger.info("\t-->" + str(raxml_cline))
	raxml_cline()
	#the resultfile will be :"RAxML_bestTree.final_tree"
	outputfiles = ["RAxML_bestTree." + outname]
	mylogger.info("\tSUCCESS")
	print "deleting temporary files"
	for delfile in os.listdir("."):
		if delfile.startswith("RAxML_") and "." + outname + ".RUN." in delfile:
			#print "deleting temp-file: " + delfile
			os.remove(delfile)
	return outputfiles

def randomnumber():
	mylogger.debug("randomnumber()")
	#returns a random integer to use as seed for rxml and pyml
	random.seed()
	
	return random.randint(1, 2000)

def print_version():
	print "\n ==PO_2_MLSA.py %s by John Vollmers==\n" % version


def main():
	aligned_filelist, unaligned_filelist = [], []
	mylogger.debug("main()")
	global args, docontinue
	
	
	try:
		checkargs()
		##move this to checkargs()
		infotext = "-keep temporary files: %s" % args.keep_temp
		if not args.keep_temp:
			infotext += " --> will delete all temporary files "
		mylogger.info(infotext)
		mylogger.info("-using aligner '%s'" % args.alignmeth)
		mylogger.info("-alignment filter: " + args.afilter)
		headers, MLSA_list, OG_number = read_PO_file(args.PO_file)
		##move this to checkargs()
		
		#read input
		record_dict = read_fasta_seqs(headers, MLSA_list)
		seq_filelist = write_seq_files(record_dict, MLSA_list)
		unaligned_filelist = write_temp_files(record_dict, MLSA_list, "unaligned_temp_fasta_OG")
		
		#alignments
		aligned_filelist = make_alignments(unaligned_filelist)
		if args.alignmeth == "muscle":
			correct_for_muscle_bug(aligned_filelist, unaligned_filelist) #Necessary because bug in muscle '-stable' option (option disabled for this reason as of muscle version 3.8)
		alignment_list = read_alignments(aligned_filelist)
		
		#concatenate and filter alignments
		#custom filter (not recommended)
		if args.afilter == "degap_all":
			mylogger.info("\n%s\nRemoving all gapped positions from all single alignments prior to concatenation", hline)
			alignment_list = remove_gaps_from_complete_alignments(alignment_list)
		elif args.afilter == "flanking":
			mylogger.info("\n%s\nRemoving only the flanking gapped positions from all single alignments prior to concatenation" % hline)
			alignment_list = remove_gaps_from_alignment_borders(alignment_list)
		elif not gblocks and verbose:
			mylogger.info("Leaving the single alignments as they are (Not removing any gapped or unconserved positions)")
		mylogger.info("\n%s\nconcatenating alignments" % hline)
		concatenated_alignment = concatenate_alignments(alignment_list, headers)
		mylogger.info("\n%s\nwriting concatenated alignment to fasta-file: %s" %(hline, outputfilename))
		write_final_alignment(outputfilename, concatenated_alignment)
		#gblocks filter(highly recommended)
		if gblocks:
			mylogger.info("running gblocks on %s" % outputfilename)
			call_Gblocks(outputfilename)
			final_check = read_alignments([outputfilename + "-gb"])
			proc_aln_length = len(final_check[0][0])
			mylogger.info("-->Processed File: %s-gb\n-->Gblocks-Logfile: %s-gb.htm" %(outputfilename, outputfilename))
		
		#give summary
		mylogger.info("=" * 50)
#		mylogger.info("Finished alignments!\nthe individual (unaligned) ortholog selections were stored as the following multifastas: ")
#		mylogger.info(", ".join(seq_filelist))
		mylogger.info("Finished alignments!\nthe individual (unaligned) ortholog selections were stored under the following basename:\n\t%s " % os.path.join(args.out_path, "OGselection_MLSA_single_unaligned_multifasta_*"))
		mylogger.info("\n%d genes concatenated to a final sequence length of %d residues per organism" %(OG_number, aln_length))
		if gblocks:
			mylogger.info("after 'gblocks' the remaining sequence length is %s residues" % proc_aln_length)
			mylogger.info("\n-->final aligned and processed sequence saved as: %s-gb" % outputfilename)
		else:
			mylogger.info("\n-->final aligned sequence saved as: %s" % outputfilename)
			mylogger.info("\nRemember to set sequence-type as 'PROTEIN' and filetype as 'fasta_wgap' when loading the MLSA sequence into Arb!\n")
			
		#RaxML
		if "raxml" in args.tree_method:
			mylogger.info("\nNow generating maximum likelihood trees from final output alignment")
			if args.seed_nr == 0:
				args.seed_nr = randomnumber()
			treefiles = []
			outtreename = "%s.%s.TREE" %(outputfilename, args.tree_method)
			treefiles.extend(eval(args.tree_method)(outputfilename))
			if docontinue:
				mylogger.info("the following trees were generated:")
				mylogger.info("\t" + "\n\t".join(treefiles))
			else:
				mylogger.info("some kind of error occured during ML calculation")
		elif "nj" in args.tree_method:
			mylogger.info("\nNow generating neighbor-joining tree from final output alignment")
			mylogger.debug("reading final alignment file: %s-gb" %outputfilename)
			final_align = AlignIO.read(outputfilename + "-gb", "fasta") 
			njtree = make_single_njtree(final_align)
			if args.tree_method == "nj_bs" and args.nr_bootstraps > 1:
				bs_list = run_multiprocess_bootstrap(final_align) #generate bootstrap trees using multiprocessing
				final_tree = get_support(njtree, bs_list) #get support values from bootstrap trees
			else:
				final_tree = njtree
			nj_filename = write_nj_tree(final_tree, outputfilename)
			mylogger.info("wrote neighbor joining tree in newick format to %s" % nj_filename)
			
	except Exception as e:
		for frame in traceback.extract_tb(sys.exc_info()[2]):
			fname,lineno,fn,text = frame
			#print "Error in %s on line %d :> %text" % (fname, lineno, text)
		mylogger.error("Error in %s, %s, on line %d :> %s" % (fname, fn, lineno, text))
		mylogger.error(str(e))
		
	finally:
		#move cleanups here!
		if len(unaligned_filelist) > 0 and len(aligned_filelist) > 0 and not args.keep_temp:
			mylogger.info("\n%s\ncleaning up temporary alignmentfiles" % hline)
			for delfile in unaligned_filelist:
				os.remove(delfile)
				if "clustalw" in args.alignmeth:
					os.remove(delfile.rstrip(".fasta") + ".dnd") #remove pesky "guide-tree" files produced by clustal aligners as well
			for delfile in aligned_filelist:
				os.remove(delfile)
		mylogger.info("cleaning up RaxML tempfiles")
		for delfile in os.listdir("."):
			if "delme_tempfile" in delfile:
				os.remove(delfile)

def make_single_njtree(protalign):
	mylogger.debug("make_single_njtree(protalign)")
	#subst_model can be ["identity", "BLOSUM", "PAM"] Default= "identity", "BLOSUM"= "BLOSUM62", "PAM"="PAM160"
#	try:
#		from Bio.Phylo.TreeConstruction import DistanceCalculator
#		from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
#	except ImportError:
#		mylogger.error("failed to import Bio.Phylo.TreeConstruction methods DistanceCalculator and DistanceTreeConstructor!\nPlease make sure you are using the newest version of BioPython")
	mylogger.debug("DistanceCalculator(%s)" % args.subst_model)
	distcalc = DistanceCalculator(args.subst_model)
	mylogger.debug("DistanceTreeConstructor(distcalc, 'nj')")
	tree_constructor = DistanceTreeConstructor(distcalc, 'nj')
	mylogger.debug("tree_constructor.build_tree(protalign)")
	tree = tree_constructor.build_tree(protalign)
	mylogger.debug("finished with nj_tree")
	return tree

#todo: phylip method (perhaps faster?)
#  convert alignment into phylip format (rename seqs first)
#  write phylip alignment to file
#  call emboss seqboot to create permutations (multiprocessing)
#  call emboss protdist to create distance-matrices for all permutations
#  call emboss fneighbor to create neighbor-joining trees for all distance_matrices
#  also call emboss protdist and fneighbor on original alignment file to create main tree
#  read in main_tree
#  read in support-trees
#  read in fuck
#  read in all trees from fneighbor results and in frer bootstrap support

def calculate_phylip_NJ_trees(temp_matrix_files):
	mylogger.debug("calculate_phylip_NJ_trees(temp_matrix_files)")
	if len(temp_matrix_files) == 0:
		mylogger.warning("No Matrix files found in list!")
	from Bio.Emboss.Applications import FNeighborCommandline
	NJtrees = []
	for mf in temp_matrix_files:
		treefile = mf + ".tree"
		cline = FNeighborCommandline(datafile = mf,\
									outtreefile = treefile,\
									jumble = True,\
									outfile = "delfile_%s" % timetag) #set outfile = /dev/null to minimize unneccesary output
		cline()
		NJtrees.append(treefile)
		if not args.no_verbose:
			sys.stdout.write("\rGenerated tree %s of %s" %(len(NJtrees), len(temp_matrix_files)))
			sys.stdout.flush()
#	temp_files.append("delfile_%s" % timetag) #delete this line when setting outfile = /dev/null
	mylogger.info("finished generating NJ trees using Phylip")
	return NJtrees

#def make_bootstrap_njtrees(protalign, nr_bootstraps, mp_output):
#	mylogger.debug("make_bootstrap_njtrees(protalign, %s)" % nr_bootstraps)
#	#subst_model can be ["identity", "BLOSUM", "PAM"] Default= "identity", "BLOSUM"= "BLOSUM62", "PAM"="PAM160"
#	try:
#		#bs_protalign = bootstrap(protalign, nr_bootstraps) 
#		distcalc = DistanceCalculator(args.subst_model)
#		tree_constructor = DistanceTreeConstructor(distcalc, 'nj')
#		trees = bootstrap_trees(protalign, nr_bootstraps, tree_constructor)
#		#return list(trees) #bootstrap_trees() returns a generator object that has to be converted
#		mylogger.debug("another set of bootstrap trees successfully generated,hu")
#		mp_output.put(list(trees))
#	except:
#		mylogger.error("Something went wrong during bootstrap tree generation")
#		mp_output.put([None])
		
def make_bootstrap_njtrees(protalign, nr_bootstraps): #obsolete soon
	mylogger.debug("make_bootstrap_njtrees(protalign, %s)" % nr_bootstraps)
	#subst_model can be ["identity", "BLOSUM", "PAM"] Default= "identity", "BLOSUM"= "BLOSUM62", "PAM"="PAM160"
	try:
		#bs_protalign = bootstrap(protalign, nr_bootstraps) 
		distcalc = DistanceCalculator(args.subst_model)
		tree_constructor = DistanceTreeConstructor(distcalc, 'nj')
		trees = bootstrap_trees(protalign, nr_bootstraps, tree_constructor)
		#return list(trees) #bootstrap_trees() returns a generator object that has to be converted
		mylogger.debug("another set of bootstrap trees successfully generated,hu")
		treelist=[]
		counter = 0
		for t in trees:
			counter += 1
			print "appending bs-tree nr %s" % counter
			treelist.append(t)
		print "finished appending %s trees" % counter
		return treelist
	except:
		mylogger.error("Something went wrong during bootstrap tree generation")
		return [None]

def get_bootstrap_support(maintree, bs_trees):
	mylogger.debug("get_bootstrap_support(maintree, bs_trees)")
#	try:
#		from Bio import Phylo
#		from Bio.Phylo.Consensus import get_support
#	except ImportError:
#		mylogger.Error("get_bootstrap_support() failed. Biopython modules \"Phylo\" and \"Phylo.consensus\"could not be imported\n\
#		Is the newest version of BioPython installed?")
#		raise ImportError
	support_tree = get_support(maintree, bs_trees)
	return support_tree

def write_nj_tree(tree, outname):
	mylogger.debug("write_nj_tree(%s)" % tree)
	
	from Bio import Phylo
	outname += ".%s" %(args.tree_method)
	if args.tree_method == "nj_bs":
		outname += ".%s" % args.nr_bootstraps
	outname += ".tree.newick"
	mylogger.info("writing treefile to %s" % outname)
	
	out_file = open(outname, "w")
	Phylo.write(tree, out_file, "newick")
	out_file.close()
	return outname

def write_distmatrix(distmatrix):
	mylogger.debug("write_distmatrix(distmatrix)")
	#optional: write distance matrix to file
	#figure out how to do that later
	
#def run_multiprocess_bootstrap(alignments):
#	#lazily copied and adapted the other mulitiprocess function for this purpose instead of adapting one function to serve both purposes
#	mylogger.debug("run_multiprocess_bootstrap(alignments)")
#	if __name__ == '__main__': #just making sure function is only called within its intended context
#		group_results = []
#	#define portions for multiprocessing
#		if args.nr_bootstraps % args.nthreads == 0:
#			bs_portions = [ args.nr_bootstraps // args.nthreads for i in range(args.nthreads)]
#		else:
#			bs_portions = [ args.nr_bootstraps // (args.nthreads - 1) for i in range(args.nthreads - 1)]
#			remainder = args.nr_bootstraps % (args.nthreads - 1)
#			bs_portions.append(remainder)
#
#			if len(bs_portions) > args.nthreads: #just an additional safety catch
#				raise RuntimeError("More alignment jobs sent to queue than specified cpus! there must be something wrong with this script!")
#			
#			#start bootstrapping
#			mp_output = multiprocessing.Queue()
#			processes = [multiprocessing.Process(target=make_bootstrap_njtrees, args=(alignments, x, mp_output)) for x in bs_portions]
#			mylogger.debug("starting %s processes")
#			for p in processes:
#				p.start()
#			mylogger.debug("all processes running")
#			mylogger.debug("now beginning to join %s processes" % len(processes))
#			print processes
#			for p in processes:
#				print p
#				p.join()
#				print "joined one!"
#			mylogger.debug(" all processed joined")
#			group_results = [mp_output.get() for p in processes]
#			
#			#now combine these results
#			bootstrap_treelist = []
#			for gr in group_results:
#				bootstrap_treelist.append(gr)
#			mylogger.debug("generated %s bootstrap trees in total" % len(bootstrap_treelist))
#			return bootstrap_treelist
#	else:
#		raise RuntimeError("FORBIDDEN TO CALL THIS (MULTIPROCESSING) FUNCTION FROM AN EXTERNAL MODULE\n-->ABORTING")

def run_multiprocess_bootstrap(alignments):
	mylogger.debug("currently ignoring multiprocessing for bootstrap_generation!\nthis may take a while")
	mylogger.debug("generating bootstraptrees")
	bs_trees = make_bootstrap_njtrees(alignments, args.nr_bootstraps)
	mylogger.debug("finished generating %s bootstraptrees" % len(bs_trees))
	return bs_trees

if args.showversion:
	print_version()
else:
	mylogger = setlogger(logfilename)
	main()
