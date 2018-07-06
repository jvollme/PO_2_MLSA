#!/usr/bin/env python
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
import numpy

#TODO: add parameters for gblocks
#base options
myparser=argparse.ArgumentParser(description="\n==PO_2_MLSA.py v1.5 by John Vollmers==\nCreates concatenated alignments of UNIQUE Genes with orthologues in comparison organisms for the construction of MLSA-based phylogenetic trees.\nOptionally the resulting concatenated alignments may contain all gapped alignmentpositions or may be stripped either of ALL gapped positions or of all gapped positions in the flanking regions of each composite ortholog\nThis script is supposed to be part of a pipeline consisting of:\n\tA.)Conversion of Genbank/Embl-Files to ANNOTATED(!) Fastas using CDS_extractor.pl by Andreas Leimbach\n\tB.)Calculation of orthologs and paralogs using proteinortho5 (WITH the '-single' and '-self' arguments!)\n\tC.)The creation of concatenated MLSA-sequences based on:\n\t\t-the fasta sequences of step A\n\t\t-the proteinortho5-results from step B\n\nThe output-file will be in fasta format (but including gapped positions, so remember to use 'fasta_wgap' when loading into Arb!). However it's absolutely no problem to include other common output-alignmentformats on request!", formatter_class=argparse.RawTextHelpFormatter)
myparser.add_argument("-po", "--proteinortho", action = "store", dest = "PO_file", help = "(String) file with proteinortho5 results", required = True)
myparser.add_argument("-f", "--fastas", action = "store", dest = "fasta_path", help = "(String) path to fasta files (produced by CDS_extractor) \nDefault = current working directory", default = "")
myparser.add_argument("-t", "--threads", action = "store", dest = "nthreads", type = int, default = 1, help = "Maximum number of threads to use for alignment steps\nDefault = 1")
myparser.add_argument("--filter", action = "store_true", dest = "filter_outliers", default = False, help = "OPTIONAL: Filter outliers based on single marker phylogenies in order to remove HGT events (Default: False)")

#output-options
myparser.add_argument("-op", "--out_path", action = "store", dest = "out_path", default = ".", help = "Path to output (will be created if it does not exist)\nDefault = current working directory")
myparser.add_argument("-tp", "--temp_path", action = "store", dest = "temp_path", default = "", help = "Path for temporary files (will be created if does not exist)\nDefault =  current working directory")
myparser.add_argument("-kt", "--keep_temp", action = "store_true", dest = "keep_temp", default = False, help = "Keep all temporary files and intermediate results\n(Multifasta_files containing the singe genes involved in the MLSA-alignments are stored in any case)\nDefault: delete all temporary files")
#myparser.add_argument("-am", "--align_method", action = "store", dest = "alignmeth", choices = ["muscle", "clustalw", "clustalw2", "clustalo"], default = "muscle", help = "The Alignmentmethod to use.\nDefault = 'muscle'")
#myparser.add_argument("-dg", "--degap", action = "store", dest = "degap", choices = ["none", "all", "flanking"], default = "all", help = "(Only meant for use, if gblocks is not installed\nSpecify if and which gaps to remove:\n\t'none' keep all gapped positions in the final MLSA-alignment\n\t'all' remove ALL gapped postitions in the alignments\n\t'flanking' remove flanking gapped positions from all individual alignments\nDefault = 'all'")

#external binary options
myparser.add_argument("-mb", "--muscle_binary", action = "store", dest = "muscle_binary", default = "muscle", help = "Muscle binary to use (including path if not in $PATH). default: assume \"muscle\" in $PATH")
myparser.add_argument("-gbb", "--gblocks", action = "store", dest = "gblocks_binary", default = "Gblocks", help = "Gblocks binaries to use (include path to binaries if not listed in $PATH). Default: \"Gblocks\" (assumed to be in PATH)")
myparser.add_argument("-rmlb", "--raxml", action = "store", dest = "raxml_binary", default = "", help = "RaxML excecutable biaries to use (Include path to binary if not listed in $PATH). Default: check for common naming of raxml binaries in $PATH, prioritizing binaries with \"PTHREADS\" in the name")

#tree-calculation options
myparser.add_argument("-mt", "--make_tree", action = "store", dest = "tree_method", choices = ["raxml", "raxml_rapidbs", "nj", "nj_bs", "none"], default = "none", help = "Generate ML phylogenetic trees using RAxML with the substitution model \"PROTGAMMAAUTO\"\n\tchoices:\t\"raxml\": single tree without bootstraps (using new rapid hill climbing)\n\t\traxml_rapidbs: rapid bootstrap analyses and search for best ML tree in one run\n\t\tnone\nDefault = none")
myparser.add_argument("-sd", "--seed", action = "store", dest = "seed_nr", type = int, default = 0, help = "Integer to provide as seed for RAxML\n0 = seed generated randomly\nDefault = random seed")
myparser.add_argument("-bs", "--bootstraps", action = "store", dest = "nr_bootstraps", type = int, default = 1000, help = "Number of bootstraps(if any)\ndefault = 1000")
#myparser.add_argument("-ctba", "--custom_tree_builder_args", action = "store", dest = "custom_tree_builder_args", default = None, help = "custom arguments for raxml. CAUTION: will overide Only use if you know EXACTLY what you are doing!")
myparser.add_argument("--nj_substmodel", action = "store", dest = "subst_model",\
 choices = ['identity', "ident", 'blosum', 'pam',\
 'blosum30', 'blosum35', 'blosum40', 'blosum45', 'blosum50', 'blosum55', 'blosum60','blosum62',\
 'blosum65', 'blosum70', 'blosum75', 'blosum80', 'blosum85', 'blosum90', 'blosum95', 'blosum100',\
 'pam30','pam60', 'pam90', 'pam120', 'pam180', 'pam250',  'pam300'],\
  default = "ident", help = "Substitution model for distance matrix calculation. Only \"identity\" and variances of \"blosum\" and \"pam\" are offered. abbreviations:\n\t\"identity\"\n\t\"blosum\" = blosum62\n\t\"pam\" = pam120, default = \"identity\"")

#other options
myparser.add_argument("-v", "--version", action = "store_true", dest = "showversion", default = False, help = "show version information and then quit (don't run complete script)")
myparser.add_argument("-s", "--silent", action = "store_true", dest = "no_verbose", help = "non-verbose mode")
myparser.add_argument("--debug", action = "store_true", dest = "debug", default = False, help = "Log extra info for debugging")
#myparser.add_argument("--existing_align", action = "store", dest = existing_align, default = None, help = "contignue calculation from existing concatenated (and filtered!) alignment generated with PO_2_MLSA.py

args = myparser.parse_args()

#TOdo: add option "return_selection" to store selection of MLSA genes as unagligned multifastas or only lists of fasta-headers

version = "v1.5.2"
available_cores = multiprocessing.cpu_count() #counts how many cores are available, to check if the user-argument for threads can be fulfilled
aln_length, proc_aln_length, OG_number = 0, 0, 0
wstrings, estrings, lstrings = [], [], [] #warning, error and log messages respectively
alignmeth = "muscle"

if not os.path.exists(args.out_path):
	os.makedirs(args.out_path)
outputfilename = os.path.join(args.out_path, "concatenated_orthologs_%s_%s.fasta" %(alignmeth, os.path.basename(args.PO_file)))
logfilename = os.path.join(args.out_path, "PO_2_MLSA_%s.log" % time.strftime("%Y%m%d%H%M%S"))
raxml_prog = "raxmlHPC"
verbose = True
docontinue = True
gblocks = True
#if args.afilter != "gblocks":
#	gblocks = False
#gblocks_binary = args.gblocks_binary
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
		
	if args.muscle_binary != "muscle":
		if os.path.exists(args.muscle_binary) and os.path.isdir(args.muscle_binary):
			raise OSError("{} was specified as binary for muscle-aligner, but is actually a directory".format(args.muscle_binary))
		elif not os.path.exists(args.muscle_binary):
				raise OSError("muscle binary \"{}\" does not exist".format(args.muscle_binary))
		elif os.path.exists(args.muscle_binary) and os.path.isfile(args.muscle_binary):
			mylogger.info("using \"{}\" for alignments".format(args.muscle_binary))
	elif args.muscle_binary == "muscle":
		test_aligner = which(alignmeth)
		if test_aligner == None:
			raise OSError("Can't find \"%s\" in any directory in the PATH variable. Please provide the complete path to the binary" % alignmeth)
			
	if not os.path.exists(args.PO_file) or not os.path.isfile(args.PO_file):
		raise OSError("cannot find proteinortho-resultfile: %s" % args.PO_file)
		
	if args.temp_path != "" and docontinue:
		if not os.path.exists(args.temp_path) or not os.path.isdir(args.temp_path):
			if verbose:
				mylogger.info("Creating directory for temporary and intermediate result files: %s" % args.temp_path)
			os.mkdir(args.temp_path)
		mylogger.info("-Will store temporary and intermediate result files in %s" % os.path.abspath(args.temp_path))
		
	if args.out_path != "" and docontinue:
		if not os.path.exists(args.out_path) or not os.path.isdir(args.out_path):
			if verbose:
				mylogger.info("Creating directory for final result files: %s" % args.out_path)
			os.mkdir(args.out_path)
		mylogger.info("-Will store final result files in %s" % os.path.abspath(args.out_path))
		
	if alignmeth == "clustalo" and docontinue: #check clustalo version
		clustalomega_cline = ClustalOmegaCommandline("clustalo", version = True)
		clustalo_version = clustalomega_cline()[0].rstrip().split(".")
		if int(clustalo_version[0]) < 1 or (int(clustalo_version[0]) == 1 and int(clustalo_version[1]) < 2) : #only accept versions 1.2 and newer
			raise OSError("found clustalo version is v%s ! Version 1.2 or higher is required!" % ".".join(clustalo_version))
			
	if args.gblocks_binary in ["Gblocks", "gblocks"]:
		for gblock_name in ["Gblocks", "gblocks"]: #apparently gblocks can be named in lower OR upper case on some systems
			test_gblocks = which(gblock_name)
			if test_gblocks != None:
				args.gblocks = test_gblocks
				break
		if test_gblocks == None:
			raise OSError("can't locate gblocks in any path in PATH variable. please provide a Path to gblocks using the '-gbp' agrument, or choose a different filtering option ('-F')")
		elif verbose:
			mylogger.info("Located gblocks executable: %s" % test_gblocks)
	else:
		if os.path.exists(args.gblocks_binary) and os.path.isfile(args.gblocks_binary):
			if verbose:
				print "Located gblocks executable: %s" % os.path.join(args.gblocks_binary, "gblocks")
		elif os.path.isdir(args.gblocks_binary):
			raise OSError("\"{}\" specified as Gblocks binary, But it is a directory not a binary!".format(args.gblocks_binary))
		else:
			raise OSError("the specified gblocks executable \"{}\" does not exist!".format(args.gblocks_binary))
			
	if args.tree_method in ["raxml", "raxml_bs", "raxml_rapidbs"]:
		#print "CHECKING raxml-binaries"
		if args.raxml_binary == "":
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
		elif os.path.exists(args.raxml_binary) and os.path.isfile(args.raxml_binary):
			if args.nthreads > 1 and not "PTHREADS" in args.raxml_binary:
				mylogger.warning("multithreading is only supported with 'PTHREADS'-versions of raxml. Not sure if your choosen binaries support this.\t\nif raxml calculations fail, recompile raxml with 'PTHREADS'-option")
			try:
				checkraxml_cline = RaxmlCommandline(version = True)
				versiontext = checkraxml_cline()[0]
				startv = versiontext.find("version ")+len("version ")
				endv = versiontext[startv:].find(" ")
				version = versiontext[startv:startv+endv].split(".")
				if int(version[0]) < 8 or (int(version[0]) == 8 and int(version[1]) == 0 and int(version[2]) < 20):
						mylogger.warning("This script was devised for and tested with RAxML v8.0.20. Your version is v%s !\n\tThis may very well still work, but if it doesn't it's YOUR fault!" % ".".join(version))
				raxml_prog = args.raxml_binary
				
			except:
				mylogger.warning("Correct raxML-version not found under %s !\nWill NOT calculate ML trees!" % args.raxml_binary)
				args.tree_method = "none"

def which(thisfile):
	mylogger.debug("which(%s)", thisfile)
	for path in os.environ["PATH"].split(":"):
		if os.path.exists(path + "/" + thisfile):
			return path + "/" + thisfile
	return None

def check_PO_format(line):
			if line.startswith("#species"):
				mylogger.warning("It seems you are using results from Proteinortho4. It's recommended to use Proteinortho5 or later!\n\t(However, the derivation of MLSA genes should still work)")
			elif line.startswith("# Species"):
				if verbose:
					print "\nProteinortho-results seem to be based on Proteinortho5 or later. Good."
			else:
				mylogger.warning("Cannot clearly recognize Format of Proteinortho-results! This may produce erroneous results!\n\t For best results use Proteinortho5!")	
	
def read_PO_file(filename):
	cutoff_perc = 99 #testing request that markers present in MOST but not ALL are included
	mylogger.debug("read_PO_file(%s)" % filename)
	open_PO_file = open(filename, 'r')
	firstline = True
	org_index = 3 #default index of beginning of organism-result-columns in proteinortho4+
	headers = []
	MLSA_gene_counter = 0
	MLSA_list = []
	
	for line in open_PO_file:
		if firstline:
			check_PO_format(line)
			headers = line.split("\t")[org_index:]
			#remove line_end_symbols:
			for h in range(len(headers)):
				headers[h] = headers[h].rstrip()
			MLSA_list = [[header, []] for header in headers]
			firstline = False
			
		elif not line.startswith("#"):
			if not "," in line: #ignore line if multiple orthologs for ANY comparison organism
				if not "*" in line: #ignore line if ortholog is missing for ANY comparison organism
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
	mylogger.debug("run_multiprocess_alignment(%s, unaligned_filelist)" % alignmeth)
	mylogger.info("\n%s\nAligning Orthologeous Groups (OGs) using %s and %d cpus" %(hline, alignmeth, args.nthreads))
	full_thread_mp_groups = len(unaligned_filelist) // args.nthreads
	remaining_mp_group_threads = len(unaligned_filelist) % args.nthreads
	aligned_filelist = []
	startindex = 0
	
	for mp_group in range(full_thread_mp_groups):
		endindex = startindex + args.nthreads
		mylogger.debug("aligned_filelist.extend(run_multiprocess_alignment(%s, unaligned_filelist[%d:%d]))" %(alignmeth, startindex, endindex))
		aligned_filelist.extend(run_multiprocess_alignment(alignmeth, unaligned_filelist[startindex:endindex]))
		sys.stdout.write("\raligned %d of %d OGs using %s" %(endindex, len(unaligned_filelist), alignmeth))
		sys.stdout.flush()
		startindex = endindex
	
	#finish off any remaining alignment jobs (in case the total number of alignment-jobs was not evenly divisible by the number of cpus)
	if remaining_mp_group_threads > 0:
		endindex = len(unaligned_filelist)
		aligned_filelist.extend(run_multiprocess_alignment(alignmeth, unaligned_filelist[startindex:endindex]))
		sys.stdout.write("\raligned %d of %d OGs using %s" %(endindex, len(unaligned_filelist), alignmeth))
		sys.stdout.flush()
	
#	print "\n".join(aligned_filelist)
#	print hline
	aligned_filelist.sort() #filist gets all jumbled up by multiprocessing --> sort back into original order according to OG number!
#	print "\n".join(aligned_filelist)
	mylogger.debug("length of aligned_filelist: %s" %len(aligned_filelist)) 
	return aligned_filelist

def clustalw(inputfile, mp_output): #deactivated
		outputfile = inputfile.replace("unaligned_temp_fasta_", "SINGLEalignment_CLUSTALW_temp_fasta_", 1)
		clustalw_cline = ClustalwCommandline(os.path.join(args.aligner_path, alignmeth), INFILE = inputfile, outfile = outputfile, type = "PROTEIN", align = True, quiet = True, OUTORDER = "INPUT")
		try:
			clustalw_cline()
			mp_output.put(outputfile)
		except Exception:
			raise RuntimeError("Your clustalw version is older than v2 (probably v1.83). You should use version 2 or newer (called clustalw2 on many systems)")

def clustalw2(inputfile, mp_output): #deactivated
		outputfile = inputfile.replace("unaligned_temp_fasta_", "SINGLEalignment_CLUSTALW2_temp_fasta_", 1)
		clustalw_cline = ClustalwCommandline(os.path.join(args.aligner_path, alignmeth), INFILE = inputfile, outfile = outputfile, type = "PROTEIN", align = True, quiet = True, OUTORDER = "INPUT")
		clustalw_cline()
		mp_output.put(outputfile)

def clustalo(inputfile, mp_output):#Todo: find out a way to check clustalo version #deactivated
		outputfile = inputfile.replace("unaligned_temp_fasta_", "SINGLEalignment_CLUSTALO_temp_fasta_", 1)
		clustalomega_cline = ClustalOmegaCommandline(os.path.join(args.aligner_path, alignmeth), infile = inputfile, outfile = outputfile, seqtype = "Protein", threads = args.nthreads, verbose = False, force = True, outputorder = "input-order")
		clustalomega_cline()
		mp_output.put(outputfile)

def muscle(inputfile, mp_output):
		#mylogger.debug("muscle(%s)" % inputfile)
		outputfile = inputfile.replace("unaligned_temp_fasta_", "SINGLEalignment_MUSCLE_temp_fasta_", 1)
		muscle_cline = MuscleCommandline(args.muscle_binary, input = inputfile, out = outputfile, quiet = True) #add 'stable = True' to the end of this list, if the stable-bug in muscle is fixed (remove the correct_for_muscle_bug() method in that case)
		muscle_cline()
		mp_output.put(outputfile)

def correct_for_muscle_bug(aligned_filelist, seq_filelist):
	#this def is a necessary workaround for the missing "stable" function of "muscle", which has been deactivated due to the discovery of a bug in this function
	#without that function, the sequence order in the alignment might change, and that would be catastrophic for MLSA-analyses!
	#This workaround method could become obsolete, when future muscle versions reimplement a correct '-stable' option (hopefully in version 3.9)
	#In that case, add the option", stable = True" to the MuscleCommandline-call in 'call_muscle()'.
	mylogger.debug("correct_for_muscle_bug(aligned_filelist, seq_filelist)")
	mylogger.info("\n%s\ncorrecting sequence order in muscle alignments" % ("-" * 50))
	mylogger.debug("alignmentfiles to read %s\n" %len(aligned_filelist))
	mylogger.debug("seqeunce files to read %s\n" %len(aligned_filelist))
	#sys.stdout.write(str(aligned_filelist))
	
	for f in range(len(aligned_filelist)): #I)open each alignmentfile, II) open corresponding unaligned seqeunce file III) iterate through seqeunce file and grab all aligned seqeunces in the exact order they are listed in the seqeunce file
		alignmentfile = aligned_filelist[f]
		alignmenthandle = open(alignmentfile, 'r')
		alignment = AlignIO.read(alignmenthandle, "fasta")
		alignmenthandle.close()
		mylogger.debug("correcting alignmentfile %d of %d (%s) containg %d seqeunces" %(f + 1, len(aligned_filelist), aligned_filelist[f], len(alignment)))
		seqs = []
		seq_file = seq_filelist[f]
		seq_handle = open(seq_file, 'r')
		
		for record in SeqIO.parse(seq_handle, "fasta"):
			seqs.append(record)
		seq_handle.close()
		corrected = None# corrected will be a alignment-object (list-like) containing all comparison-organisms orthologs of ONE markergene (iterating through markergenes) in odentical order (of comparison-org)
		
		for x in range(len(seqs)):
			for y in range(len(alignment)):
				#print str(x) + " " + seqs[x].id + " == " + alignment[y].id
				if seqs[x].id == alignment[y].id:
					#print "True"
					if x == 0:
						corrected = alignment[y:y + 1]# if it is the FIRST sequence --> create new alignment object from template using slice annotation
					else:
						corrected.append(alignment[y])#otherwise just append to the existing alignment object using list-append function
				#else:
				#	print "False"
		mylogger.debug("after correction corrected file %d of %d (%s) has %d seqeunces" %(f + 1, len(aligned_filelist), aligned_filelist[f], len(corrected)))				
		alignmenthandle = open(alignmentfile, 'w')
		AlignIO.write([corrected], alignmenthandle, "fasta")
		alignmenthandle.close()
	sys.stdout.write("\ncorrecting alignmentfiles completed")
	mylogger.info("corrected %s alignmentfiles" %f)  #PROBLEMS MAY OCCUR IF THE SAME LOCUS_TAG OCCURS TWICE
	#return corrected #NOT needed! calling function does not process this anyway

def run_multiprocess_alignment(targetfunction, unaligned_files_portion):
	#mylogger.debug("run_multiprocess_alignment(%s, %s)" %(targetfunction, unaligned_files_portion))
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
	if alignmeth == "muscle" or alignmeth == "clustalo" or len(input_filelist) == 1: #Last condition assumes that ONLY final result files would be passed as a filelist of only one file
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

def remove_gaps_from_alignment_borders(alignmentlist): #optional. Better to use gblocks if available
	mylogger.debug("remove_gaps_from_alignment_borders(alignmentlist)")
	global outputfilename
	global proc_aln_length
	outputfilename = outputfilename.replace("concatenated_orthologs_", "concatenated_orthologs_FlankingBordersDegapped_Gblocked")
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
		
	if args.keep_temp == True:
		processed_alignmentfiles = write_temp_alignments(processed_alignmentlist, "Processed_removedFLANKINGgaps_SINGLEalignment_MUSCLE_temp_fasta_OG")
	else:
		processed_alignmentfiles = write_temp_alignments(processed_alignmentlist, "delme_tempfile_Processed_removedFLANKINGgaps_SINGLEalignment_MUSCLE_temp_fasta_OG")
	#return processed_alignmentlist
	return processed_alignmentfiles
	

def remove_gaps_from_complete_alignments(alignmentlist): #optional. Better to use gblocks if available
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
	
	#if args.keep_temp == True:
		#print "keep_temp is 'True' --> will produce prossessed single alignment files"
	#	write_temp_alignments(processed_alignmentlist, "Processed_removedALLgaps_SINGLEalignment_MUSCLE_temp_fasta_OG")
	write_temp_alignments(processed_alignmentlist, "Processed_removedALLgaps_SINGLEalignment_MUSCLE_temp_fasta_OG")
	return processed_alignmentlist


def rungblocks_on_single_alignments(alignmentlist, OG_number, headers):
	mylogger.debug("rungblocks_on_single_alignments(alignmentlist)")
	gblocked_filelist = []
	for afile in alignmentlist:
		renamesingleseqs(afile, headers)
		mylogger.info("\tGblocking {}".format(afile))
		gblocked_filelist.append(call_Gblocks(afile, len(headers)))
		mylogger.info("\tGblocked file: {}".format(gblocked_filelist[-1]))
	return gblocked_filelist

def renamesingleseqs(afile, headers):
	mylogger.info("renaming seqeunces in {}".format(afile))
	suffix_list = ["_cds_aa.fasta", ".fasta", ".fas", ".faa", ".fa"]#list of most probable input-sequence suffixes for removal from sequence identifier
	alignment = read_alignments([afile])[0]
	for index in range(len(headers)):
		newid = headers[index]
		for suffix in suffix_list: #this removes the most probable suffixes from the sequence names in the final alignments
				if headers[index].endswith(suffix):
					newid = headers[index].rstrip(suffix)
					break
		alignment[index].id = newid
		alignment[index].description = ""
	write_final_alignment(afile, alignment)

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

#CHANGED THE FUNXTION BELOW TO COMBINE GBOLCKS WITH OWN METHOD:
#	-gapped positions have already been previously removed from FLANKING REGIONS of each individual markergene alignment (removing positions not covered by e.g. partial genes)
#	-THEN all all alignments have been concatenated
#	-THEN gblocks is called on this concatement here to remove bad alignment positions (BUT CHANGED flanking cutoff to 50% and KEEP ALL GAPPED POSITIONS!)
def call_Gblocks(file_name, ORG_number): #this calls Gblocks with standard settings. I tried not to overload the argument list for this python script
	mylogger.debug("call_Gblocks(%s)" % file_name)
	gb_cutoff_quotient = 0.5 #fraction (0-9) of residues having to be identical/cimilar to count as "conserved". !!CONSIDER CHANGING THIS TO 0.25 IF POSSIBLE!!
	gb_cutoff_value = int(ORG_number * gb_cutoff_quotient) + 1
	tempfile_name, temp_name_dict = rename_for_gblocks(file_name)
	gblocks_args = ['-t=p', '-e=-gb', '-d=n', '-b1=%s' %gb_cutoff_value, '-b2=%s' %gb_cutoff_value, '-b3=8', '-b4=10', '-b5=a']
	
	gblocks_command = [args.gblocks_binary, tempfile_name] + gblocks_args
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
	#print "\n----------fuuuuu"
	#print alignment_list
	#print "\n---------fuuuuu"
	concatenated_alignment = alignment_list[0][:, 0:0] #simply create an empty an empty alignment with pre-defined placeholders for every comparison-genome (numpy annotation, as alignment-objects are based on numpy)
	
	for alignment in alignment_list:
		counter += 1
		mylogger.debug("concatenating alignment %d of %d with %d comparison-markers" %(counter, maxnumber, len(alignment)))
		concatenated_alignment += alignment
	del alignment_list #not needed anymore --> free up memory
	
	if len(concatenated_alignment) == len(headers): #len(concatenated_alignment) counts the NUMBER of comparison sequences, not their length
		mylogger.debug("GOOD: number_concat_seqs: %s == number headers: %s" %(len(concatenated_alignment),len(headers)))
		for index in range(len(headers)): #renaming the concatenated sequences (otherwise they will just be named after the first ortholog)
			newid = headers[index]
			for suffix in suffix_list: #this removes the most probable suffixes from the sequence names in the final alignments
					if headers[index].endswith(suffix):
						newid = headers[index].rstrip(suffix)
						break
			concatenated_alignment[index].id = newid
			concatenated_alignment[index].description = "concatenated_unique_coregenome"
	else:
		mylogger.debug("number_concat_seqs: %s != number headers: %s" %(len(concatenated_alignment),len(headers)))
		mylogger.debug("first header: %s" % headers[0])
		mylogger.debug("last header: %s" % headers[-1])
		mylogger.debug("first seq: %s " % concatenated_alignment[0])
		mylogger.debug("last seq: %s \n" % concatenated_alignment[-1])
		sys.stdout.write(str(concatenated_alignment))
		sys.stdout.write("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
		raise RuntimeError("Records in concatenated_alignment and headers are not synchronized!")
		
	if len(concatenated_alignment) != 0:
		 aln_length = len(concatenated_alignment[0])
		
	return concatenated_alignment

def write_final_alignment(outputfilename, concatenated_alignment):
	mylogger.debug("write_final_alignment(%s, concatenated_alignment)" % outputfilename)
	outputfile = open(outputfilename, 'w')
	AlignIO.write(concatenated_alignment, outputfilename, "fasta")
	outputfile.close()
#
#IMPORTANT!!!
#DOUBLE CHECK IF TREE-CAL_METHODS DO NOT STILL USE UNFILTERED ALIGNMENT(VEFORE GBLOCKS) !
#IMPORTANT!!!
#
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
	mylogger.debug("deleting temporary files")
	for delfile in os.listdir("."):
		if delfile.startswith("RAxML_") and "." + outname + ".RUN." in delfile:
			#print "deleting temp-file: " + delfile
			os.remove(delfile)
	return outputfiles

def fasttree(alignmentfile):
	mylogger.debug("fasttree(%s)" % alignmentfile)
	outfile = alignmentfile + "_fastree.newick"
	fasttree_command = ["FastTree", "-lg", "-gamma", alignmentfile]
	mylogger.debug("\t{}>{}".format(" ".join(fasttree_command), outfile))
	with open(outfile, "w") as f:
		exitstat = call(fasttree_command, stdout=f)
	mylogger.info(" created {}".format(outfile))
	if exitstat == 1:
		raise Exception("an error occured during fasttree. Maybe you have duplicate identifiers?")
	return outfile
	
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
		mylogger.info("-using aligner '%s'" % alignmeth)
		#mylogger.info("-alignment filter: " + args.afilter)
		headers, MLSA_list, OG_number = read_PO_file(args.PO_file)
		##move this to checkargs()
		
		#read input
		record_dict = read_fasta_seqs(headers, MLSA_list)
		seq_filelist = write_seq_files(record_dict, MLSA_list)
		unaligned_filelist = write_temp_files(record_dict, MLSA_list, "unaligned_temp_fasta_OG")
		
		#alignments
		aligned_filelist = make_alignments(unaligned_filelist)
		mylogger.debug("MAIN: length of aligned_filelist = %s" % len(aligned_filelist))
		if alignmeth == "muscle":
			correct_for_muscle_bug(aligned_filelist, unaligned_filelist) #Necessary because bug in muscle '-stable' option (option disabled for this reason as of muscle version 3.8)
		alignment_list = read_alignments(aligned_filelist)
		mylogger.debug("length alignment_list after reading corrected files back in : %s" %len(alignment_list))
		
		###Now removing all flanking gapped positions of all single alignments prior to gblocks and concatenating
		mylogger.info("\n%s\nRemoving only the flanking gapped positions from all single alignments prior to concatenation" % hline)
		alignment_list = remove_gaps_from_alignment_borders(alignment_list)
		gblocked_alignment_list = rungblocks_on_single_alignments(alignment_list, OG_number, headers)
		#create fastrre phylogenies for each single markergene in order to filter out potentially falsly identified homologs (e.g. horizontal gentransfer) that may distort the phylogeny 
		if args.filter_outliers:
			mylogger.info("\nNOW FILTERING OUT POTENTIAL OUTLIERS")
			singlemarkertreelist = []
			for a in gblocked_alignment_list:
				singlemarkertreelist.append(fasttree(a))
			pf_treelist = [ pf_read_tree(infilename) for infilename in singlemarkertreelist ]
			#if args.rootlabels != None: #TODO: add this option
			#	mylogger.info("\nrerooting trees\n") #TODO: add this option
			#	root_trees(singlemarkertreelist, pf_treelist, args.rootlabels) #TODO: add this option
			pf_depthlist = [ pf_get_depth_sum(tree) for tree in pf_treelist ]
			mylogger.info("total number of single marker trees :\t{}\n".format(len(pf_treelist)))
			mylogger.info(" max depth:\t{}\n".format(max(pf_depthlist)))
			mylogger.info(" min depth:\t{}\n".format(min(pf_depthlist)))
			outlier_indices, mean, stdev = pf_find_outliers_meanbased(pf_depthlist) #only using mean here
			mylogger.info(" mean depth:\t{}\n".format(mean))
			mylogger.info(" standard deviation:\t{}\n".format(stdev))
			mylogger.info("upper cutoff:\t{}\n".format(mean + (2 * stdev)))
			mylogger.info("The following marker-alignments have been identified as outliers and will NOT be used for MLSA:\n{}".format([ " - " + gblocked_alignment_list[x] for x in outlier_indices]))
			gblocked_alignment_list = [ gblocked_alignment_list[x] for x in range(len(gblocked_alignment_list)) if x not in outlier_indices ]
			mylogger.info("--> Removed {} potential outliers from the alignment list prior to concatenating!".format(len(outlier_indices)))
		
		mylogger.info("\n%s\nconcatenating alignments" % hline)
		concatenated_alignment = concatenate_alignments(read_alignments(alignment_list), headers)
		mylogger.info("\n%s\nwriting concatenated alignment to fasta-file: %s" %(hline, outputfilename))
		write_final_alignment(outputfilename, concatenated_alignment)

		
		#give summary
		mylogger.info("=" * 50)
#		mylogger.info("Finished alignments!\nthe individual (unaligned) ortholog selections were stored as the following multifastas: ")
#		mylogger.info(", ".join(seq_filelist))
		mylogger.info("Finished alignments!\nthe individual (unaligned) ortholog selections were stored under the following basename:\n\t%s " % os.path.join(args.out_path, "OGselection_MLSA_single_unaligned_multifasta_*"))
		mylogger.info("\n%d genes concatenated degapped and gblocked to a final sequence length of %d residues per organism" %(OG_number, aln_length))
		mylogger.info("\n-->final aligned and processed sequence saved as: %s" % outputfilename)

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
			mylogger.debug("reading final alignment file: %s" %outputfilename)
			final_align = AlignIO.read(outputfilename, "fasta") 
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
				mylogger.debug("deleting {}".format(delfile))
				os.remove(delfile)
				if "clustalw" in alignmeth:
					os.remove(delfile.rstrip(".fasta") + ".dnd") #remove pesky "guide-tree" files produced by clustal aligners as well
			for delfile in aligned_filelist:
				mylogger.debug("deleting {}".format(delfile))
				os.remove(delfile)
			#if singlemarkertreelist:
			#	for delfile in singlemarkertreelist:
			#		mylogger.debug("deleting {}".format(delfile))
			#		os.remove(delfile)
			
		mylogger.info("cleaning up RaxML tempfiles")
		for delfile in os.listdir("."):
			if "delme_tempfile" in delfile:
				mylogger.debug("deleting {}".format(delfile))
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

def pf_read_tree(infilename): #phylofilter
	infile = open(infilename, "r")
	intree = list(Phylo.NewickIO.parse(infile))
	#print infilename
	#print intree
	if len(intree) > 1:
		sys.stderr.out.write("\nWARNING. {} contains {} trees. Will only consider the first. If these are not permutations of the same tree, please split them into individual files".format(infilenme, len(intree)))
	elif len(intree) == 0:
		raise IOError("ERROR: No tree in {}".format(infilename))
	return intree[0]

def pf_get_depth_sum(tree): #phylofilter
	return tree.total_branch_length()

def pf_get_mean_and_stdev(depthlist): #phylofilter
	#import numpy
	#depthmax = max(depthlist)
	#depthmin = min(depthlist)
	depthmean = numpy.mean(depthlist)
	depthstdev = numpy.std(depthlist, ddof = 1)
	return depthmean, depthstdev

def pf_get_median_and_percentiles(depthlist): #phylofilter
	#import numpy
	median = numpy.median(depthlist)
	upper95perc = numpy.percentile(depthlist, 95)
	#upper75perc = numpy.percentile(depthlist, 75)
	lower25perc = numpy.percentile(depthlist, 25)
	#dmax = max(depthlist
	#dmin = min(depthlist)
	return median, upper95perc
	
def pf_find_outliers_meanbased(depthlist): #phylofilter
	#depthlist = [ get_depth_sum(tree) for tree in treelist ]
	mean, stdev = pf_get_mean_and_stdev(depthlist)
	#lowerbound = mean - (stdev) #not 
	upperbound = mean + (2 * stdev)
	outlier_indices = []
	for i in range(len(depthlist)):
		if depthlist[i] > upperbound:
			outlier_indices.append(i)
	return outlier_indices, mean, stdev

def pf_root_trees(infilenames, treelist, rootlabels): #phylofilter
	for i in range(len(infilenames)):
		outfilename = infilenames[i] + "_rerooted.newick"
		outfile = open(outfilename, "w")
		treelist[i].root_with_outgroup(rootlabels)
		Phylo.NewickIO.write([treelist[i]], outfile)
		outfile.close()

if args.showversion:
	print_version()
else:
	mylogger = setlogger(logfilename)
	main()
