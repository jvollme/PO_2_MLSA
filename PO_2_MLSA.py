#!/usr/bin/python
#created 21.06.2014 by John Vollmers
import os, sys, argparse, time, multiprocessing, random
from Bio.Align.Applications import ClustalwCommandline, MuscleCommandline, ClustalOmegaCommandline
from Bio.Phylo.Applications import RaxmlCommandline
from Bio import AlignIO, SeqIO
from Bio.Alphabet import IUPAC
from subprocess import call

myparser=argparse.ArgumentParser(description="\n==PO_2_MLSA.py v1.2 by John Vollmers==\nCreates concatenated alignments of UNIQUE Genes with orthologues in comparison organisms for the construction of MLSA-based phylogenetic trees.\nOptionally the resulting concatenated alignments may contain all gapped alignmentpositions or may be stripped either of ALL gapped positions or of all gapped positions in the flanking regions of each composite ortholog\nThis script is supposed to be part of a pipeline consisting of:\n\tA.)Conversion of Genbank/Embl-Files to ANNOTATED(!) Fastas using CDS_extractor.pl by Andreas Leimbach\n\tB.)Calculation of orthologs and paralogs using proteinortho5 (WITH the '-single' and '-self' arguments!)\n\tC.)The creation of concatenated MLSA-sequences based on:\n\t\t-the fasta sequences of step A\n\t\t-the proteinortho5-results from step B\n\nThe output-file will be in fasta format (but including gapped positions, so remember to use 'fasta_wgap' when loading into Arb!). However it's absolutely no problem to include other common output-alignmentformats on request!", formatter_class=argparse.RawTextHelpFormatter)
myparser.add_argument("-po", "--proteinortho", action="store", dest="po_resultfile", help="(String) file with proteinortho5 results", required=True)
myparser.add_argument("-f", "--fastas", action="store", dest="fasta_path", help="(String) path to fasta files (produced by CDS_extractor) \nDefault=current working directory", default="")
myparser.add_argument("-tp", "--temp_path", action="store", dest="temp_path", default="", help="Path for temporary files (will be created if does not exist)\nDefault= current working directory")
myparser.add_argument("-kt", "--keep_temp", action="store_true", dest="keep_temp", default=False, help="Keep all temporary files and intermediate results\n(Multifasta_files containing the singe genes involved in the MLSA-alignments are stored in any case)\nDefault: delete all temporary files")
myparser.add_argument("-am", "--alignmethod", action="store", dest="alignmeth", choices=["muscle","clustalw","clustalw2","clustalo"], default="muscle", help="The Alignmentmethod to use.\nDefault='muscle'")
myparser.add_argument("-ap", "--aligner_path", action="store", dest="aligner_path", default="", help="(OPTIONAL: set path to the aligner of choice IF not included in global PATH variable")
myparser.add_argument("-dg", "--degap", action="store", dest="degap",choices=["none", "all", "flanking"], default="all", help="(Only meant for use, if Gblocks is not installed\nSpecify if and which gaps to remove:\n\t'none' keep all gapped positions in the final MLSA-alignment\n\t'all' remove ALL gapped postitions in the alignments\n\t'flanking' remove flanking gapped positions from all individual alignments\nDefault='all'")
myparser.add_argument("-s", "--silent", action="store_true", dest="no_verbose", help="non-verbose mode")
#todo: implement an own "mutlrithreading" for alignment-step (start nthreads alignment-processes simultaniously. wait for all processes to finish before carrying on)
myparser.add_argument("-t", "--threads", action="store", dest="nthreads", type=int, default=1, help="Maximum number of threads to use for alignment steps\nDefault=1")
myparser.add_argument("-gb", "--gblocks", action="store", dest="gblocks", choices=["n","no","f","false","y","yes", "t", "true"], default="true", help="calls gblocks (if installed) to remove gapped positions and poorly aligned regions\n(Overrides '-dg'|'--degap'\nchoices:\n\t[n|no|f|false]: will NOT use Gblocks\n\t[y|yes|t|true]: WILL use Gblocks\nDefault=true (WILL use Gblocks)")
myparser.add_argument("-gbp", "--gblocks_path", action="store", dest="gblocks_path", default="", help="(OPTIONAL: set path to Gblocks IF not included in global PATH variable)")
myparser.add_argument("-op", "--out_path", action="store", dest="out_path", default="", help="Path to output (will be created if it does not exist)\nDefault=current working directory")
#V still to implement:
myparser.add_argument("-mt", "--make_tree", action="store", dest="tree_method", choices=["raxml", "raxml_bs", "raxml_rapidbs", "none"], default="none", help="Generate ML phylogenetic trees using RAxML with the substitution model \"PROTGAMMAAUTO\"\n\tchoices:\t\"raxml\": single tree without bootstraps (using new rapid hill climbing)\n\t\traxml_bs: thorough bootstrap analyses and search for best ML tree\n\t\traxml_rapidbs: rapid bootstrap analyses and search for best ML tree in one run\n\t\tnone\nDefault=none")
myparser.add_argument("-tbp","--tree_builder_path", action="store", dest="treebuilder_path", default="", help="Path to treebuilder (currently only raxml supported) if not listed in $PATH")
myparser.add_argument("-sd", "--seed", action="store", dest="seed_nr", type=int, default=0, help="Integer to provide as seed for RAxML\n0=seed generated randomly\nDefault=random seed")
myparser.add_argument("-bs", "--bootstraps", action="store", dest="nr_bootstraps", type=int, default=1000, help="Number of bootstraps(if any)\ndefault=1000")
#myparser.add_argument("-ctba", "--custom_tree_builder_args", action="store", dest="custom_tree_builder_args", default=None, help="custom arguments for raxml. CAUTION: will overide Only use if you know ExACTLY what you are doing!")
#^ still to implement
myparser.add_argument("-v", "--version", action="store_true", dest="showversion", default=False, help="show version information and then quit (don't run complete script)")
args=myparser.parse_args()

#TOdo: add option "return_selection" to store selection of MLSA genes as unagligned multifastas or only lists of fasta-headers

version="v1.04"
available_cores=multiprocessing.cpu_count() #counts how many cores are available, to check if the user-argument for threads can be fulfilled
aln_length, proc_aln_length, OG_number=0,0,0
wstrings, estrings, lstrings=[], [], [] #warning, error and log messages respectively
out_path, temp_path, fasta_path, PO_file=args.out_path, args.temp_path, args.fasta_path, args.po_resultfile
alignmeth, aligner_path, treebuilder_path=args.alignmeth, args.aligner_path, args.treebuilder_path
tree_method, seed=args.tree_method, args.seed_nr
degap=args.degap
nthreads=args.nthreads
tree_method, bootstraps, seed = args.tree_method, args.nr_bootstraps, args.seed_nr
keep_temp=args.keep_temp
outputfilename=os.path.join(out_path,"concatenated_orthologs_"+alignmeth+"_"+os.path.basename(PO_file)+".fasta")
logfilename=os.path.join(out_path,"PO_2_MLSA_"+time.strftime("%Y%m%d%H%M%S")+".log")
raxml_prog="raxmlHPC"
verbose=True
docontinue=True
gblocks=True
if args.gblocks in ["n", "no", "f", "false"]:
	gblocks=False
gblocks_path=args.gblocks_path

def checkargs(args):
	#print "CHECKING ARGS"
	global verbose, nthreads, tree_method, raxml_prog
	if args.no_verbose:
		verbose=False
	if available_cores<nthreads:
		datsAwarning("WARNING: less than "+str(nthreads)+" cores/threads available! Setting nthreads to "+str(available_cores))
		nthreads=available_cores
	if fasta_path!="" and (not os.path.exists(fasta_path) or not os.path.isdir(fasta_path)):
		datsANerror("ERROR: fasta_path: '"+fasta_path+" does not exist or is no directory!")
	if aligner_path!="":
		if os.path.exists(aligner_path) and os.path.isdir(aligner_path):
			if not os.path.exists(os.path.join(aligner_path, alignmeth)):
				datsANerror("ERROR: can't find '"+alignmeth+"' at "+aligner_path)
		else:
			datsANerror("ERROR: aligner_path: '"+aligner_path+"' does not exist or is not a directory!")
	elif aligner_path=="":
		test_aligner=which(alignmeth)
		if test_aligner==None:
			datsANerror("Can't find "+alignmeth+" in any directory in the PATH variable. Please provide a path")
	if not os.path.exists(PO_file) or not os.path.isfile(PO_file):
		datsANerror("ERROR: cannot find proteinortho-resultfile: "+PO_file)
	if temp_path!="" and docontinue: 
		if not os.path.exists(temp_path) or not os.path.isdir(temp_path):
			if verbose: 
				print "Creating directory for temporary and intermediate result files: "+temp_path
			os.mkdir(temp_path)
		datsAlogmessage("-Will store temporary and intermediate result files in "+os.path.abspath(temp_path))
	if out_path!="" and docontinue: 
		if not os.path.exists(out_path) or not os.path.isdir(out_path):
			if verbose: 
				print "Creating directory for final result files: "+out_path
			os.mkdir(out_path)
		datsAlogmessage("-Will store final result files in "+os.path.abspath(out_path))
	if alignmeth=="clustalo" and docontinue: #check clustalo version
		clustalomega_cline=ClustalOmegaCommandline("clustalo", version=True)
		clustalo_version=clustalomega_cline()[0].rstrip().split(".")
		if int(clustalo_version[0])<1 or (int(clustalo_version[0])==1 and int(clustalo_version[1])<2) : #only accept versions 1.2 and newer
			datsANerror("ERROR: installed clustalo version is v"+".".join(clustalo_version)+"! Version 1.2 or higher is required!")
	if gblocks and docontinue:
		if gblocks_path=="":
			test_gblocks=which("Gblocks")
			if test_gblocks==None:
				datsANerror("ERROR: can't locate Gblocks in any path in PATH variable. please provide a Path to Gblocks using the '-gbp' agrument!")
			elif verbose:
				print "Located Gblocks executable: "+test_gblocks
		else:
			if os.path.exists(gblocks_path) and os.path.isdir(gblocks_path):
				if os.path.exists(os.path.join(gblocks_path, "Gblocks")) and os.path.isfile(os.path.join(gblocks_path, "Gblocks")):
					if verbose:
						print "Located Gblocks executable: "+os.path.join(gblocks_path, "Gblocks")
			elif gblocks_path.endswith("Gblocks") and os.path.isfile(gblocks_path):
				if verbose:
					print "Located Gblocks executable: "+gblocks_path
			else:
				datsANerror("ERROR: Gblocks executable could not be found in the specified path: "+gblocks_path)
	if tree_method!="none":
		#print "CHECKING raxml-binaries"
		if treebuilder_path=="":
			if which("raxmlHPC")==None and which("raxmlHPC-PTHREADS")==None and which("raxmlHPC-PTHREADS-SSE3")==None and which("raxmlHPC-SSE3")==None:
				datsANerror("ERROR: No raxmlHPC binaries found in any directory within $PATH! please provide a PATH to raxml binaries!")
				tree_method="none"
			else:
				if which("raxmlHPC-PTHREADS-SSE3")!=None:
					raxml_prog=which("raxmlHPC-PTHREADS-SSE3")
					datsAlogmessage("found "+raxml_prog+ "! will use this tool") 
				elif which("raxmlHPC-PTHREADS")!=None:
					raxml_prog=which("raxmlHPC-PTHREADS")
					datsAlogmessage("found "+raxml_prog+ "! will use this tool")
				else:
					datsAwarning("Warning: multithreading is only supported with 'PTHREADS'-versions of raxml. Not sure if your raxml binaries support this.\t\nif raxml calculations fail, recombile raxml with 'PTHREADS'-option")
					if which("raxmlHPC-SSE3")!=None:
						raxml_prog=which("raxmlHPC-SSE3")
						datsAlogmessage("found "+raxml_prog+ "! will use this tool")
					else:
						raxml_prog=which("raxmlHPC")
						datsAlogmessage("found "+raxml_prog+ "! will use this tool")
				try:
					checkraxml_cline=RaxmlCommandline(version=True)
					versiontext=checkraxml_cline()[0]
					startv=versiontext.find("version ")+len("version ")
					endv=versiontext[startv:].find(" ")
					version=versiontext[startv:startv+endv].split(".")
					if int(version[0])<8 or (int(version[0])==8 and int(version[1])==0 and int(version[2])<20):
						datsAwarning("Warning: This script was devised for and tested with RAxML v8.0.20. Your version is v"+".".join(version)+" !\n\tThis may very well still work, but if it doesn't it's YOUR fault!")
				except:
					datsAwarning("Warning: This script was devised for and tested with RAxML v8.0.20.\n\tNot sure which version of RAxML you're using, but it sure as hell isn't v7 or v8!\n\tThis may very well still work, but if it doesn't it's YOUR fault!")
		elif os.path.exists(treebuilder_path) and os.path.isfile(treebuilder_path):
			if nthreads>1 and not "PTHREADS" in treebuilder_path:
				datsAwarning("Warning: multithreading is only supported with 'PTHREADS'-versions of raxml. Not sure if your choosen binaries support this.\t\nif raxml calculations fail, recombile raxml with 'PTHREADS'-option") 
			try:
				checkraxml_cline=RaxmlCommandline(version=True)
				versiontext=checkraxml_cline()[0]
				startv=versiontext.find("version ")+len("version ")
				endv=versiontext[startv:].find(" ")
				version=versiontext[startv:startv+endv].split(".")
				if int(version[0])<8 or (int(version[0])==8 and int(version[1])==0 and int(version[2])<20):
						datsAwarning("Warning: This script was devised for and tested with RAxML v8.0.20. Your version is v"+".".join(version)+" !\n\tThis may very well still work, but if it doesn't it's YOUR fault!")
				raxml_prog=treebuilder_path
			except:
				datsAwarning("Warning: Correct raxML-version not found under "+treebuilder_path+"!\nWill NOT calculate ML trees!")
				tree_method="none" 
				
def datsAwarning(wstring):
	global wstrings
	wstrings.append(wstring)
	if verbose:
		print "\n"+wstring+"\n"
	
def datsANerror(estring):
	global estrings
	global docontinue
	estrings.append(estring)
	docontinue=False
	if verbose:
		print >>sys.stderr, "\n"+estring+"\n"

def datsAlogmessage(lstring):
	global lstrings
	lstrings.append(lstring)
	if verbose:
		print lstring

def which(file):
    for path in os.environ["PATH"].split(":"):
        if os.path.exists(path + "/" + file):
                return path + "/" + file
    return None

def read_PO_file(filename):
	open_PO_file=open(filename, 'r')
	firstline=True
	org_index=3 #default index of beginning of organism-result-columns in proteinortho4+
	headers=[]
	MLSA_gene_counter=0
	MLSA_list=[]
	for line in open_PO_file:
		if firstline:
			if line.startswith("#species"):
				datsAwarning("WARNING! It seems you are using results from Proteinortho4. It's recommended to use Proteinortho5 or later!\n\t(However, the derivation of MLSA genes should still work)")
			elif line.startswith("# Species"):
				if verbose:
					print "\nProteinortho-results seem to be based on Proteinortho5 or later. Good."
			else:
				datsAwarning("WARNING! Cannot clearly recognize Format of Proteinortho-results! This may produce erroneous results!\n\t For best results use Proteinortho5!") 
			headers=line.split("\t")[org_index:]
			#remove line_end_symbols:
			for h in range(len(headers)):
				headers[h]=headers[h].rstrip()
			MLSA_list=[[header,[]] for header in headers]
			firstline=False
		elif not line.startswith("#"):
			if not "," in line:
				if not "*" in line:
					zeilentokens=line.split("\t")[org_index:]
					#remove line end symbols:
					for zt in range(len(zeilentokens)):
						zeilentokens[zt]=zeilentokens[zt].rstrip()
					if len(zeilentokens)!=len(headers):
						datsANerror("ERROR: Different numbers of columns in header_line and Locus_tag_lines in "+filename)
						break
					for h in range(len(headers)):
						MLSA_list[h][1].append(zeilentokens[h].rstrip())
	open_PO_file.close()
	#check if everything is ok before continuing
	for org in MLSA_list:
		if len(org[1])!=len(MLSA_list[0][1]):
			datsANerror("ERROR: something went wrong while reading proteinortho-results. Not the same number of comparison genes for all Organisms!")
	datsAlogmessage("\nrecognized "+str(len(MLSA_list[0][1]))+" unambigious (=unique) 'Orthologeous Groups' (OGs) shared by the comparison organisms")
	OG_number=len(MLSA_list[0][1])
	return headers, MLSA_list, OG_number

def read_fasta_seqs(headers, MLSA_list):
	record_dict=dict.fromkeys(headers, {})
	for org in record_dict:
		record_dict[org]={}
	for h in range(len(headers)):
		if MLSA_list[h][0]!=headers[h]:
			datsANerror("ERROR: MLSA_list and headers are not synchronized! WTF!?!?!?!?")
			break
		if not os.path.exists(os.path.join(fasta_path, headers[h])) or not os.path.isfile(os.path.join(fasta_path, headers[h])):
			datsANerror("ERROR: Can't find fasta file: "+(os.path.join(fasta_path, headers[h]))+"!\n\tPlease provide a path to a directory, that contains all fasta-files listed in "+PO_file)
			return
		open_fastafile=open(os.path.join(fasta_path, headers[h]), 'r')
		temprecord_dict=dict.fromkeys(MLSA_list[h][1], None)
		for record in SeqIO.parse(open_fastafile, "fasta", alphabet=IUPAC.protein):
			if record.id in temprecord_dict:
				record_dict[headers[h]][record.id]=record
		#add another check_function here!!
		open_fastafile.close()
	return record_dict
	
def write_seq_files(record_dict, MLSA_list):
	selection_filelist=[]
	for orgindex in range(len(MLSA_list)):
		tempfilename=os.path.join(out_path, "OGselection_MLSA_single_unaligned_multifasta_"+MLSA_list[orgindex][0])
		selection_filelist.append(tempfilename)
		open_tempfile=open(tempfilename, 'w')
		for locus in MLSA_list[orgindex][1]:
			SeqIO.write([record_dict[MLSA_list[orgindex][0]][locus]], open_tempfile, "fasta")
		open_tempfile.close()
	return selection_filelist

def write_temp_files(record_dict, MLSA_list, prefix):
	unaligned_filelist=[]
	OG_counter=0 #numbering of the used 'Orthologeous Groups'(OGs), to differentiate between the different temp_files 
	for locusindex in range(len(MLSA_list[0][1])):
		OG_counter+=1
		tempfilename=os.path.join(temp_path, prefix+str(OG_counter).zfill(5)+".fasta")
		unaligned_filelist.append(tempfilename)
		open_tempfile=open(tempfilename, 'w')
		for orgindex in range(len(MLSA_list)):
			SeqIO.write([record_dict[MLSA_list[orgindex][0]][MLSA_list[orgindex][1][locusindex]]], open_tempfile, "fasta")
		open_tempfile.close()
	return unaligned_filelist

def write_temp_alignments(alignmentlist, prefix):
	proc_aligned_filelist=[]
	OG_counter=0 #numbering of the used 'Orthologeous Groups'(OGs), to differentiate between the different temp_files 
	for index in range(len(alignmentlist)):
		OG_counter+=1
		alignfilename=os.path.join(temp_path, prefix+str(OG_counter).zfill(5)+".fasta")
		proc_aligned_filelist.append(alignfilename)
		open_tempfile=open(alignfilename, 'w')
		AlignIO.write([alignmentlist[index]], open_tempfile, "fasta")
		open_tempfile.close()
	return proc_aligned_filelist

def make_alignments(unaligned_filelist):
	if verbose:
		print "\n----------------------------------"
		print "Aligning Orthologeous Groups (OGs)"
	aligned_filelist=[]
	counter=0
	for uf in unaligned_filelist:
		counter+=1
		if verbose:
			sys.stdout.write("\raligning OG "+str(counter)+" from "+str(len(unaligned_filelist))+" using "+alignmeth)
			sys.stdout.flush()
		aligned_filelist.append(eval(alignmeth)(uf))
	return aligned_filelist
			 
def clustalw(inputfile):#Todo: change number of threads
		outputfile=inputfile.replace("unaligned_temp_fasta_", "SINGLEalignment_CLUSTALW_temp_fasta_",1)
		clustalw_cline= ClustalwCommandline(os.path.join(aligner_path, alignmeth), INFILE=inputfile, outfile=outputfile, type="PROTEIN", align=True, quiet=True, OUTORDER="INPUT")
		try:
			clustalw_cline()
			return outputfile
		except Exception:
			datsANerror("ERROR: clustalw version is older than v2 (probably v1.83). You should use version 2 or newer (called clustalw2 on many systems)")

def clustalw2(inputfile):#Todo: change number of threads
		outputfile=inputfile.replace("unaligned_temp_fasta_", "SINGLEalignment_CLUSTALW2_temp_fasta_",1)
		clustalw_cline= ClustalwCommandline(os.path.join(aligner_path, alignmeth), INFILE=inputfile, outfile=outputfile, type="PROTEIN", align=True, quiet=True, OUTORDER="INPUT")
		clustalw_cline()
		return outputfile

def clustalo(inputfile):#Todo: find out a way to check clustalo version
		outputfile=inputfile.replace("unaligned_temp_fasta_", "SINGLEalignment_CLUSTALO_temp_fasta_",1)
		clustalomega_cline= ClustalOmegaCommandline(os.path.join(aligner_path, alignmeth), infile=inputfile, outfile=outputfile, seqtype="Protein", threads=nthreads, verbose=False, force=True, outputorder="input-order")
		clustalomega_cline()
		return outputfile

def muscle(inputfile):#changing the number of threads aparently not possible
		outputfile=inputfile.replace("unaligned_temp_fasta_", "SINGLEalignment_MUSCLE_temp_fasta_",1)
		muscle_cline= MuscleCommandline(os.path.join(aligner_path, alignmeth), input=inputfile, out=outputfile, quiet=True) #add 'stable=True' to the end of this list, if the stable-bug in muscle is fixed (remove the correct_for_muscle_bug() method in that case)
		muscle_cline()
		return outputfile

def read_alignments(input_filelist):
	alignmentlist=[]
	if alignmeth=="muscle" or alignmeth=="clustalo" or len(input_filelist)==1: #Last condition assumes that ONLY final result files would be passed as a filelist of only one file
		aformat="fasta"# just a workaround for this method
	else:
		aformat="clustal"
	for alignmentfile in input_filelist:
		currentalignfile=open(alignmentfile, 'r')
		alignmentlist.append(AlignIO.read(currentalignfile, aformat))
		currentalignfile.close()
	if len(alignmentlist)!=len(input_filelist):
		datsANerror("ERROR: could not read in all temporary alignmentfiles!")
	return alignmentlist

def remove_gaps_from_alignment_borders(alignmentlist): #optional. Better to use Gblocks if available
	global outputfilename
	global proc_aln_length
	outputfilename=outputfilename.replace("concatenated_orthologs_", "concatenated_orthologs_FlankingBordersDegapped_")
	processed_alignmentlist=[]
	for currentalignment in alignmentlist:
		start, stop = 0, len(currentalignment[0])
		for c in range(stop): #identify start of alignable region
			if not "-" in currentalignment[:,c]:
				start=c
				break
		for c in range(1, stop): #identify end of alignable region
			nc=c*(-1)
			if not "-" in currentalignment[:,nc]:
				stop+=nc #=maxlength minus endposition
				break
		processed_alignmentlist.append(currentalignment[:,start:stop+1])
	#print "keep_temp: "+str(keep_temp)
	if keep_temp==True:
		#print "keep_temp is 'True' --> will produce prossessed single alignment files"
		write_temp_alignments(processed_alignmentlist, "Processed_removedFLANKINGgaps_SINGLEalignment_MUSCLE_temp_fasta_OG")
	return processed_alignmentlist

def remove_gaps_from_complete_alignments(alignmentlist): #optional. Better to use Gblocks if available
	global outputfilename
	global proc_aln_length
	outputfilename=outputfilename.replace("concatenated_orthologs_", "concatenated_orthologs_AllDegapped_OG")
	processed_alignmentlist=[]
	counter=0
	maxnumber=len(alignmentlist)
	for currentalignment in alignmentlist:
		counter+=1
		sys.stdout.write("\rprocessing alignment "+str(counter)+" of "+str(maxnumber))
		sys.stdout.flush()
		processed_alignment=currentalignment[:,0:0] #initialize alignment without sequence
		for c in range(len(currentalignment[0])):
			if not "-" in currentalignment[:, c]:
				processed_alignment+=currentalignment[:,c:(c+1)]
		processed_alignmentlist.append(processed_alignment)
	#print "keep_temp: "+str(keep_temp)
	if keep_temp==True:
		#print "keep_temp is 'True' --> will produce prossessed single alignment files"
		write_temp_alignments(processed_alignmentlist, "Processed_removedALLgaps_SINGLEalignment_MUSCLE_temp_fasta_OG")
	return processed_alignmentlist
	
def concatenate_alignments(alignment_list, headers):
	global aln_length
	suffix_list=["_cds_aa.fasta", ".fasta", ".fas", ".faa", ".fa"]#list of most probable input-sequence suffixes for removal from sequence identifier
	counter=0
	maxnumber=len(alignment_list)
	concatenated_alignment=alignment_list[0][:,0:0]
	for alignment in alignment_list:
		counter+=1
		sys.stdout.write("\rconcatenating alignment "+str(counter)+" of "+str(maxnumber))
		sys.stdout.flush()
		concatenated_alignment+=alignment
	del alignment_list #not needed anymore --> free up memory
	if len(concatenated_alignment)==len(headers): #len(concatenated_alignment) counts the NUMBER of comparison sequences, not their length
		for index in range(len(headers)): #renaming the concatenated sequences (otherwise they will just be named after the first ortholog)
			newid=headers[index]
			for suffix in suffix_list: #this removes the most probable suffixes from the sequence names in the final alignments
					if headers[index].endswith(suffix):
						newid=headers[index].rstrip(suffix)
						break
			concatenated_alignment[index].id=newid
			concatenated_alignment[index].description="concatenated_unique_coregenome"
	else:
		datsANerror("ERROR: Records in concatenated_alignment and headers are not synchronized!")
	if len(concatenated_alignment)!=0:
		 aln_length=len(concatenated_alignment[0])
	return concatenated_alignment

def write_final_alignment(outputfilename, concatenated_alignment):
	outputfile=open(outputfilename, 'w')
	AlignIO.write(concatenated_alignment, outputfilename, "fasta")
	outputfile.close() 

def write_logfile(logfilename):
	logfile=open(logfilename,'w')
	logfile.write("PO_2_MLSA.py logfile")
	logfile.write("\n"+ time.strftime("%c"))
	logfile.write("\nMessages\n:")
	for l in lstrings:
		logfile.write(l+"\n")
	logfile.write("\n-----------\nWarnings:\n")
	for w in wstrings:
		logfile.write(w+"\n")
	logfile.write("\n------------\nErrors:\n")
	for e in estrings:
		logfile.write(e+"\n")
	logfile.close()

def correct_for_muscle_bug(aligned_filelist, seq_filelist):
	#this def is a necessary workaround for the missing "stable" function of "muscle", which has been deactivated due to the discovery of a bug in this function
	#whithout that function, the sequence order in the alignment might change, and that would be catastrophic for MLSA-analyses!
	#This workaround method could become obsolete, when future muscle versions reimplement a correct '-stable' option (hopefully in version 3.9)
	#In that case, add the option", stable=True" to the MuscleCommandline-call in 'call_muscle()'.
	print "\n----------------------------------"
	print "correcting sequence order in muscle alignments"
	for f in range(len(aligned_filelist)):
		sys.stdout.write("\rcorrecting alignmentfile "+str(f+1)+" of "+str(len(aligned_filelist)))
		alignmentfile=aligned_filelist[f]
		alignmenthandle=open(alignmentfile, 'r')
		alignment=AlignIO.read(alignmenthandle, "fasta")
		alignmenthandle.close()
		seqs=[]
		seq_file=seq_filelist[f]
		seq_handle=open(seq_file, 'r')
		for record in SeqIO.parse(seq_handle, "fasta"):
			seqs.append(record)
		seq_handle.close()
		corrected=None
		for x in range(len(seqs)):
			for y in range(len(alignment)):
				#print str(x)+" "+seqs[x].id+ " == "+alignment[y].id
				if seqs[x].id==alignment[y].id:
					#print "True"
					if x==0:
						corrected=alignment[y:y+1]
					else:
						corrected.append(alignment[y])
				#else:
				#	print "False"
		alignmenthandle=open(alignmentfile, 'w')
		AlignIO.write([corrected], alignmenthandle, "fasta")
		alignmenthandle.close()	
	return corrected

def rename_for_gblocks(align_file):
	#temporarily rename sequences, so that gblocks doesn't freak out
	print "renaming sequences before gblocks"
	alignment = read_alignments([align_file])[0]
	index = 1
	temp_name_dict = {}
	tempfile_name = align_file + "_" + "deltemp_indexed_alignments_" + time.strftime("%Y%m%d%H%M%S")
	for al in alignment:
		print "renaming " + al.id + " to " + str(index)
		temp_name_dict[str(index)] = al.id
		al.id = str(index)
		al.description = ""
		index += 1
	write_final_alignment(tempfile_name, alignment)
	#print "rename_dict:"
	#print temp_name_dict
	return tempfile_name, temp_name_dict
	
def rename_after_gblocks(align_file, temp_name_dict, final_filename):
	print "renaming sequences back after gblocks"
	alignment = read_alignments([align_file])[0]
	for al in alignment:
		#print "renaming " + al.id + " to " + temp_name_dict[al.id]
		al.id = temp_name_dict[al.id]
		al.description = "gblocks-filtered concatenated coregenome"
	write_final_alignment(final_filename, alignment)

def call_Gblocks(file_name): #this calls Gblocks with standard settings. I tried not to overload the argument list for this python script 
	tempfile_name, temp_name_dict = rename_for_gblocks(file_name)
	gblocks_args=['-t=p', '-e=-gb', '-d=n']
	gblocks_command=[os.path.join(gblocks_path,"Gblocks"), tempfile_name]+gblocks_args
	#print "TEST: "+" ".join(gblocks_command)
	call(gblocks_command)
	rename_after_gblocks(tempfile_name+"-gb", temp_name_dict, file_name + "-gb")
	os.remove(tempfile_name)
	os.remove(tempfile_name+"-gb")
	return file_name+"-gb"
	
def randomnumber():
	#returns a random integer to use as seed for rxml and pyml
    random.seed()
    return random.randint(1,2000)

def raxml_rapidbs(alignmentfile, outputfilename, seed, parameters): #parameters should be a dictionary (This dictionary thing was introduced, so that the script can be more easily adapted to accept custom commandline-parameters for raxml by the user)
	nr_threads=4
	if "-T" in parameters:
		nr_threads=parameters["-T"]
	bootstraps=1000
	if "-N" in parameters:
		bootstraps=parameters["-N"]
	elif "-#" in parameters:
		bootstraps=parameters["-#"]
	
	datsAlogmessage("Calculating phylogenies: 'rapid bootstrap analyses and search for best-scoring ML Tree in one run' using raxmlHPC and "+ str(nr_threads)+" threads")
	try:
		outname="rapidBS_"+time.strftime("%Y%m%d%H%M%S")+"_"+"final_tree"	
		raxml_cline=RaxmlCommandline(raxml_prog, sequences=alignmentfile, algorithm="a", model="PROTGAMMAAUTO", name=outname, parsimony_seed=seed, rapid_bootstrap_seed=seed, num_replicates=bootstraps, threads=nr_threads) 
		datsAlogmessage("-->"+str(raxml_cline))
		raxml_cline()
		datsAlogmessage("-->SUCCESS")
		#the resultfiles will be: "RAxML_bipartitions.rapidBS_final_tree" and "RAxML_bipartitionsBranchLabels.rapidBS_final_tree"
		#Labels on nodes or branches, respectively
		outputfiles=["RAxML_bipartitions."+outname, "RAxML_bipartitionsBranchLabels."+outname]
	except Exception as e:
		#note for future versions: create a custom exception "raxml_except" instead of using this construct over and over!
		datsAlogmessage("-->FAILURE")
		datsAwarning("WARNING: rapid bootstrap analyses and search for best-scoring ML tree using raxmlHPC has failed!\n\t"+str(e))	
		return None
	return outputfiles
				
def raxml_bs(alignmentfile, outputfilename, seed, parameters):
	nr_threads=4
	if "-T" in parameters:
		nr_threads=parameters["-T"]
	bootstraps=1000
	if "-N" in parameters:
		bootstraps=parameters["-N"]
	elif "-#" in parameters:
		bootstraps=parameters["-#"]
	datsAlogmessage("Calculating phylogenies: Thorough bootstrap analyses with raxml")
		
	datsAlogmessage("\tDetermining best ML tree of 20 raxmlHPC runs using "+str(nr_threads)+" threads") 
	try:
		raxml_cline=RaxmlCommandline(raxml_prog, model="PROTGAMMAAUTO", name="best_delme_tempfile", parsimony_seed=seed, num_replicates=20, sequences=alignmentfile, threads=nr_threads) 
		datsAlogmessage("\t-->"+str(raxml_cline))
		raxml_cline()
		#the resultfile will be :"RAxML_bestTree.best_delme_tempfile"
		datsAlogmessage("\t-->SUCCESS")
	except Exception as e:
		datsAlogmessage("\t--FAILURE")
		datsAwarning("WARNING: thorough bootstrap analyses using raxmlHPC failed!\n\t"+str(e))
		return None
	
	datsAlogmessage("\tDoing bootstrap analyses with "+str(bootstraps)+" runs using raxmlHPC using "+str(nr_threads)+" threads")
	try:
		raxml_cline=RaxmlCommandline(raxml_prog, model="PROTGAMMAAUTO", sequences=alignmentfile, name="boot_delme_tempfile", parsimony_seed=seed, bootstrap_seed=seed, num_replicates=bootstraps, threads=nr_threads)
		datsAlogmessage("\t-->"+str(raxml_cline))
		raxml_cline()
		#the resultfile will be: "RAxML_bootstrap.boot_delme_tempfile"
		datsAlogmessage("\t-->SUCCESS")
	except Exception as e:
		datsAlogmessage("\t-->FAILURE")
		datsAwarning("WARNING: thorough bootstrap analyses using raxmlHPC failed!\n\t"+str(e))
		return None
		
	datsAlogmessage("\tDrawing bipartitions of bootstrap trees onto best ML tree using raxmlHPC using "+str(nr_threads)+" threads")
	try:
		outname="BS"+str(bootstraps)+"_"+time.strftime("%Y%m%d%H%M%S")+"_"+"final_tree"
		raxml_cline=RaxmlCommandline(raxml_prog, model="PROTGAMMAAUTO", parsimony_seed=seed, algorithm="b", starting_tree="RAxML_bestTree.best_delme_tempfile", bipartition_filename="RAxML_bootstrap.boot_delme_tempfile", name=outname)
		datsAlogmessage("\t-->"+str(raxml_cline))
		raxml_cline()
		#The resultfiles will be: RAxML_bipartitions.final_tree" and "RAxML_bipartitionsBranchLabels.final_tree"
		outputfiles=["RAxML_bipartitions."+outname,"RAxML_bipartitionsBranchLabels."+outname]
		datsAlogmessage("\t-->SUCCESS")	
	except Exception as e:
		datsAlogmessage("\t-->FAILURE")
		datsAwarning("WARNING: thorough bootstrap analyses using raxmlHPC failed!\n\t"+str(e))
		return None
	return outputfiles
			
def raxml(alignmentfile, outputfilename, seed, parameters):
	nr_threads=4
	if "-T" in parameters:
		nr_threads=parameters["-T"]
	try:
		outname="raxml_"+time.strftime("%Y%m%d%H%M%S")+"_"+"final_tree"
		datsAlogmessage("Calculating phylogeny: Determining best ML tree of 20 raxmlHPC runs using "+str(nr_threads)+" threads")
		raxml_cline=RaxmlCommandline(raxml_prog, sequences=alignmentfile, model="PROTGAMMAAUTO", name=outname,  parsimony_seed=seed, num_replicates=20, threads=nr_threads)
		datsAlogmessage("\t-->"+str(raxml_cline))
		raxml_cline()
		#the resultfile will be :"RAxML_bestTree.final_tree"
		outputfiles=["RAxML_bestTree."+outname]
		datsAlogmessage("\tSUCCESS")
		print "deleting temporary files"
		for delfile in os.listdir("."):
			if delfile.startswith("RAxML_") and "."+outname+".RUN." in delfile:
				#print "deleting temp-file: "+delfile
				os.remove(delfile)
	except Exception as e:
		datsAlogmessage("\t-->FAILURE")
		datsAwarning("WARNING: searching for best ML tree using raxmlHPC failed!\n\t"+str(e))
		return None
	return outputfiles
	
def print_version():
	print "\n ==PO_2_MLSA.py "+version+" by John Vollmers==\n"

def main():
	global seed, docontinue
	if verbose:
		print_version()
	checkargs(args)
	if docontinue:
		datsAlogmessage("-keep temporary files: "+str(keep_temp))
		if not keep_temp and verbose:
			print "-->will delete all temporary files"
		datsAlogmessage("-using aligner '"+alignmeth+"'")
		if gblocks in ["no","n","false","f"]:
				datsAlogmessage("-remove gaps: "+degap)
		else:
			datsAlogmessage("-remove gaps: <OVERRIDDEN>\n-->will use Gblocks (with standard settings) on final result files for removal of gapped positions and poorly conserved regions")
		if docontinue:
			headers, MLSA_list, OG_number=read_PO_file(PO_file)
			if docontinue:
				record_dict=read_fasta_seqs(headers, MLSA_list)
				if docontinue:
						seq_filelist=write_seq_files(record_dict, MLSA_list)			
						unaligned_filelist=write_temp_files(record_dict, MLSA_list,"unaligned_temp_fasta_OG")
						aligned_filelist=make_alignments(unaligned_filelist)
						if alignmeth=="muscle":
							correct_for_muscle_bug(aligned_filelist, unaligned_filelist) #Necessary because bug in muscle '-stable' option (option disabled for this reason as of muscle version 3.8)
						alignment_list=read_alignments(aligned_filelist)
						if docontinue:
							if not gblocks and degap=="all":
								if verbose:
									print "\n----------------------------------"
									print "Removing all gapped positions from all single alignments prior to concatenation"
								alignment_list=remove_gaps_from_complete_alignments(alignment_list)
							elif not gblocks and degap=="flanking":
								if verbose:
									print "\n----------------------------------"
									print "Removing only the flanking gapped positions from all single alignments prior to concatenation"
								alignment_list=remove_gaps_from_alignment_borders(alignment_list)
							elif not gblocks and verbose:
								print "\nLeaving the single alignments as they are (Not removing any gapped or unconserved positions)"
							if docontinue:
								if verbose:
									print "\n----------------------------------"
									print "concatenating alignments"
								concatenated_alignment=concatenate_alignments(alignment_list, headers)
								if docontinue:
									if verbose:
										print "\n----------------------------------"
										print "writing concatenated alignment to fasta-file: "+outputfilename
									write_final_alignment(outputfilename, concatenated_alignment)
									if gblocks:
										print "running gblocks on "+outputfilename
										call_Gblocks(outputfilename)
										try:
											final_check=read_alignments([outputfilename+"-gb"])
											proc_aln_length=len(final_check[0][0])
										except Exception as e:
											datsANerror("ERROR: Something wrong with the gblocks output. Please check version and executable rights of gblocks")	
											datsANerror(str(e))
										if docontinue:
											datsAlogmessage("-->Processed File: "+outputfilename+"-gb\n-->Gblocks-Logfile: "+outputfilename+"-gb.htm")

		if keep_temp!=True:
			print "\n----------------------------------"
			print "cleaning up..."
			try:
				for delfile in unaligned_filelist:
					os.remove(delfile)
					if "clustalw" in alignmeth:
						os.remove(delfile.rstrip(".fasta")+".dnd") #remove pesky "guide-tree" files produced by clustal aligners as well	
				for delfile in aligned_filelist:
					os.remove(delfile)
			except:
				print "Nothing to delete, due to aborted run"
		if docontinue:
			if verbose:
				print "=================================="
			datsAlogmessage("FINISHED!\nthe individual (unaligned) ortholog selections were stored as the following multifastas: ")
			datsAlogmessage(", ".join(seq_filelist))
			datsAlogmessage(str(OG_number)+" genes concatenated to a final sequence length of "+str(aln_length)+ " residues per organism")
		if gblocks:
			datsAlogmessage("after 'gblocks' the remaining sequence length is "+str(proc_aln_length)+" residues")
			datsAlogmessage("\n-->final aligned and processed sequence saved as: "+outputfilename+"-gb")
		else:
			datsAlogmessage("\n-->final aligned sequence saved as: "+outputfilename)
			datsAlogmessage("\nRemember to set sequence-type as 'PROTEIN' and filetype as 'fasta_wgap' when loading the MLSA sequence into Arb!\n")
		
		if tree_method!="none":
			datsAlogmessage("\nWAIT, I wasn't finished after all!\nNow generating maximum likelihood trees from final output alignment")
			if seed==0:
				seed=randomnumber()
			treefiles=[]
			outtreename=outputfilename+"."+tree_method+".TREE"
			treefiles.extend(eval(tree_method)(outputfilename, outtreename, seed, {"-N":bootstraps, "-T":nthreads}))
			if docontinue:
				datsAlogmesage("the following trees were generated:")
				datsAlogmessage("\t"+"\n\t".join(treefiles))
			else:
				datsAlogmessage("some kind of error occured during ML calculation")
			datsAlogmessage("cleaning up again...")
			for delfile in os.listdir("."):
				if "delme_tempfile" in delfile:
					os.remove(delfile)
			
	if len(wstrings)>0 or len(estrings)>0:
		if len(wstrings)>0:
			print "\nThere were Warnings!"
		if len(estrings)>0:
			print "\nThere were Errors!"
	print "See "+logfilename+" for details\n"
	write_logfile(logfilename)

if args.showversion:
	print_version()
else:
	main()
