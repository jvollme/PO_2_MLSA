#!/usr/bin/python
#created 21.06.2014 by John Vollmers
import os, sys, argparse
from Bio.Align.Applications import ClustalwCommandline, MuscleCommandline
from Bio import AlignIO, SeqIO
from Bio.Alphabet import IUPAC

myparser=argparse.ArgumentParser(description="\n==PO_2_MLSA.py v.1.0 by John Vollmers==\nCreates concatenated alignments of UNIQUE Genes with orthologues in comparison organisms for the construction of MLSA-based phylogenetic trees.\nOptionally the resulting concatenated alignments may contain all gapped alignmentpositions or may be stripped either of ALL gapped positions or of all gapped positions in the flanking regions of each composite ortholog\nThis script is supposed to be part of a pipeline consisting of:\n\tA.)Conversion of Genbank/Embl-Files to ANNOTATED(!) Fastas using CDS_extractor.pl by Andreas Leimbach\n\tB.)Calculation of orthologs and paralogs using proteinortho5 (WITH the '-single' and '-self' arguments!)\n\tC.)The creation of concatenated MLSA-sequences based on:\n\t\t-the fasta sequences of step A\n\t\t-the proteinortho5-results from step B\n\nThe output-file will be in fasta format (but including gapped positions, so remember to use 'fasta_wgap' when loading into Arb!). However it's absolutely no problem to include other common output-alignmentformats on request!", formatter_class=argparse.RawTextHelpFormatter)
myparser.add_argument("-po", "--proteinortho", action="store", dest="po_resultfile", help="(String) file with proteinortho5 results", required=True)
myparser.add_argument("-f", "--fastas", action="store", dest="fasta_path", help="(String) path to fasta files (produced by CDS_extractor) \nDefault=current working directory", default="")
myparser.add_argument("-tp", "--temp_path", action="store", dest="temp_path", default="", help="Path for temporary files\nDefault= current working directory")
myparser.add_argument("-kt", "--keep_temp", action="store_true", dest="keep_temp", default=False, help="Keep all temporary files and intermediate results\n(Multifasta_files containing the singe genes involved in the MLSA-alignments are stored in any case)\nDefault: delete all temporary files")
myparser.add_argument("-am", "--alignmethod", action="store", dest="alignmeth", choices=["muscle","clustal"], default="clustal", help="The Alignmentmethod to use.\nDefault='clustal'")
myparser.add_argument("-ap", "--aligner_path", action="store", dest="aligner_path", default="", help="(set path to the aligner of choice IF not included in global PATH variable")
myparser.add_argument("-dg", "--degap", action="store", dest="degap",choices=["none", "all", "flanking"], default="all", help="Specify if and which gaps to remove:.\n\t'none' keep all gapped positions in the final MLSA-alignment\n\t'all' remove ALL gapped postitions in the alignments\n\t'flanking' remove flanking gapped positions from all individual alignments\nDefault='all'")
myparser.add_argument("-s", "--silent", action="store_true", dest="no_verbose", help="non-verbose mode")
args=myparser.parse_args()

wstrings, estrings=[], []
temp_path, fasta_path, PO_file=args.temp_path, args.fasta_path, args.po_resultfile
alignmeth, aligner_path=args.alignmeth, args.aligner_path
degap=args.degap
keep_temp=args.keep_temp
outputfilename="concatenated_orthologs_"+os.path.basename(PO_file)+".fasta"
verbose=True
docontinue=True

def checkargs(args):
	global verbose
	if args.no_verbose:
		verbose=False
	if fasta_path!="" and (not os.path.exists(fasta_path) or not os.path.isdir(fasta_path)):
		datsANerror("ERROR: fasta_path: '"+fasta_path+" does not exist or is no directory!")
	if aligner_path!="":
		if os.path.exists(aligner_path) and os.path.isdir(aligner_path):
			if not os.path.exists(os.path.join(aligner_path, "muscle")) and not os.path.exists(os.path.join(aligner_path, "clustalw")):
				datsANerror("ERROR: can't find 'muscle' or 'clustalw' at "+aligner_path)
		else:
			datsANerror("ERROR: aligner_path: '"+aligner_path+"' does not exist or is not a directory!")
	if not os.path.exists(PO_file) or not os.path.isfile(PO_file):
		datsANerror("ERROR: cannot find proteinortho-resultfile: "+PO_file)
	if temp_path!="": 
		if not os.path.exists(temp_path) or not os.path.isdir(temp_path):
			if verbose: 
				print "Creating directory for temporary and intermediate result files: "+temp_path
			os.mkdir(temp_path)
		else:
			if verbose:
				print "Will store temporary and intermediate result files in "+temp_path
	
def datsAwarning(wstring):
	global wstrings
	wstrings.append(wstring)
	if verbose:
		print wstring
	
def datsANerror(estring):
	global estrings
	global docontinue
	estrings.append(estring)
	docontinue=False
	if verbose:
		print >>sys.stderr, estring

def read_PO_file(filename):
	open_PO_file=open(filename, 'r')
	firstline=True
	org_index=3 #default index of beginning of organism-result-columns in proteinortho4+5
	headers=[]
	MLSA_gene_counter=0
	MLSA_list=[]
	for line in open_PO_file:
		if firstline:
			if line.startswith("#species"):
				datsAwarning("WARNING! It seems you are using results from Proteinortho4. It's recommended to use Proteinortho5!\n\t(However, this annotation transfer should still work)")
			elif line.startswith("# Species"):
				if verbose:
					print "\nProteinortho-results are based on Proteinortho5 or later. Good."
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
			datsANerror("ERROR: something went wron when reading proteinortho-results. Not the same number of comparison genes for all Organisms!")
	if verbose:
		print "\nrecognized "+str(len(MLSA_list[0][1]))+" unambigious (=unique) 'Orthologeous Groups' (OGs) shared by the comparison organisms"
	return headers, MLSA_list

def read_fasta_seqs(headers, MLSA_list):
	record_dict=dict.fromkeys(headers, {})
	for org in record_dict:
		record_dict[org]={}
	for h in range(len(headers)):
		if MLSA_list[h][0]!=headers[h]:
			datsANerror("ERROR: MLSA_list and headers are not synchronized! WTF!?!?!?!?")
			break
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
		tempfilename=os.path.join(temp_path, "selection_OG_fasta_"+MLSA_list[orgindex][0])
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
		if alignmeth=="clustal":
			aligned_filelist.append(call_clustalw(uf))
		elif alignmeth=="muscle":
			aligned_filelist.append(call_muscle(uf))
	return aligned_filelist
			 
def call_clustalw(inputfile):
		outputfile=inputfile.replace("unaligned_temp_fasta_", "SINGLEalignment_CLUSTALW_temp_fasta_",1)
		clustalw_cline= ClustalwCommandline("clustalw", INFILE=inputfile, outfile=outputfile, type="PROTEIN", align=True, quiet=True, OUTORDER="INPUT")
		clustalw_cline()
		return outputfile

def call_muscle(inputfile):
		outputfile=inputfile.replace("unaligned_temp_fasta_", "SINGLEalignment_MUSCLE_temp_fasta_",1)
		muscle_cline= MuscleCommandline(input=inputfile, out=outputfile, quiet=True) #add 'stable=True' to the end of this list, if the stable-bug in muscle is fixed (remove the correct_for_muscle_bug() method in that case)
		muscle_cline()
		return outputfile

def read_alignments(input_filelist):
	alignmentlist=[]
	if alignmeth=="muscle":
		format="fasta"# just a workaround for this method
	else:
		format=alignmeth
	for alignmentfile in input_filelist:
		currentalignfile=open(alignmentfile, 'r')
		alignmentlist.append(AlignIO.read(currentalignfile, format))
		currentalignfile.close()
	if len(alignmentlist)!=len(input_filelist):
		datsANerror("ERROR: could not read in all temporary alignmentfiles!")
	return alignmentlist

def remove_gaps_from_alignment_borders(alignmentlist):
	global outputfilename
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

def remove_gaps_from_complete_alignments(alignmentlist):
	global outputfilename
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
			concatenated_alignment[index].id="concatenated_orthologs_"+headers[index]
	else:
		datsANerror("ERROR: Records in concatenated_alignment and headers are not synchronized!")
	return concatenated_alignment

def write_final_alignment(outputfilename, concatenated_alignment):
	outputfile=open(outputfilename, 'w')
	AlignIO.write(concatenated_alignment,outputfile, "fasta")
	outputfile.close() 

def write_logfile():
	logfile=open("PO_2_MLSA.log",'w')
	logfile.write("Warnings:\n")
	for w in wstrings:
		logfile.write(w+"\n")
	logfile.write("Errors:\n")
	for e in estrings:
		logfile.write(e+"\n")

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
	
#main body:
if verbose:
	print "\n ==PO_2_MLSA.py by John Vollmers==\n"
checkargs(args)
if verbose:
	print "keep temporary files: "+str(keep_temp)
	if not keep_temp:
		print "will delete all temporary files"
	print "using aligner '"+alignmeth+"'"
	print "remove gaps: "+degap
if docontinue:
	headers, MLSA_list=read_PO_file(PO_file)
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
					if degap=="all":
						if verbose:
							print "\n----------------------------------"
							print "Removing all gapped positions from all single alignments prior to concatenation"
						alignment_list=remove_gaps_from_complete_alignments(alignment_list)
					elif degap=="flanking":
						if verbose:
							print "\n----------------------------------"
							print "Removing only the flanking gapped positions from all single alignments prior to concatenation"
						alignment_list=remove_gaps_from_alignment_borders(alignment_list)
					elif verbose:
						print "\nLeaving the single alignments as they are (Not removing any gapped positions)"
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

if keep_temp!=True:
	print "\n----------------------------------"
	print "cleaning up..."
	for delfile in unaligned_filelist:
		os.remove(delfile)
	for delfile in aligned_filelist:
		os.remove(delfile)

if docontinue:
	print "=================================="
	print "the individual ortholog selections were stored as: "
	print seq_filelist
	print "FINISHED: final aligned sequence saved as: "+outputfilename
	print "Remember to set sequence-type as 'PROTEIN' and filetype as 'fasta_wgap' when loading the MLSA sequence into Arb!\n"

				
if len(wstrings)>0 or len(estrings)>0:
	if len(wstrings)>0:
		print "\nThere were Warnings!"
	if len(estrings)>0:
		print "\nThere were Errors!"
	print "See PO_2_MLSA.log for further details\n"
	write_logfile()
