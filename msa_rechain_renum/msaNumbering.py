#!/usr/bin/env python

from __future__ import print_function
import sys
import os
from itertools import imap
from itertools import chain
import subprocess
import re
#import urllib2
from string import rjust
# This module is from Biopython but that takes a lot of installing so I just
# plucked out the one module I need
#from Bio import pairwise2
import pairwise2
import string

# Given PDB structures, generate an MSA from relevant chains and then rename
# and renumber all the relevant chains accoringly.

# A fasta is required to properly generate the MSA. Generating MSA straight from
# PDB files could leave ambiguities. Suppose every sequence contains two
# alanines in a row and one sequence has missing density at one of the alanines.
# Which alanine is missing? The first or second? Fastas don't have missing density.

# This approach requires that missing denisty be matched to the fasta. This is
# first attempted with REMARK 465. As a fall back Calpha distances are tested.
# The fasta can come from SEQRES annotation inside a PDB or by querying
# RCSB for a fasta file. Using SEQRES requires annotated PDBs and using
# RCSB requires pdbcodes in filenames. If your PDB is properly annotated, it
# also has the PDB code so the second method is less restrictive and is
# what this script uses.

longer_names={'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D',
              'CYS': 'C', 'GLU': 'E', 'GLN': 'Q', 'GLY': 'G',
              'HIS': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K',
              'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S',
              'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'}

# from clean_pdb.py in Rosetta/tools/protein_tools/scripts
modres={'0CS':'ALA', '1AB':'PRO', '1LU':'LEU', '1PA':'PHE', '1TQ':'TRP', '1TY':'TYR', '23F':'PHE', '23S':'TRP', '2BU':'ALA', '2ML':'LEU', '2MR':'ARG', '2MT':'PRO', '2OP':'ALA', '2TY':'TYR', '32S':'TRP', '32T':'TRP', '3AH':'HIS', '3MD':'ASP', '3TY':'TYR', '4DP':'TRP', '4F3':'ALA', '4FB':'PRO', '4FW':'TRP', '4HT':'TRP', '4IN':'TRP', '4PH':'PHE', '5CS':'CYS', '6CL':'LYS', '6CW':'TRP', 'A0A':'ASP', 'AA4':'ALA', 'AAR':'ARG', 'AB7':'GLU', 'ABA':'ALA', 'ACB':'ASP', 'ACL':'ARG', 'ACY':'GLY', 'AEI':'THR', 'AFA':'ASN', 'AGM':'ARG', 'AGT':'CYS', 'AHB':'ASN', 'AHO':'ALA', 'AHP':'ALA', 'AIB':'ALA', 'AKL':'ASP', 'ALA':'ALA', 'ALC':'ALA', 'ALG':'ARG', 'ALM':'ALA', 'ALN':'ALA', 'ALO':'THR', 'ALS':'ALA', 'ALT':'ALA', 'ALY':'LYS', 'AME':'MET', 'AP7':'ALA', 'APH':'ALA', 'API':'LYS', 'APK':'LYS', 'AR2':'ARG', 'AR4':'GLU', 'ARG':'ARG', 'ARM':'ARG', 'ARO':'ARG', 'ASA':'ASP', 'ASB':'ASP', 'ASI':'ASP', 'ASK':'ASP', 'ASL':'ASP', 'ASN':'ASN', 'ASP':'ASP', 'AYA':'ALA', 'AYG':'ALA', 'AZK':'LYS', 'B2A':'ALA', 'B2F':'PHE', 'B2I':'ILE', 'B2V':'VAL', 'B3A':'ALA', 'B3D':'ASP', 'B3E':'GLU', 'B3K':'LYS', 'B3S':'SER', 'B3X':'ASN', 'B3Y':'TYR', 'BAL':'ALA', 'BBC':'CYS', 'BCS':'CYS', 'BCX':'CYS', 'BFD':'ASP', 'BG1':'SER', 'BHD':'ASP', 'BIF':'PHE', 'BLE':'LEU', 'BLY':'LYS', 'BMT':'THR', 'BNN':'ALA', 'BOR':'ARG', 'BPE':'CYS', 'BTR':'TRP', 'BUC':'CYS', 'BUG':'LEU', 'C12':'ALA', 'C1X':'LYS', 'C3Y':'CYS', 'C5C':'CYS', 'C6C':'CYS', 'C99':'ALA', 'CAB':'ALA', 'CAF':'CYS', 'CAS':'CYS', 'CCS':'CYS', 'CGU':'GLU', 'CH6':'ALA', 'CH7':'ALA', 'CHG':'GLY', 'CHP':'GLY', 'CHS':'PHE', 'CIR':'ARG', 'CLB':'ALA', 'CLD':'ALA', 'CLE':'LEU', 'CLG':'LYS', 'CLH':'LYS', 'CLV':'ALA', 'CME':'CYS', 'CML':'CYS', 'CMT':'CYS', 'CQR':'ALA', 'CR2':'ALA', 'CR5':'ALA', 'CR7':'ALA', 'CR8':'ALA', 'CRK':'ALA', 'CRO':'THR', 'CRQ':'TYR', 'CRW':'ALA', 'CRX':'ALA', 'CS1':'CYS', 'CS3':'CYS', 'CS4':'CYS', 'CSA':'CYS', 'CSB':'CYS', 'CSD':'CYS', 'CSE':'CYS', 'CSI':'ALA', 'CSO':'CYS', 'CSR':'CYS', 'CSS':'CYS', 'CSU':'CYS', 'CSW':'CYS', 'CSX':'CYS', 'CSY':'ALA', 'CSZ':'CYS', 'CTH':'THR', 'CWR':'ALA', 'CXM':'MET', 'CY0':'CYS', 'CY1':'CYS', 'CY3':'CYS', 'CY4':'CYS', 'CY7':'CYS', 'CYD':'CYS', 'CYF':'CYS', 'CYG':'CYS', 'CYJ':'LYS', 'CYQ':'CYS', 'CYR':'CYS', 'CYS':'CYS', 'CZ2':'CYS', 'CZZ':'CYS', 'DA2':'ARG', 'DAB':'ALA', 'DAH':'PHE', 'DAL':'ALA', 'DAM':'ALA', 'DAR':'ARG', 'DAS':'ASP', 'DBU':'ALA', 'DBY':'TYR', 'DBZ':'ALA', 'DCL':'LEU', 'DCY':'CYS', 'DDE':'HIS', 'DGL':'GLU', 'DGN':'GLN', 'DHA':'ALA', 'DHI':'HIS', 'DHL':'SER', 'DIL':'ILE', 'DIV':'VAL', 'DLE':'LEU', 'DLS':'LYS', 'DLY':'LYS', 'DMH':'ASN', 'DMK':'ASP', 'DNE':'LEU', 'DNG':'LEU', 'DNL':'LYS', 'DNM':'LEU', 'DPH':'PHE', 'DPL':'PRO', 'DPN':'PHE', 'DPP':'ALA', 'DPQ':'TYR', 'DPR':'PRO', 'DSE':'SER', 'DSG':'ASN', 'DSN':'SER', 'DTH':'THR', 'DTR':'TRP', 'DTY':'TYR', 'DVA':'VAL', 'DYG':'ALA', 'DYS':'CYS', 'EFC':'CYS', 'ESB':'TYR', 'ESC':'MET', 'FCL':'PHE', 'FGL':'ALA', 'FGP':'SER', 'FHL':'LYS', 'FLE':'LEU', 'FLT':'TYR', 'FME':'MET', 'FOE':'CYS', 'FOG':'PHE', 'FOR':'MET', 'FRF':'PHE', 'FTR':'TRP', 'FTY':'TYR', 'GHG':'GLN', 'GHP':'GLY', 'GL3':'GLY', 'GLH':'GLN', 'GLN':'GLN', 'GLU':'GLU', 'GLY':'GLY', 'GLZ':'GLY', 'GMA':'GLU', 'GMU':'ALA', 'GPL':'LYS', 'GT9':'CYS', 'GVL':'SER', 'GYC':'CYS', 'GYS':'GLY', 'H5M':'PRO', 'HHK':'ALA', 'HIA':'HIS', 'HIC':'HIS', 'HIP':'HIS', 'HIQ':'HIS', 'HIS':'HIS', 'HLU':'LEU', 'HMF':'ALA', 'HMR':'ARG', 'HPE':'PHE', 'HPH':'PHE', 'HPQ':'PHE', 'HRG':'ARG', 'HSE':'SER', 'HSL':'SER', 'HSO':'HIS', 'HTI':'CYS', 'HTR':'TRP', 'HY3':'PRO', 'HYP':'PRO', 'IAM':'ALA', 'IAS':'ASP', 'IGL':'ALA', 'IIL':'ILE', 'ILE':'ILE', 'ILG':'GLU', 'ILX':'ILE', 'IML':'ILE', 'IPG':'GLY', 'IT1':'LYS', 'IYR':'TYR', 'KCX':'LYS', 'KGC':'LYS', 'KOR':'CYS', 'KST':'LYS', 'KYN':'ALA', 'LA2':'LYS', 'LAL':'ALA', 'LCK':'LYS', 'LCX':'LYS', 'LDH':'LYS', 'LED':'LEU', 'LEF':'LEU', 'LET':'LYS', 'LEU':'LEU', 'LLP':'LYS', 'LLY':'LYS', 'LME':'GLU', 'LNT':'LEU', 'LPD':'PRO', 'LSO':'LYS', 'LYM':'LYS', 'LYN':'LYS', 'LYP':'LYS', 'LYR':'LYS', 'LYS':'LYS', 'LYX':'LYS', 'LYZ':'LYS', 'M0H':'CYS', 'M2L':'LYS', 'M3L':'LYS', 'MAA':'ALA', 'MAI':'ARG', 'MBQ':'TYR', 'MC1':'SER', 'MCL':'LYS', 'MCS':'CYS', 'MDO':'ALA', 'MEA':'PHE', 'MEG':'GLU', 'MEN':'ASN', 'MET':'MET', 'MEU':'GLY', 'MFC':'ALA', 'MGG':'ARG', 'MGN':'GLN', 'MHL':'LEU', 'MHO':'MET', 'MHS':'HIS', 'MIS':'SER', 'MLE':'LEU', 'MLL':'LEU', 'MLY':'LYS', 'MLZ':'LYS', 'MME':'MET', 'MNL':'LEU', 'MNV':'VAL', 'MPQ':'GLY', 'MSA':'GLY', 'MSE':'MET', 'MSO':'MET', 'MTY':'PHE', 'MVA':'VAL', 'N10':'SER', 'NAL':'ALA', 'NAM':'ALA', 'NBQ':'TYR', 'NC1':'SER', 'NCB':'ALA', 'NEP':'HIS', 'NFA':'PHE', 'NIY':'TYR', 'NLE':'LEU', 'NLN':'LEU', 'NLO':'LEU', 'NMC':'GLY', 'NMM':'ARG', 'NPH':'CYS', 'NRQ':'ALA', 'NVA':'VAL', 'NYC':'ALA', 'NYS':'CYS', 'NZH':'HIS', 'OAS':'SER', 'OBS':'LYS', 'OCS':'CYS', 'OCY':'CYS', 'OHI':'HIS', 'OHS':'ASP', 'OLT':'THR', 'OMT':'MET', 'OPR':'ARG', 'ORN':'ALA', 'ORQ':'ARG', 'OSE':'SER', 'OTY':'TYR', 'OXX':'ASP', 'P1L':'CYS', 'P2Y':'PRO', 'PAQ':'TYR', 'PAT':'TRP', 'PBB':'CYS', 'PBF':'PHE', 'PCA':'PRO', 'PCS':'PHE', 'PEC':'CYS', 'PF5':'PHE', 'PFF':'PHE', 'PG1':'SER', 'PG9':'GLY', 'PHA':'PHE', 'PHD':'ASP', 'PHE':'PHE', 'PHI':'PHE', 'PHL':'PHE', 'PHM':'PHE', 'PIA':'ALA', 'PLE':'LEU', 'PM3':'PHE', 'POM':'PRO', 'PPH':'LEU', 'PPN':'PHE', 'PR3':'CYS', 'PRO':'PRO', 'PRQ':'PHE', 'PRR':'ALA', 'PRS':'PRO', 'PSA':'PHE', 'PSH':'HIS', 'PTH':'TYR', 'PTM':'TYR', 'PTR':'TYR', 'PYA':'ALA', 'PYC':'ALA', 'PYR':'SER', 'PYT':'ALA', 'PYX':'CYS', 'R1A':'CYS', 'R1B':'CYS', 'R1F':'CYS', 'R7A':'CYS', 'RC7':'ALA', 'RCY':'CYS', 'S1H':'SER', 'SAC':'SER', 'SAH':'CYS', 'SAR':'GLY', 'SBD':'SER', 'SBG':'SER', 'SBL':'SER', 'SC2':'CYS', 'SCH':'CYS', 'SCS':'CYS', 'SCY':'CYS', 'SDP':'SER', 'SEB':'SER', 'SEC':'ALA', 'SEL':'SER', 'SEP':'SER', 'SER':'SER', 'SET':'SER', 'SGB':'SER', 'SGR':'SER', 'SHC':'CYS', 'SHP':'GLY', 'SIC':'ALA', 'SLZ':'LYS', 'SMC':'CYS', 'SME':'MET', 'SMF':'PHE', 'SNC':'CYS', 'SNN':'ASP', 'SOC':'CYS', 'SOY':'SER', 'SUI':'ALA', 'SUN':'SER', 'SVA':'SER', 'SVV':'SER', 'SVX':'SER', 'SVY':'SER', 'SVZ':'SER', 'SXE':'SER', 'TBG':'GLY', 'TBM':'THR', 'TCQ':'TYR', 'TEE':'CYS', 'TH5':'THR', 'THC':'THR', 'THR':'THR', 'TIH':'ALA', 'TMD':'THR', 'TNB':'CYS', 'TOX':'TRP', 'TPL':'TRP', 'TPO':'THR', 'TPQ':'ALA', 'TQQ':'TRP', 'TRF':'TRP', 'TRN':'TRP', 'TRO':'TRP', 'TRP':'TRP', 'TRQ':'TRP', 'TRW':'TRP', 'TRX':'TRP', 'TTQ':'TRP', 'TTS':'TYR', 'TY2':'TYR', 'TY3':'TYR', 'TYB':'TYR', 'TYC':'TYR', 'TYI':'TYR', 'TYN':'TYR', 'TYO':'TYR', 'TYQ':'TYR', 'TYR':'TYR', 'TYS':'TYR', 'TYT':'TYR', 'TYX':'CYS', 'TYY':'TYR', 'TYZ':'ARG', 'UMA':'ALA', 'VAD':'VAL', 'VAF':'VAL', 'VAL':'VAL', 'VDL':'VAL', 'VLL':'VAL', 'VME':'VAL', 'X9Q':'ALA', 'XX1':'LYS', 'XXY':'ALA', 'XYG':'ALA', 'YCM':'CYS', 'YOF':'TYR'}
#'

# The sequence of 1UBQ.pdb is taken as the true sequence of ubiquitin.
# First, use this sequence to find the relevant chains in the given PDBs.
oneubq_seq = "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG"
start_resi = 1
#cutoff = 0.7 * len(oneubq_seq)
cutoff = 0.25 * len(oneubq_seq)
# cutoff = for 1ubq seq, 0.7*len(seq) = 53.2

# alignment function that will be used
# global alignment
# align score = 1 pt for match, -0.5 for mismatch, -0.25 for gap, -0.1 for gap extension
# With these numbers most ubiquitins will score ~ 76.0, non-ubiquitins ~ 6
alignScore = lambda x: max( a[2] for a in pairwise2.align.globalms(oneubq_seq, x, 1.0, -0.4, -0.1, -0.1))
# Explanation: globalxx(seq1, seq2) does a pairwise alignment and returns a list of
# tuples in which element 2 of each tuple is the alignment score. globalms is like
# globalxx but takes parameters for how to score the alignment.

def fastaSieve(fname, filterFunc, cutoff):
	atleast_1_match = False
	matches = []
	with open(fname) as f:
		nameline = None
		pdb = None
		ch = None
		seq = ""
		for l in f:
			if l.startswith('>'):
				if nameline:
					alignScore = filterFunc(seq)
					if alignScore > cutoff:
						atleast_1_match = True
						# TODO use a namedtuple
						matches.append( nameline + seq )
				nameline = l
				seq = ""
				#pdb, ch = l.split('|')[0][1:].split(':')
				continue
			#elif not pdb or not ch:
			elif not nameline:
				print("Input file", fname, "contains an unlabeled sequence.")
				exit(0)
			seq = seq + l.strip()
		alignScore = filterFunc(seq)
		if alignScore > cutoff:
			matches.append( nameline + seq )
			atleast_1_match = True
	if not atleast_1_match:
		print("Input file \"" + fname + "\" does not contain any sequence with "
			"sufficient alignment score to 1ubq.")
	return( matches )

if (len(sys.argv) < 2):
	print("Need fastas.")
	exit(0)
if (len(sys.argv) < 3):
	print("Need more than 1 fastas. *Multiple* Sequence Alignment needs multiple sequences.")
	exit(0)

infiles = sys.argv[1:]

# For every input file check for sequence matches
fnames = sys.argv[1:]

matches = chain.from_iterable( fastaSieve(fname, alignScore, cutoff) for fname in fnames )

# the matches objects has contents like this:
#[
#">2LVO:A|PDBID|CHAIN|SEQUENCE",
#"MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG",
#">1P3Q:U|PDBID|CHAIN|SEQUENCE",
#"MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG",
#">1P3Q:V|PDBID|CHAIN|SEQUENCE",
#"MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG"
#]

# Add the template sequence, oneubq, to include it in the MSA.
# After MSA, the renumbering will be such that oneubq starts with 1.
# All leading gaps in the MSA of the oneubq sequence will have negative indices.

oneubq_name = "ONEUBQ_TEMPLATE"
oneubq_chain = "U"
matches = [ ''.join([">", oneubq_name, ":", oneubq_chain, "|PDBID|CHAIN|SEQUENCE"]), oneubq_seq ] + list(matches)

devnull = open(os.devnull, 'w')
muscle = subprocess.Popen(["muscle", "-clw"], stdin=subprocess.PIPE, \
		stdout=subprocess.PIPE, stderr=devnull)
musclein = '\n'.join(matches)
muscle.stdin.write(musclein)
muscle.stdin.close()
muscle.wait()

pdbAlignments = {}
conservations = ""

muscleMSA = open("muscleMSA.txt", "wt")
for l in muscle.stdout.readlines():
	muscleMSA.write(l)
	if l.startswith("MUSCLE"):
		continue
	l = l.strip()
	if not l:
		continue
	# save conservation lines when we find them and continue
	if l.strip().startswith( ("*", ":", ".") ):
		conservations += l
		continue
	l = l.split('|')
	pdbcode, chain = l[0].split(':')
	alignment = l[-1].split()[-1]
	# previously seen neither pdbcode nor chain
	if pdbcode not in pdbAlignments:
		pdbAlignments[ pdbcode ] = {chain : alignment}
		continue
	# previously seen pdbcode but not chain
	elif chain not in pdbAlignments[pdbcode]:
		pdbAlignments[pdbcode][chain] = alignment
		continue
	# previously seen pdbcode and chain
	elif chain in pdbAlignments[pdbcode]:
		pdbAlignments[pdbcode][chain] += alignment
		continue
	# wtf? This shouldn't happen, print out some stuff for debugging and stop
	print(pdbAlignments)
	print("pdb:", pdbcode)
	print(len(pdbcode))
	print("chain:", chain)
	print(len(chain))
	assert(False)
muscleMSA.close()

# Clustal W alignment symbols key
# *  -- all residues or nucleotides in that column are identical
#         :  -- conserved substitutions have been observed
#         .  -- semi-conserved substitutions have been observed
#            -- no match.

# This script numbers ':', '.', and ' ' positions the same. Identical
# sequence numbers generated by this script do not necessarily
# indicate perfect conservation.

def parseMissDens(pdblines):
	missDens = {}
	pdblines = iter(pdblines)
	for l in pdblines:
		if re.match("REMARK 465   [M ] RES C SSSEQI", l):
			break
	for l in pdblines:
		if not l.startswith("REMARK 465"):
			break
		l = l.strip().split()
		if len(l) != 5:
			break
		ch = l[3]
		resNum = int(l[4])
		if ch in missDens:
			missDens[ch].append(resNum)
		else:
			missDens[ch] = [resNum]
	return(missDens)


# handle oneubq numbering here
oneubq_alignment = pdbAlignments[oneubq_name][oneubq_chain]
nongap = re.search(r'[^-]', oneubq_alignment)
# this index should be resi #1 for all sequences
oneubq_leading_gaps = nongap.start()
global_start_resi = 1 - oneubq_leading_gaps

#start_resi = global_start_resi
#for i in oneubq_alignment:
#	print(''.join([str(start_resi), ": ", i]))
#	start_resi += 1
#print(pdbAlignments)
#exit(0)

# Needs chainIDmap and returns (chainIDmap, renumMap)
def genRenumMap( resNamesNums, alignment, missDens):
	nongap = re.search(r'[^-]', alignment)
	# If entire seq is gap, nongap == None
	assert(nongap != None)
	# align_i is the index where we are in the alignment string
	align_i  = nongap.start()
	# numATOM is the name we're writing to the PDB
	numATOM = nongap.start() + global_start_resi

	#print("alignment: " + alignment)
	#print("nongap.start(): " + str(nongap.start()))
	#print("global_start_resi: " + str(global_start_resi))
	#exit(0)

	# very first HETATM number
	numHETATM = len(alignment)

	renumMap = {}

	# case 1:
	#	name3 matches alignment and resNum not in missDens
	#
	#	increment numATOM and alignment, add renumMap entry, advance to next pdbRes
	#
	# case2:
	#	missing density: name3 doesn't match alignment and resNum in missDens
	#
	#	increment numATOM and alignment, don't advance pdbRes
	#
	# case3:
	#	alignment gap: and resNum NOT in missDens
	#
	#	increment numATOM and alignment, don't advance pdbRes
	#
	# otherwise raise an error

	#print("resNamesNums:", resNamesNums)
	missDens_i = 0
	firstResNum = resNamesNums[0][1]
	skip = 0
	for m in missDens:
		if firstResNum > m:
			missDens_i += 1
		else:
			break
	
	align_i += missDens_i

	#print("alignment:", alignment)
	#print("missDens:", missDens)
	#print("resNamesNums:", [ r[0] for r in resNamesNums ])
	lastResNum = resNamesNums[0][1] - 1
	missDensResNum = missDens[missDens_i] if missDens_i < len(missDens) else lastResNum
	for rNN in resNamesNums:
		#print("rNN:", rNN)
		name3, resNum = rNN
		assert( resNum not in renumMap )
		caa = name3 in longer_names
		name3 = modres[name3] if name3 in modres else name3
		#print("name3:", name3)
		#print("name3 in longer_names:", name3 in longer_names)
		hetatm = True if name3 not in longer_names else False
		if hetatm:
			#print("HETATM numbering")
			renumMap[resNum] = numHETATM
			numHETATM += 1
			continue
		# ATOM records beyond this point
		name1 = longer_names[name3]
		while( missDensResNum > lastResNum and missDensResNum < resNum ):
			#print("skipping missing dens")
			#print("alignment[align_i]:", alignment[align_i])
			numATOM += 1
			align_i += 1
			missDens_i += 1
			missDensResNum = missDens[missDens_i] if missDens_i < len(missDens) \
					else lastResNum
		#print("name1:", name1)
		#print("align:", alignment[align_i])
		#print("align_i:", align_i)
		# if NCAA and alignment mismatch, assume HETATM
		if not caa and alignment[align_i] == '-':
			#print("Not CAA and no alignment, continue")
			renumMap[resNum] = numATOM
			numATOM += 1
			continue
		# otherwise this may be an alignment gap
		while(alignment[align_i] == '-'):
			#print("Skipping sequence gap")
			numATOM += 1
			align_i += 1
			# TODO: fix the handling of noncanonical aminoacids
			#print("CAA:", caa)
			#print("align_i == len(alignment)", align_i == len(alignment))
			if align_i == len(alignment) - 1:
				break
		#python list index != output resi, align_i anymore
		#assert(name1 == alignment[align_i - global_start_resi])

		#print("Writing renumMap entry")
		renumMap[resNum] = numATOM
		#print(renumMap)
		align_i += 1
		numATOM += 1
		lastResNum = resNum

	#print(renumMap)

	return(renumMap)


def writeModel(modelBuffer, chainAlignments, chainIDmap, missDens, pdbout):

	chainResNamesNums = {}
	#print("Collecting res names and numbers in each chain", file=sys.stderr)
	for l in modelBuffer:
		ch = l[21]
		if ch not in chainAlignments:
			continue
		name3 = l[17:20]
		resNum = int( l[22:26] )
		resNameNum = (name3, resNum)
		if ch in chainResNamesNums:
			if resNameNum not in chainResNamesNums[ch]:
				chainResNamesNums[ch].append( resNameNum )
		else:
			chainResNamesNums[ch] = [ resNameNum ]
	
	renumMaps = {}
	for c in chainAlignments:
		#print("on chain:", c)
		chainMissDens = missDens[c] if c in missDens else []
		renumMaps[c] = genRenumMap(chainResNamesNums[c], chainAlignments[c], chainMissDens)
	
	ubqChIDs = list("UWXYZuwxyz")
	for c in chainAlignments:
		if c in chainIDmap:
			continue
		for u in ubqChIDs:
			if u not in chainIDmap.values():
				chainIDmap[c] = u
				break
		chainIDmap[c] = u
		assert(c in chainIDmap)

	#print("chainIDmap:", chainIDmap)
	
	# every chain we need to edit is now in chainIDmap and renumMaps
	for l in modelBuffer:
		ch = l[21]
		if ch not in chainIDmap:
			pdbout.write(l)
			continue
		newCh = chainIDmap[ch]
		resNum = int(l[22:26])
		renum = rjust(str(renumMaps[ch][resNum]), 4)
		pdbout.write( ''.join( (l[:21], newCh, renum, l[26:]) ) )
	
	return(chainIDmap)

for pdb in pdbAlignments:
	if(pdb == oneubq_name):
		continue
	print("Writing renumbered pdb:" + pdb)

	# Assume pdbs have the same naming as given inside the fastas
	pdbname = pdb + ".pdb"
	#print("pdb:", pdbname)
	pdblines = open(pdbname).readlines()
	missDens = parseMissDens(pdblines)
	pdboutname = pdb + "_MSArenum.pdb"
	pdbout = open(pdboutname, 'wt')
	chainAlignments = pdbAlignments[pdb]
	modelBuffer = []
	# chainIDs may need to be coordinated across multiple models,
	# store the renamings
	chainIDmap = {}
	for l in pdblines:
		# If this is a new model forget all you thought you knew
		if re.match(r"ATOM  |HETATM|ANISOU|TER   ", l):
			modelBuffer.append(l)
		else:
			if modelBuffer:
				#print("Writing modelBuffer", file=sys.stderr)
				chainIDmap = writeModel(modelBuffer, chainAlignments, \
						chainIDmap, missDens, pdbout)
				modelBuffer = []
			pdbout.write(l)
	pdbout.close()
