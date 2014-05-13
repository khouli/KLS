#!/usr/bin/env python

import sys
import argparse
import gzip

def get_parser():
	parser = argparse.ArgumentParser(description="Add Rosetta scores to"\
			" b-factor and charge columns in PDB files.")
	parser.add_argument("--bfactor", "-b", type=str,
		help="Residue energy term to place in bfactor column.")
	parser.add_argument("--occupancy", "-q", type=str,
		help="Residue energy term to place in occupancy column.")
	parser.add_argument("--inplace", "-i", action="store_true",
		help="Modify the input PDB files. Don't create new ones.")
	parser.add_argument("--suffix", "-s", type=str,
		help="Suffix for output PDB files. Input 1a0o.pdb gives output"\
				" 1a0o<SUFFIX>.pdb.")
	parser.add_argument("--pdb_gz", "-z", action="store_true",
		help="Output is a gzipped pdb, just like Rosetta.")
	parser.add_argument("pdbfnames", nargs='+', type=str,
						help="Names of PDB files")
	return parser


def main(argv=None):

	if argv is None:
		aparser = get_parser()
		argv = aparser.parse_args()

	if argv.inplace and argv.suffix:
		print "Concurrent use of inplace and suffix options is invalid."
		exit(1)
	inplace = argv.inplace

	bfac = argv.bfactor
	occ = argv.occupancy
	if not bfac and not occ:
		bfac = "total"

	suffix = argv.suffix
	if suffix is None and not inplace:
		bfacStr = "_%s_in_b" % (bfac,) if bfac else ""
		occStr = "_%s_in_q" % (occ,) if occ else ""
		suffix = ''.join([bfacStr, occStr])
	
	for f in argv.pdbfnames:
		gzIn = (f[-3:] == ".gz")
		pdb = gzip.open(f) if gzIn else open(f)
		hasScores = False
		for l in pdb:
			if l.startswith("#BEGIN_POSE_ENERGIES_TABLE"):
				hasScores = True
				break
		if not hasScores:
			print "PDB file", f, "does not have score information. Skipping"
			continue
		# Keep and split this line
		#label fa_atr fa_rep fa_sol fa_intra_rep fa_elec pro_close hbond_sr_bb hbond_lr_bb hbond_bb_sc hbond_sc dslf_fa13 rama
		labels = next(pdb)
		assert(labels.startswith("label"))
		labels = labels.strip().split()
		bndx = labels.index(bfac) if bfac else None
		qndx = labels.index(occ) if occ else None
		# Skip these two lines
		#weights ...
		#pose ...
		next(pdb); next(pdb)
		indices = []
		outcols = []
		outcolfirst = 80
		outcolafter = 0
		#COLUMNS     DATA TYPE       CONTENTS
		#----------------------------------------------------------------------
		# 55 - 60    Real(6.2)       Occupancy.
		# 61 - 66    Real(6.2)       Temperature factor (Default = 0.0).
		if qndx:
			indices.append(qndx)
			outcolfirst = min(outcolfirst, 54)
			outcolafter = max(outcolafter, 60)
		if bndx:
			indices.append(bndx)
			outcolfirst = min(outcolfirst, 60)
			outcolafter = max(outcolafter, 66)
		scores = []
		for l in pdb:
			if l.startswith("#END_POSE_ENERGIES_TABLE"): break
			scores.append( ["%6.2f" % (float(l.strip().split()[x]),) for x in indices] )
		pdb.seek(0)
		lines = pdb.readlines()
		if inplace and gzIn: argv.pdb_gz = True
		gzSuffix = ".gz" if argv.pdb_gz else ""
		outname = f if inplace else ''.join([ f.partition('.')[0], suffix, ".pdb", gzSuffix])
		outf = gzip.open(outname, 'w') if argv.pdb_gz else open(outname, 'w')
		scIter = iter(scores)
		s = next(scIter)
		lastResi = lines[0][22:26]
		for l in lines:
			outline = l
			if l.startswith("ATOM"):
				resi = l[22:26]
				if resi != lastResi:
					try:
						s = next(scIter)
					except StopIteration:
						break
					lastResi = resi
				sstr = ''.join(s)
				outline = ''.join([l[:outcolfirst], sstr, l[outcolafter:]])
			outf.write(outline)
		for l in lines: outf.write(l)
		pdb.close()
		outf.close()

if __name__ == "__main__":
	sys.exit(main())
