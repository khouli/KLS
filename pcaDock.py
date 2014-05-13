#!/usr/bin/env python

import sys
import argparse
from string import ljust, rjust
import pandas as pd
import numpy as np
from numpy import NaN
import scipy as sp
import scipy.spatial
from sklearn.decomposition import PCA
from math import cos, sin, pi

#class Structure:
#	sequences
#	models
#	chains
#	atoms

def genATOMlines(lines, hetatm=False):
	#http://www.wwpdb.org/documentation/format33/sect9.html#ATOM
	#COLUMNS        DATA  TYPE    FIELD        DEFINITION
	#-------------------------------------------------------------------------------------
	# 1 -  6        Record name   "ATOM  "
	# 7 - 11        Integer       serial       Atom  serial number.
	#13 - 16        Atom          name         Atom name.
	#17             Character     altLoc       Alternate location indicator.
	#18 - 20        Residue name  resName      Residue name.
	#22             Character     chainID      Chain identifier.
	#23 - 26        Integer       resSeq       Residue sequence number.
	#27             AChar         iCode        Code for insertion of residues.
	#31 - 38        Real(8.3)     x            Orthogonal coordinates for X in Angstroms.
	#39 - 46        Real(8.3)     y            Orthogonal coordinates for Y in Angstroms.
	#47 - 54        Real(8.3)     z            Orthogonal coordinates for Z in Angstroms.
	#55 - 60        Real(6.2)     occupancy    Occupancy.
	#61 - 66        Real(6.2)     tempFactor   Temperature  factor.
	#77 - 78        LString(2)    element      Element symbol, right-justified.
	#79 - 80        LString(2)    charge       Charge  on the atom.
	for line in lines:
		if line.startswith("ATOM  ") or (hetatm and line.startswith("HETATM")):
			record    = line[0:6].strip()
			atomNum   = np.int_(line[6:11])
			atomName  = line[12:16]
			altLoc    = line[16]
			resName   = line[17:20]
			chainID   = line[21]
			resNum    = np.int_(line[22:26])
			iCode     = line[26]
			x         = np.float_(line[30:38])
			y         = np.float_(line[38:46])
			z         = np.float_(line[46:54])
			occupancy = np.float_(line[54:60]) if line[54:60] else NaN
			bFactor   = np.float_(line[60:66]) if line[60:66] else NaN
			element   = line[76:78]
			charge    = line[78:80] #type is LString(2), so no cast
			yield (record, atomNum, atomName, altLoc, resName, chainID, \
				   resNum, iCode, x, y, z, occupancy, bFactor, element, charge)


def writePDB(struct, fname, writemode='w'):
	testout = open(fname, writemode)
	for v in struct.values:
		record, atomNo, atomName, altLoc, resName, chainID, resNum, iCode, \
				x, y, z, occupancy, bFactor, element, charge = v
		# 1 -  6        Record name   "ATOM  "
		record = ljust(record, 6)
		# 7 - 11        Integer       serial       Atom  serial number.
		atomNo = rjust(str(atomNo), 5)
		#13 - 16        Atom          name         Atom name.
		#atomName
		#17             Character     altLoc       Alternate location indicator.
		#altLoc
		#18 - 20        Residue name  resName      Residue name.
		#resName
		#22             Character     chainID      Chain identifier.
		#chainID
		#23 - 26        Integer       resSeq       Residue sequence number.
		resNum = rjust(str(resNum), 4)
		#27             AChar         iCode        Code for insertion of residues.
		#iCode
		#31 - 38        Real(8.3)     x            Orthogonal coordinates for X in Angstroms.
		#39 - 46        Real(8.3)     y            Orthogonal coordinates for Y in Angstroms.
		#47 - 54        Real(8.3)     z            Orthogonal coordinates for Z in Angstroms.
		xyz = "%8.3f%8.3f%8.3f" % (x, y, z)
		#55 - 60        Real(6.2)     occupancy    Occupancy.
		occupancy = "%6.2f" % (occupancy,)
		#61 - 66        Real(6.2)     tempFactor   Temperature  factor.
		bFactor = "%6.2f" % (bFactor,)
		#77 - 78        LString(2)    element      Element symbol, right-justified.
		#element
		#79 - 80        LString(2)    charge       Charge  on the atom.
		#charge

		testout.write( ''.join((record,atomNo," ",atomName,altLoc,resName,' ',chainID,resNum,\
				      ' '*3,iCode,xyz,occupancy,bFactor,' '*10,element,charge,'\n')) )
	testout.write( "TER" + ' '*77 + '\n' )
	testout.close()


def readPDB(fname):
	with open(fname) as f:
		try:
			lines = f.readlines()
		except IOError:
			print "Error opening file."
		atoms = genATOMlines(lines)
		struct = pd.DataFrame.from_records(atoms,
				columns=["record", "atomno", "atomName", "altLoc", "resName", \
						"chainID", "resNum", "iCode", "x", "y", "z", \
						"occupancy", "bFactor", "element", "charge"])
	return(struct)


# return a ndarray of point COM
# not weighted by element masses
def xyzCOM(xyz):
	#xyz = struct.loc[:,["x","y","z"]]
	return(xyz.mean(axis=0).values)


def xyzsets_min_dist(xyz1, xyz2):
	#xyz = ["x","y","z"]
	#chain_cdist = sp.spatial.distance.cdist( atoms1.loc[:,xyz], atoms2.loc[:,xyz] )
	chain_cdist = sp.spatial.distance.cdist( xyz1, xyz2 )
	cdist_argmin = np.argmin(chain_cdist)
	struct1_atom, struct2_atom = np.unravel_index(cdist_argmin, chain_cdist.shape)
	cdist_min = chain_cdist[struct1_atom, struct2_atom]
	return(struct1_atom, struct2_atom, cdist_min)


def xyzsets_mean_dist(xyz1, xyz2):
	#xyz = ["x","y","z"]
	#chain_cdist = sp.spatial.distance.cdist( atoms1.loc[:,xyz], atoms2.loc[:,xyz] )
	chain_cdist = sp.spatial.distance.cdist( xyz1, xyz2 )
	return(np.mean(chain_cdist))


def rot2align(mobileFrame, fixedFrame):
	V, S, Wt = np.linalg.svd( np.dot(np.transpose(mobileFrame), fixedFrame) )
	reflect = float(np.linalg.det(V) * np.linalg.det(Wt))
	if abs(reflect + 1) < 1e-6:
		S[-1]   = -S[-1]
		V[:,-1] = -V[:,-1]
	U = np.dot(V, Wt)
	return(U)


def slideAlongVec(vec, mobile, fixed):
	mobile = mobile + 500*vec
	dist = xyzsets_min_dist(fixed, mobile)[2]
	#print "dist:", dist
	while(dist > 2.8):
		#if dist > 1000: break
		# 0.3 seems stable but it still possible to
		# "miss" and loop forever
		mobile = mobile - max(0.25*dist, 0.25)*vec
		dist = xyzsets_min_dist(fixed, mobile)[2]
		#print "dist:", dist
	return(mobile)


# this works but don't use this, just have docking_protocol do the spining
def rotation(theta):
   tx,ty,tz = theta

   Rx = np.array([[1,0,0], [0, cos(tx), -sin(tx)], [0, sin(tx), cos(tx)]])
   Ry = np.array([[cos(ty), 0, -sin(ty)], [0, 1, 0], [sin(ty), 0, cos(ty)]])
   Rz = np.array([[cos(tz), -sin(tz), 0], [sin(tz), cos(tz), 0], [0,0,1]])

   return np.dot(Rx, np.dot(Ry, Rz))


def get_parser():
	parser = argparse.ArgumentParser(
		description="Generate an initial dock of two PDBs using principle axes.\n"
			"Example usage:\n./pcaDock.py -f 1ubq.pdb -m myUBD.pdb -v \"1, 0, 0\" -o boundUbq.pdb")
	parser.add_argument("--fixed", "-f", type=str, required=True,
		help="Structure to place at origin and then keep fixed")
	parser.add_argument("--mobile", "-m", type=str, required=True,
		help="PDB cotaining structure to move about fixed structure")
	parser.add_argument("--vector", "-v", type=str,
		help="Displacement vector to apply to COM of mobile structure."
			"Ex: -v \"1, 0, 0\"")
	parser.add_argument("--outfile", "-o", type=str, required=True,
		help="Output file name. Contains both input structures.")
	parser.add_argument("--show", "-s", action="store_true",
		help="Generate pdb with fixed structure and print COM"
			" displacement of mobile from that fixed structure.")
	return parser


def main(argv=None):

	if argv is None:
		aparser = get_parser()
		argv = aparser.parse_args()


	fname1 = argv.fixed
	fname2 = argv.mobile
	outname = argv.outfile
	
	struct1 = readPDB(fname1)
	struct2 = readPDB(fname2)
	
	cID1 = struct1.loc[1, "chainID"]
	cID2 = struct2.loc[1, "chainID"]
	if cID1 == cID2:
		delta_ord = 1 if cID2 not in ['z', 'Z', ' ', '9'] else -1
		cID2 = chr(ord(cID2) + delta_ord)
		struct2.loc[:,"chainID"] = cID2
	
	struct1_xyz = struct1.loc[:,["x","y","z"]]
	struct2_xyz = struct2.loc[:,["x","y","z"]]

	struct1_xyz = struct1_xyz - xyzCOM(struct1_xyz)

	if argv.show:
		print "Mobile structure COM xyz displacement:", xyzCOM(struct2_xyz)
		print "Displacement is relative to output PDB with fixed structure."
	struct2_xyz = struct2_xyz - xyzCOM(struct2_xyz)
	
	# Put the most siginificant principle axis along z, this way binding tends
	# to be in XY plane
	xyzFrame = np.array( [[0.0, 0.0, 1.0], 
			              [1.0, 0.0, 0.0],
			              [0.0, 1.0, 0.0]] )
	
	pca = PCA(n_components=3)
	pca.fit(struct1_xyz)
	struct1_paxis = pca.components_
	pca.fit(struct2_xyz)
	struct2_paxis = pca.components_
	
	U = rot2align(struct1_paxis, xyzFrame)
	struct1_xyz = np.dot( struct1_xyz, U )
	
	struct1.loc[:, ["x","y","z"]] = struct1_xyz
	# weirdness happened when I rotated the principle axes
	# whatev, just resolve
	pca.fit(struct1_xyz)
	struct1_paxis = pca.components_
	
	# change this so both are alinged to xyz
	U = rot2align(struct2_paxis, xyzFrame)
	struct2_xyz = np.dot( struct2_xyz, U )
	
	vec = np.array([[-1,0,0]]) # default value
	if argv.vector:
		vec = np.fromstring(argv.vector, dtype=np.float_, sep=',')
	
	struct2.loc[:, ["x","y","z"]] = slideAlongVec(vec, struct2_xyz, struct1_xyz)
	writePDB( struct1, outname, 'w' )
	if argv.show: return
	writePDB( struct2, outname, 'a' )

if __name__ == "__main__":
	sys.exit(main())
