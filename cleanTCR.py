from Bio.PDB import *
import argparse
import os

def args():
	parser = argparse.ArgumentParser(description='Clean TCR structure')
	parser.add_argument('-i', '--infile', type=str, help='name of input raw pdb file', required=True)
	parser.add_argument('-a', '--achain', type=str, help='alpha chain label', required=True)
	parser.add_argument('-b', '--bchain', type=str, help='beta chain label', required=True)
	parser.add_argument('-m1', '--mchain1', type=str, help='mhc chain 1 label', required=True)
	parser.add_argument('-m2', '--mchain2', type=str, help='mhc chain 2 label', required=True)
	parser.add_argument('-p', '--pchain', type=str, help='peptide/antigen chain label', required=True)
	parser.add_argument('-o', '--outfile', type=str, help='name of output clean pdb file', required=False)
	return parser.parse_args()

# truncate TCR to variable domain only (alpha)
# truncate TCR to variable domain only (beta)
# truncate MHC to g-domain only (correct terminology? Hold for all MHC?)
# remove water and other non-protein molecules 
# for model with multiple TCRs etc, keep only first set
# keep old header

class CleanPDB(Select):

	def __init__(self, achain, bchain, mchain1, mchain2, pchain):
		self.achain  = achain
		self.bchain  = bchain
		self.mchain1 = mchain1
		self.mchain2 = mchain2
		self.pchain  = pchain

	def accept_model(self, model):
		# keep only the first model if more than one model exists
		if model.id == 0:
			return 1
		else:
			return 0

	def accept_chain(self, chain):
		# reject anything not labelled as first TCRpMHC complex
		# for now, reject mhc chain 2
		if chain.id in [self.achain, self.bchain, self.mchain1, self.pchain]:
			return 1
		else:
			return 0

	def accept_residue(self, residue):
		# assumes chains are already imgt numbered
		# truncate a and b chains to variable region (imgt numbers 1 to 127)
		if residue.parent.id == self.achain or residue.parent.id == self.bchain:
			if residue.id[1] < 1 or residue.id[1] > 127:
				return 0

		# assumes chains are already imgt numbered
		# truncate mhc chains to G domain only (imgt numbers 1 to 92)
		if residue.parent.id == self.mchain1:
			if residue.id[1] not in range(1,91) and residue.id[1] not in range(1001,1091):
				return 0
		
		# accept only amino acids, reject water molecules etc
		if is_aa(residue):
			return 1
		else:
			return 0

# tidy up to make this general

def clean(pdb, pdbout, achain, bchain, mchain1, mchain2, pchain):
	parser = PDBParser()
	structure = parser.get_structure('pdb', pdb)
	io = PDBIO()
	io.set_structure(structure)
	clean_pdb = CleanPDB(achain = achain, bchain = bchain, 
		mchain1 = mchain1, mchain2 = mchain2, pchain = pchain)
	io.save(pdbout, clean_pdb)
	print "Clean structure saved to", pdbout


if __name__ == "__main__":

	args = args()
	if not args.outfile:
		args.outfile = "clean_" + os.path.basename(args.infile)

	clean(args.infile, args.outfile, args.achain, args.bchain, 
		args.mchain1, args.mchain2, args.pchain)


