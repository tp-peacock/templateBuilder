import pandas as pd
import math
import os
from Bio.PDB import *
import cleanTCR

def getSequence(chain):
	ppb = PPBuilder()
	polypeps = ppb.build_peptides(chain)
	return "".join(map(lambda x: x.get_sequence()._data,polypeps))

def checkChain(chain):
	if chain == None:
		return None  # if no chain, return none
	if type(chain) == float and math.isnan(chain):
		return None # if no chain, return none
	else:
		return chain[0] # if a string, return first letter in string (in case of case: chain == "C | C")

df = pd.DataFrame.from_csv('mhcsearch_summary.tsv', sep='\t')

df_out = pd.DataFrame(columns=["template_file", "alpha_seq",
	"beta_seq", "antigen_seq", "mhc_seq"])

struct_dir = "template_structures/imgt"
# keep track of analysed pdbs to remove duplicates in input file
analysed = []

for index, row in df.iterrows():
	if index in analysed:
		continue

	achain = checkChain(row.Achain)
	bchain = checkChain(row.Bchain)
	pchain = checkChain(row.antigen_chain)
	mchain1 = checkChain(row.mhc_chain1)
	mchain2 = checkChain(row.mhc_chain2)

	if achain == None or bchain == None or mchain1 == None:
		continue

	infile = os.path.abspath(struct_dir + os.sep + index + ".pdb")
	# clean tcr/pmhc pdb file
	cleanfile = os.path.abspath(struct_dir +"_clean"+ os.sep + index + "_clean.pdb")
	cleanTCR.clean(infile, cleanfile, achain, bchain, mchain1, mchain2, pchain)
	# extract sequences from clean file
	parser = PDBParser()
	model = parser.get_structure('pdb', cleanfile)[0]

	sequences = {}
	for chain in [achain, bchain, pchain, mchain1]:
		if chain:
			sequences[chain] = getSequence(model.child_dict[chain])
		else:
			sequences[chain] = ""

	# add seqs to output dataframe
	df_out.loc[index] = [cleanfile, sequences[achain], sequences[bchain], sequences[pchain], sequences[mchain1]]

	analysed.append(index)

df_out.to_csv("templates.csv", index_label ="pdb_code")
print "Template data saved to templates.csv"
