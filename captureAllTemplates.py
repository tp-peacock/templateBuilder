import os
import argparse
from pymol import cmd, util

if __name__ == "__main__":

	for pdb in os.listdir("template_structures/imgt"):
		pdbcode = os.path.splitext(pdb)[0]

		raw = "template_structures/imgt/"+pdb
		clean = "template_structures/imgt_clean/"+pdbcode+"_clean.pdb"
		raw_pymol = os.path.splitext(os.path.basename(raw))[0]
		clean_pymol = os.path.splitext(os.path.basename(clean))[0]
		out = pdbcode+"_raw_to_clean.png"

		cmd.reinitialize()
		cmd.load(raw)
		cmd.load(clean)
		cmd.set("ray_opaque_background",1)
		cmd.bg_color("white")
		cmd.orient(clean_pymol)
		cmd.rotate(axis="z", angle=90)
		cmd.zoom(raw_pymol)
		cmd.color("grey70", raw_pymol)
		util.cbc(clean_pymol)

		cmd.png("raw_to_clean_images/"+out) 
		print "image saved to", out
		
