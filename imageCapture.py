import os
import argparse
from pymol import cmd, util

def args():
	parser = argparse.ArgumentParser(description='Rotate, colour and capture image of protein structure.')
	parser.add_argument('-i', '--infile', type=str, help='Input pdb file', required=True)
	parser.add_argument('-o', '--outfile', type=str, help='Output png file', required=False)
	return parser.parse_args()

def capture(pdb, pngout):
	cmd.load(pdb)
	cmd.set("ray_opaque_background",1)
	cmd.bg_color("white")
	util.cbc()
	cmd.orient()
	cmd.rotate(axis="z", angle=-90)
	cmd.png(pngout) 
	print "image saved to", pngout

if __name__ == "__main__":
	args = args()
	pdb = args.infile
	if args.outfile:
		png = args.outfile
	else:
		png = os.path.splitext(args.infile)[0] + ".png"

	capture(pdb, png)
