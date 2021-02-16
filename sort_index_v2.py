# Sort sam files, convert to bam

import pysam
from argparse import ArgumentParser

def main():

	# parser
	parser = ArgumentParser()
	parser.add_argument("--InputSam", action="store", dest="input_sam", help="Input SAM file", required=True)
	parser.add_argument("--OutputBam", action="store", dest="output_bam", help="Output BAM file", required=True)
	o = parser.parse_args()
	sortindex(o.input_sam, o.output_bam)


def sortindex(input_sam, output_bam):
	pysam.sort('-o', output_bam, input_sam)
	pysam.index(output_bam)

main()

