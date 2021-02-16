# Put UMI in read1 name; r1 leave "cutlen" nt, r2 discard

from argparse import ArgumentParser

def main():

	# parser

	parser = ArgumentParser()

	parser.add_argument("--FastqR1", action="store", dest="input_name1", help="Input FASTQ R1 file", required=True)

	parser.add_argument("--FastqR2",  action="store",  dest="input_name2", help="Input FASTQ R2 file", required=True)

	parser.add_argument("--OutputName",  action="store",  dest="output_name1", help="Output trimmed FASTQ file", required=True)

	o = parser.parse_args()

	clean_fastq_wPrimer(o.input_name1, o.input_name2, o.output_name1)


def reverse_complement(seq):
	complement = {'A':'T', 'T':'A', 'C':'G', 'G':'C', 'a':'t', 't':'a', 'c':'g', 'g':'c'}
	bases = list(seq)
	bases = [complement[base] for base in bases]
	return ''.join(bases[::-1])


def clean_fastq_wPrimer(input_name1, input_name2, output_name1):
	from itertools import izip
	with open(input_name1, "r") as r1, open(input_name2, "r") as r2, open(output_name1, "w") as w1:
		counter = 0
		last = False
		bunchmax = 1000000
		bunchsize = 0
		bunchnum = 0
		bunch1 = []
		#bunch2 = []
		minlen = 60
		cutlen = 60
		for l1, l2 in izip(r1, r2):
			if counter%4 == 0: #name
				#create two read elements and 
				read1 = ['','','','']
				#read2 = ['','','','']
				read1[0] = l1
				#read2[0] = l2
			elif counter%4 == 1:#seq
				seq1 = l1[:cutlen] +'\n'
				umi_seq1 = l2[4:19] #UMI
				to_write = 1
				if 'N'in seq1 or 'N' in umi_seq1:
					to_write = 0
				else:
					read1[1] = seq1
					spacepos = read1[0].find(' ')
					read1[0] = read1[0][:spacepos]+'UMI'+umi_seq1+'\n'
					#read2[1] = new_seq2
			elif counter%4 == 2:#+
				read1[2] = l1
				#read2[2] = l2
			elif counter%4 == 3:#Qscore
				if to_write == 1:
					read1[3] = l1[:cutlen] +'\n'
					bunch1.extend(read1)
					#bunch2.extend(read2)
					bunchsize += 1

			if bunchsize > bunchmax:
				w1.writelines(bunch1)
				#w2.writelines(bunch2)
				bunchsize = 0 
				bunchnum += 1
				bunch1 = []
				#bunch2 = []
			counter += 1 
		w1.writelines(bunch1)
		#w2.writelines(bunch2)
		print counter/4
		print bunchnum*bunchmax+bunchsize
		print ''


main()
