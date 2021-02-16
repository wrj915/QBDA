#!/usr/bin/env python
# coding=utf-8


"""
Modified by PD XW 09/21/2020 added dynamic cutoff based on amplicon (not requiring same genotype), majority vote no threshold
UMIs with family size > max(3, ceil(mean(Top3_familysize)*0.05)) are kept
previous: UMIs with family size > 3 are kept


UMI_counter.py



Version 1 (June 6 2019) By Sherry Chen



UPDATED 1 (July 3 2019) By Sherry Chen



	1. fixed a mini bug on de_novo not handleing the last base difference

	2. update the maximum size of each pileupcolumn to 1M  (max max_depth = 1000000)

	3. added level2 files and final files



Updated ver 1_2 (Aug 9 2019) By Sherry Chen

	1. changed majority vote -> if wt exists, consensus seq is wt; else majority vote --> no wt veto in this verion

	2. aligned read has to start from pos 1 (or pos end) of fasta amplicon



Updated ver 1_3 (Oct 29 2019) By Sherry Chen, can handle indel now

	1. require fp perfect match and 10 nt after ER match: determine ER this way

Updated ver 1_3 (Oct 30 2019) By Sherry Chen, fix inconsistency with File_2 folder, outputs has “WTcorrection” as suffix

Updated ver 1_3 20201109 by Lucia Wu, require 9 nt after ER match

Written for Python 2.7.3

Required modules: Pysam, Samtools

Inputs:

	1. A position sorted SAM/BAM file containing reads with a UMI in the header.

		make sure its index file (.bai) is in the same path

		sort: 	samtools sort my.sam > my_sorted.bam

		index:	samtools index my_sorted.bam

	2. The fasta file used for alignment in the sam/bam file described in 1

	3. The enrichment region sequence (for wt) for each amplicon

	4. UMI length

	example to run in this folder:

	python UMI_counter.py --AlignmentFile lib1sorted.bam --FastaFile BDA_Mela.fasta --EnrichFile BDA_Mela_EnR.txt --UMIlen 15



Outputs:

	FileFolder1 (bam file name + _File1)

	a folder of files (each reference sequence has one file) with UMI and enrichment region count information 

	with name: outfile1_refID_EnrichmentSeq.txt

	with format:

		UMI:[UMI_sequence1]\n



		\t[read_enrichment_seq1]\t[de_novo call of read_enrichment_seq1]\t[count]\n



		\t[read_enrichment_seq2]\t[de_novo call of read_enrichment_seq2]\t[count]\n



		...



		UMI:[UMI_sequence2]\n



		\t[read_enrichment_seq1]\t[de_novo call of read_enrichment_seq1]\t[count]\n



		\t[read_enrichment_seq2]\t[de_novo call of read_enrichment_seq2]\t[count]\n



	



	FileFolder2 (bam file name + _File1)



	a folder of files (each reference sequence has one file) with UMI and majority_vote, family_size and frequency



	with name: outfile2_refID_EnrichmentSeq.txt	



	with format:



		[UMI_sequence]\t[most_pop_read_enrichment]\t[de_novo call of most_pop_read_enrichment]\t[family_size]\t[frequency]







	FinalSummaryFile (bam file name + _ResultSummary.txt)



	a file contains VAF AND VRF information for all loci



	with format



	VAF_refID\t seq1:var_call:count:frequency\tseq2:var_call:count:frequency\t.....\n



	VRF_refID\t seq1:var_call:count:frequency\tseq2:var_call:count:frequency\t.....\n



	...







usage: UMI_counter.py [-h] [--AlignmentFile] [--FastaFile] [--EnrichFile] [--UMIlen ]



Example: python UMI_counter.py --AlignmentFile lib2sorted.bam --FastaFile BDA_Mela.fasta --EnrichFile BDA_Mela_EnR.txt --UMIlen 15



arguments:



	-h, --help			show this help message and exit



	--AlignmentFile		input SAM/BAM file



	--FastaFile			input FASTA file



	--EnrichFile 		input enrichment region sequence file



	--UMIlen 			UMI length







implementation:



	Use pileupcolumn to grab all reads that maps over the region we want



	For each ref_id(pileupcolumn), create a dictionary to store UMI, and a Counter object to count same sequence



	ex: {UMI1:Counter(Seq1, 2; Seq2, 1; ...), UMI2:Counter(Seq3, 1; Seq2, 3; ...)}



	At the end of the pileupcolumn, the above dictionary is send to write_SNP_file function to save to file 



	in the output format described above 







some examples of calling de_novo mutations: (can be modified the de_novo function below)



	Enrichment Region	Var in ER (gBLock) 	Mutation calling



	GGTCAAGAGGAGTACAGTG	GGaCAAGAGGAGTACAGTG	3T>A



	GGTCAAGAGGAGTACAGTG	GGaaAAGAGGAGTACAGTG	3TC>AA



	GGTCAAGAGGAGTACAGTG	GGaatAGAGGAGTACAGTG	3TCA>AAT



	GGTCAAGAGGAGTACAGTG	GGCAAGAGGAGTACAGTG	3DelT



	GGTCAAGAGGAGTACAGTG	GGAGAGGAGTACAGTG	3DelTCA



	GGTCAAGAGGAGTACAGTG	GaaaGTCAAGAGGAGTACAGTG	2InsAAA





python UMI_counter3.py --AlignmentFile alignlib2.bam --FastaFile 14gene_tube1.fasta --EnrichFile 14gene_tube1_EnR.txt --UMIlen 15



"""

import pysam

from collections import Counter

import difflib 

from argparse import ArgumentParser

import os 
from math import ceil # added by Xiangjiang
#min_family = 3 # family size need to be larger than (not equal to) min_familty
min_family_static = 3

min_majority = 0



def main():

	# parser

	parser = ArgumentParser()

	parser.add_argument("--AlignmentFile", action="store", dest="alignment_file", help="input SAM/BAM file", required=True)

	parser.add_argument("--FastaFile",  action="store",  dest="fasta_file", help="Fasta file used for alignment", required=True)

	parser.add_argument("--EnrichFile",  action="store",  dest="enrich_file", help="Enrichment region sequence file", required=True)

	parser.add_argument("--UMIlen",  type=int, default=15,  dest="umi_len", help="length of UMI saved in header", required=True)

	o = parser.parse_args()

	caller(o.alignment_file, o.umi_len, o.fasta_file,o.enrich_file)



def generate_out_file_name(new_folder_name, ref_id, enrichment_seq):

	return new_folder_name + '/outfile_' + str(ref_id) + '_' + enrichment_seq + '.txt'



def find_loc_enrichment(amplicon_filename, enrichment_filename):

	## both of them in fasta format 

	enrichment_list = []

	amplicon_list = []

	amplicon_name_list = []

	enrichment_loc_list = []



	with open(amplicon_filename) as af:

		lines = af.readlines()

		for line in lines:

			if line[0] != '>' and len(line)> 5:

				amplicon_list.append(''.join(line.split()))

			elif line[0] == '>':

				amplicon_name_list.append(line[1:-1])



	with open(enrichment_filename) as ef:

		lines = ef.readlines()

		for line in lines:

			if line[0] != '>' and len(line)> 5:

				enrichment_list.append(''.join(line.split()))



	if len(amplicon_list) != len(enrichment_list):

		print ("ERROR: number of amplicons != number of enrichment region")

		return None



	for i in range(len(amplicon_list)):

		temp_amplicon = amplicon_list[i]

		temp_enrich = enrichment_list[i]

		loc1 = temp_amplicon.find(temp_enrich)



		if loc1 == -1:

			print(("ERROR: cannot find enrichment region in amplicon for input number" + str(i+1) ))

			print(("amplicon sequence: " + temp_amplicon))

			print(("enrichment sequence: " + temp_enrich))

			return None

		loc2 = loc1 + len(temp_enrich)

		enrichment_loc_list .append([loc1, loc2])

	return enrichment_loc_list, amplicon_name_list, amplicon_list, enrichment_list







def write_SNP_file(file_name, umi_list, wt):

	# create a dict for mut_call 

	mut_call_dict = {}

	with open(file_name, 'w') as f:

		for umi in umi_list:

			f.write('UMI: ' + umi + '\n')

			if sum(umi_list[umi].values()) > 0:

				for seq, seq_count in umi_list[umi].most_common():

					if seq not in mut_call_dict:

						mut_call = de_novo(seq, wt)

						mut_call_dict[seq] = mut_call

					else:

						mut_call = mut_call_dict[seq]

					#print "{}\t{}\t{}\n".format(seq, mut_call, seq_count)

					f.write("\t{}\t{}\t{}\n".format(seq, mut_call, seq_count))

	return mut_call_dict







def de_novo(seq, wt):

	if seq == wt:

		return "wt"

	# NO wildtype for now

	seq = seq + '0'

	wt = wt + '0'

	d = difflib.Differ()

	diff = list(d.compare(seq, wt))

	wt_loc_counter = 0 

	working_loc = 0

	working_str1 = ''

	working_str2 = ''

	all_mut_str = []

	for temp_d in diff:

		if temp_d[0] == ' ':

			if working_str1 != '' or working_str2 !='':

				if working_str1 == '':

					mut_str = str(working_loc) + 'Del' + working_str2

				elif working_str2 == '':

					mut_str = str(working_loc) + 'Ins' + working_str1

				else:

					mut_str = str(working_loc) + working_str2 + '>' + working_str1

				all_mut_str.append(mut_str)

			working_str1, working_str2 = '',''

			wt_loc_counter += 1

		if temp_d[0] == '-':

			wt_loc_counter += 1

			if working_str1 == '':

				working_loc = wt_loc_counter

			working_str1+=temp_d[2]

		if temp_d[0] == '+':

			if working_str1 == '' and working_str2 == '':

				working_loc = wt_loc_counter+1

			working_str2+=temp_d[2]

	return ','.join(all_mut_str)







def caller(samfile_name, umi_len, amplicon_filename, enrichment_filename):

	# find location of each enrichment region THEN find the name of each seuqnece for alignment

	enrichment_loc_list, fasta_name_list, fasta_seq_list, enrichment_seq_list = find_loc_enrichment(amplicon_filename, enrichment_filename)

	# create a new folder to save result 

	new_folder_name = samfile_name.replace('.','_')+'_File1'

	new_folder_name2 = samfile_name.replace('.','_')+'_File2'

	os.mkdir(new_folder_name)

	os.mkdir(new_folder_name2)

	all_VAF = []

	all_VRF = []

	final_outfile_name = samfile_name.replace('.','_')+'_ResultSummary.txt'

	# open samfile and use pileupcolumn to get all reads with valid enrichment region



	samfile = pysam.AlignmentFile(samfile_name, "rb" )



	for ref_id in range(len(fasta_name_list)):

		umi_list = {}

		VRF_collection  = Counter()

		VAF_collection = Counter()

		file2wirte = {}

		#require match full fp, and 9nts after ER, not the same as CNV code which requires full fp and full rp
		fp_seq = fasta_seq_list[ref_id][:enrichment_loc_list[ref_id][0]]

		#search_seq =  fasta_seq_list[ref_id][enrichment_loc_list[ref_id][1]:enrichment_loc_list[ref_id][1]+10]
		search_seq =  fasta_seq_list[ref_id][enrichment_loc_list[ref_id][1]:enrichment_loc_list[ref_id][1]+9]

		for pileupcolumn in samfile.pileup(fasta_name_list[ref_id], enrichment_loc_list[ref_id][0],enrichment_loc_list[ref_id][0]+1, max_depth = 1000000):

			if pileupcolumn.pos == enrichment_loc_list[ref_id][0]:

				for pileupread in pileupcolumn.pileups:

					read_umi = pileupread.alignment.query_name[-umi_len:]

					temp_loc0 = pileupread.alignment.query_sequence.find(fp_seq)

					temp_loc1 = pileupread.alignment.query_sequence.find(search_seq)

					if temp_loc0 == -1 or temp_loc1 == -1:

						continue



					read_nt = pileupread.alignment.query_sequence[temp_loc0+len(fp_seq):temp_loc1]



					# if pileupread.query_position != None:

					# 	read_nt = pileupread.alignment.query_sequence[pileupread.query_position:pileupread.query_position+enrichment_loc_list[ref_id][1]-enrichment_loc_list[ref_id][0]]

					# if pileupread.alignment.query_sequence[:10] == fasta_seq_list[ref_id][:10]: # from fp side

					# if pileupread.alignment.query_sequence[-10:] == fasta_seq_list[ref_id][-10:]: # from rp side

					VRF_collection[read_nt] += 1

					if read_umi not in umi_list:

						umi_list[read_umi] = Counter()

					umi_list[read_umi][read_nt] += 1



		outfile_name = new_folder_name + '/outfile1_' + str(ref_id) + '_' + enrichment_seq_list[ref_id]+ '.txt'

		mut_call_dict = write_SNP_file(outfile_name, umi_list, enrichment_seq_list[ref_id])	

		############ Added by Xiangjiang 04/24/20 #####################
		all_family_size = [sum(u.values()) for u in umi_list.values()]
		all_family_size.sort(reverse=True) # sort from large to small
		min_family_dynamic = ceil(0.05*sum(all_family_size[0:3])/3)
		min_family = max([min_family_static, min_family_dynamic])

		for curr_umi in umi_list:

			c = umi_list[curr_umi]

			most_pop_seq = c.most_common(1)[0][0] # the most common element

			vals = list(c.values())

			family_size = sum(vals)

			family_ratio = float(max(vals))/family_size







############ MOLECULE FILTER IS HERE, FAMILY SIZE, FAMILY RATIO, VALID UMI SEQUENCE ##############################







			if family_size > min_family and family_ratio > min_majority and 'G' not in curr_umi:

				mut_call = 'None'

				if most_pop_seq in mut_call_dict:

					mut_call = mut_call_dict[most_pop_seq]

					VAF_collection[most_pop_seq] += 1

				line2write = curr_umi + '\t' + most_pop_seq + '\t' + mut_call + '\t' + str(family_size) + '\t' + str(family_ratio) + '\n'

				if most_pop_seq not in file2wirte:

					file2wirte[most_pop_seq] = []

				file2wirte[most_pop_seq].append(line2write)

		outfile_name2 = new_folder_name2 + '/outfile2_' + str(ref_id) + '_' + enrichment_seq_list[ref_id]+ '.txt'

		with open(outfile_name2,'w') as fout2:

			for ele in file2wirte:

				fout2.writelines(file2wirte[ele])

		with open(final_outfile_name,'a') as ff:

			total_VAF_count = float(sum(VAF_collection.values()))

			total_VRF_count = float(sum(VRF_collection.values()))

			ff.write('VAF_' + str(ref_id)+'\t')

			for ele, count in VAF_collection.most_common():

				mut_call = "None"

				if ele in mut_call_dict:

					mut_call = mut_call_dict[ele]

				ff.write('' + ele +':' +mut_call+':' +str(count)+':'+ str(round(count/total_VAF_count,5))+'\t')

			ff.write('\n')

			ff.write('VRF_' + str(ref_id)+'\t')

			for ele, count in VRF_collection.most_common():

				mut_call = "None"

				if ele in mut_call_dict:

					mut_call = mut_call_dict[ele]

				ff.write(''+ ele +':' +mut_call+':' +str(count)+':'+ str(round(count/total_VRF_count,5))+'\t')

			ff.write('\n')



























# print find_loc_enrichment("BDA_Mela.fasta", "BDA_Mela_EnR.txt")



# caller("test.sorted.sam", 15, "BDA_Mela.fasta", "BDA_Mela_EnR.txt")











if __name__ == "__main__":



	main()



