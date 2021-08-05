# QBDA
QBDA NGS data analysis for AML panel

Input: R1 130nt + R2 21nt paired-end fastq file; example: "aml1_S3_L001_R1_001.fastq" and "aml1_S3_L001_R2_001.fastq"

Output: Spreadsheets containing information of all called mutations: "lib1_varlist_genomic.xls" (coordinates in human genome) and "lib1_varlist.xls" (coordinates in enrichment regions)

Code:

AML_panel_run_all.sh: Main MacOS shell script running all the analysis code

adapter_trim_f60_v2.py: Python code Illumina adapter trimming

sort_index_v2.py: Python code for sorting SAM file

UMI_counter3_Vote_0_dyAmplicon_20201109.py: Python code for for de novo variant call

run_filter_v2.m: Matlab code for variant filtering

LoadVarSummary.m, FilterVarSummaryv20210127.m, WriteFilterVar2File.m, comp_str.m, EnR2GenomePos.m: Matlab functions

Other files:

BDA_Leukemia__v20201109.fasta: Amplicon sequences

BDA_Leukemia_EnR_v20201109.txt: Enrichment region sequences

BDA_Leukemia*.bt2: Index files built by Bowtie2

EnRpos_dir_AMLQBDA.csv: Coordinates and directions of enrichment regions in the human genome

homopolymer_in_EnR_AMLQBDA.csv: Homopolymer positions

Requirement:

The software has been tested on MacOS (10.13.3 or 10.15.7) with the following installed:
Python 2.7.16 or Python 2.7.18
Matlab R2019a or Matlab R2020b
Bowtie 2 version 2.3.5.1
pysam

Instructions for use:
In terminal, navigate to the current folder containing codes and the fastq file to analyze, then type in "sh AML_panel_run_all.sh"
