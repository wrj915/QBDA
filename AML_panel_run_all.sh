# Run all analysis code for R1 130nt + R2 21nt paired-end fastq file

# Please change "aml1_S3_L001_R1_001.fastq" and "aml1_S3_L001_R2_001.fastq" to your R1 and R2 fastq file names
# Please change "$PATH:/Applications/MATLAB_R2019a.app/bin/" to your Matlab path

export PATH=$PATH:/Applications/MATLAB_R2019a.app/bin/

echo Trim...
python adapter_trim_f60_v2.py --FastqR1 aml1_S3_L001_R1_001.fastq --FastqR2 aml1_S3_L001_R2_001.fastq --OutputName trim_Lib_1.fastq

echo Alignment...
bowtie2 -x BDA_Leukemia -U trim_Lib_1.fastq -S alignlib1.sam

echo Sorting...
python sort_index_v2.py --InputSam alignlib1.sam --OutputBam alignlib1.bam

echo VariantCall...
python UMI_counter3_Vote_0_dyAmplicon_20201109.py --AlignmentFile alignlib1.bam --FastaFile BDA_Leukemia__v20201109.fasta --EnrichFile BDA_Leukemia_EnR_v20201109.txt --UMIlen 15

echo VariantFilter...
matlab -nodisplay -nosplash -nodesktop -r "run('run_filter_v2.m');exit;"