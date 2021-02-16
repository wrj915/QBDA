% analyze "alignlibxx_bam_ResultSummary.txt" files
% output "libxx_varlist.xls" and "libxx_varlist_genomic.xls"(genomic coordinates)
% 
% close all
% clear all

%% File names and input amount
filename = 'alignlib1_bam_ResultSummary.txt';

input_ng = 1000;
fix_thresh = 6;
vaf_thresh = 0.003/100;


%% Run analysis
[TotalMol VarMolList VarMol WTMol TotalReads VarReadList VarReads WTReads] = LoadVarSummary(filename);
[VarMolList_filter VarMol_filter] = FilterVarSummaryv20210127(VarMolList, VarMol, input_ng, fix_thresh, vaf_thresh);
WriteFilterVar2File(VarMolList_filter, VarMol_filter, VarMolList, VarMol, input_ng, filename);

