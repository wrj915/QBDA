function [TotalMol VarMolList VarMol WTMol TotalReads VarReadList VarReads WTReads] = LoadVarSummary(filename)
% load file "alignlibxx_bam_ResultSummary.txt"

% example: filename = 'alignlib11_bam_ResultSummary.txt';

Tarcount = 1;

fid = fopen(filename);
tline=fgets(fid);
VAFlines{Tarcount,1} = strsplit(tline);
while ischar(tline)
    tline = fgets(fid);
    if tline(1) == 'V' && tline(2) == 'A'
        Tarcount = Tarcount + 1;
        VAFlines{Tarcount,1} = strsplit(tline);
    end
end
fclose(fid);

WTMol = zeros(size(VAFlines,1),1);
VarMol = zeros(size(VAFlines,1),1);
TotalMol = zeros(size(VAFlines,1),1);
VarMolList = cell(size(VAFlines,1),1);

Tarcount = 1;

fid = fopen(filename);
tline=fgets(fid);
tline=fgets(fid);
VRFlines{Tarcount,1} = strsplit(tline);
while ischar(tline)
    tline = fgets(fid);
    if tline(1) == 'V' && tline(2) == 'R'
        Tarcount = Tarcount + 1;
        VRFlines{Tarcount,1} = strsplit(tline);
    end
    
    
end
fclose(fid);

WTReads = zeros(size(VRFlines,1),1);
VarReads = zeros(size(VRFlines,1),1);
TotalReads = zeros(size(VRFlines,1),1);
TopVarReads = zeros(size(VRFlines,1),1);
VarReadList = cell(size(VRFlines,1),1);
%%
for i = 1:size(VAFlines,1) % each target
    VarFlag = 0;
    for j = 2:size(VAFlines{i,1},2)-1 % genotypes
        Genoinfo = strsplit(VAFlines{i,1}{1,j},':');
        if strcmp(Genoinfo{1,2},'wt') % wt found
            WTMol(i,1) = str2num(Genoinfo{1,3});
        elseif str2num(Genoinfo{1,3})>0
            VarFlag = VarFlag +1;
            VarMol(i,VarFlag) = str2num(Genoinfo{1,3});
            VarMolList{i,VarFlag} = Genoinfo{1,2};
        end
    end
    
    TotalMol(i,1) = WTMol(i,1) + sum(VarMol(i,:));
    
end


for i = 1:size(VRFlines,1) % each target
    VarFlag = 0;
    for j = 2:size(VRFlines{i,1},2)-1 %first 2 genotypes, WT and Mut
        Genoinfo = strsplit(VRFlines{i,1}{1,j},':');
        if strcmp(Genoinfo{1,2},'wt') % wt found
            WTReads(i,1) = str2num(Genoinfo{1,3});
        
        elseif str2num(Genoinfo{1,3})>0
            VarFlag = VarFlag +1;
            VarReads(i,VarFlag) = str2num(Genoinfo{1,3});
            VarReadList{i,VarFlag} = Genoinfo{1,2};  
        end
 % Record the Reads for Variant called as TOP Var in Molecule layer       
        if strcmp(Genoinfo{1,2},VarMolList{i,1})  
            TopVarReads(i,1) = str2num(Genoinfo{1,3});
        end
        
    end
    
    TotalReads(i,1) = WTReads(i,1) + sum(VarReads(i,:));
    
end
