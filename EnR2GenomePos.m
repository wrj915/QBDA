function mutstrnew = EnR2GenomePos(mutstr, chrnum, EnRpos, dir)
% translate mutation string to genomic positions
% input examples:
% 13A>G
% 5DelCC
% 10InsC
% 19C>T,24G>A

mutlist = strsplit(mutstr,',');
mutnum = size(mutlist,2);
newmutlist = cell(mutnum,1);

for curmutnum = 1:mutnum
    inddel = strfind(mutlist{curmutnum},'Del');
    indins = strfind(mutlist{curmutnum},'Ins');
    indmut = strfind(mutlist{curmutnum},'>');
    
    if length(inddel) == 1
        mutpos = str2num(mutlist{curmutnum}(1:inddel-1));
        delseq = mutlist{curmutnum}(inddel+3:end);
        if dir == 1
            newmutlist{curmutnum} = sprintf('%dDel%s',mutpos+EnRpos-1,delseq);
        elseif dir == 2
            newmutlist{curmutnum} = sprintf('%dDel%s',EnRpos-mutpos-length(delseq)+2,comp_str(delseq));
        else
            error('Wrong dir number.\n')
        end
        
    elseif length(indins) == 1
        mutpos = str2num(mutlist{curmutnum}(1:indins-1));
        insseq = mutlist{curmutnum}(indins+3:end);
        if dir == 1
            newmutlist{curmutnum} = sprintf('%dIns%s',mutpos+EnRpos-1,insseq);
        elseif dir == 2
            newmutlist{curmutnum} = sprintf('%dIns%s',EnRpos-mutpos+2,comp_str(insseq));
        else
            error('Wrong dir number.\n')
        end
        
    elseif length(indmut) == 1
        for tempind = fliplr([1:indmut-2])
            mutpos = str2num(mutlist{curmutnum}(1:tempind));
            if ~isempty(mutpos)
                break
            end
        end
        wtseq = mutlist{curmutnum}(tempind+1:indmut-1);
        mutseq = mutlist{curmutnum}(indmut+1:end);
        
        if dir == 1
            newmutlist{curmutnum} = sprintf('%d%s>%s',mutpos+EnRpos-1,wtseq,mutseq);
        elseif dir == 2
            newmutlist{curmutnum} = sprintf('%d%s>%s',EnRpos-mutpos-length(wtseq)+2,comp_str(wtseq),comp_str(mutseq));
        else
            error('Wrong dir number.\n')
        end
        
    else
        error('Wrong mutation name formatting.\n')
    end
    
    if chrnum == 23
        mutstrnew = 'ChrX:';
    elseif chrnum == 24
        mutstrnew = 'ChrY:';
    else
        mutstrnew = sprintf('Chr%d:',chrnum);
    end
    
    mutstrnew = [mutstrnew newmutlist{1}];
    if mutnum > 1
        for curmutnum = 2:mutnum
            mutstrnew = [mutstrnew ',' newmutlist{curmutnum}];
        end
    end
end
