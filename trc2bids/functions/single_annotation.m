function [ch_parsed]=single_annotation(annots,keyWord,ch)


ch_idx=cellfun(@(x) contains(x,{keyWord}),annots(:,2));

if(sum(ch_idx)<1)
    error('Missing annotation (example "%s;Gr01;Gr[3:5]")',keyWord)
end
ch_parsed=zeros(size(ch));
if(sum(ch_idx))
    str2parse={annots{ch_idx,2}};
    for i=1:numel(str2parse)
        C=strsplit(str2parse{i},';');
        C=C(~cellfun(@isempty,C));
        if(numel(C)>1)%TODO better check
            ch_parsed= ch_parsed | parse_annotation(str2parse{i},ch);
        end
    end
end
