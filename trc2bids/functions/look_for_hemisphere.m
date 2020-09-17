function metadata = look_for_hemisphere(metadata,hemisphere_idx,annots)

annots_hemi = annots(hemisphere_idx,2);

% putting format in order ECoG,strip,depth
annotshemisplit = strsplit([annots_hemi{:}],{';','Hemisphere;'});
annotshemisplit = annotshemisplit(~cellfun(@isempty,annotshemisplit));
rightloc = find(cellfun('length',regexpi(lower(annotshemisplit),'right')) == 1);
leftloc = find(cellfun('length',regexpi(lower(annotshemisplit),'left')) == 1);

locs_all = sort([rightloc, leftloc,size(annotshemisplit,2)+1]);

% right
if ~isempty(rightloc)
    righthemi = cell(1,size(rightloc,2));
    for i=1:size(rightloc,2)
        righthemi{i} = annotshemisplit(rightloc(i)+1: locs_all(find(locs_all==rightloc(i))+1)-1);
    end
    righthemiall =  strcat([righthemi{:}],';');
    righthemifin = ['right;' righthemiall{:}];
else
    righthemifin = [];
end

% left
if ~isempty(leftloc)
    lefthemi = cell(1,size(leftloc,2));
    for i=1:size(leftloc,2)
        lefthemi{i} = annotshemisplit(leftloc(i)+1: locs_all(find(locs_all==leftloc(i))+1)-1);
    end
    lefthemiall = strcat([lefthemi{:}],';');
    lefthemifin = ['left;' lefthemiall{:}];
else
    lefthemifin = [];
end

metadata.hemisphere = [righthemifin, lefthemifin];
end
