function metadata = look_for_format(metadata,format_idx,annots)

annots_format = annots(format_idx,2);

% putting format in order ECoG,strip,depth
annotsformatsplit = strsplit([annots_format{:}],{';','Format;'});
annotsformatsplit = annotsformatsplit(~cellfun(@isempty,annotsformatsplit));
ecogloc = find(cellfun('length',regexp(lower(annotsformatsplit),'ecog')) == 1);
striploc = find(cellfun('length',regexp(lower(annotsformatsplit),'strip')) == 1);
depthloc = find(cellfun('length',regexp(lower(annotsformatsplit),'depth')) == 1);
seegloc = find(cellfun('length',regexp(lower(annotsformatsplit),'seeg')) == 1);

locs_all = sort([ecogloc, striploc, depthloc,seegloc,size(annotsformatsplit,2)+1]);

% ECoG
if ~isempty(ecogloc)
    ecogformat = cell(1,size(ecogloc,2));
    for i=1:size(ecogloc,2)
        ecogformat{i} = annotsformatsplit(ecogloc(i)+1: locs_all(find(locs_all==ecogloc(i))+1)-1);
    end
    ecogformatall =  strcat([ecogformat{:}],';');
    ecogformatfin = ['ECoG;' ecogformatall{:}];
else
    ecogformatfin = [];
end

% strip
if ~isempty(striploc)
    stripformat = cell(1,size(striploc,2));
    for i=1:size(striploc,2)
        stripformat{i} = annotsformatsplit(striploc(i)+1: locs_all(find(locs_all==striploc(i))+1)-1);
    end
    stripformatall = strcat([stripformat{:}],';');
    stripformatfin = ['strip;' stripformatall{:}];
else
    stripformatfin = [];
end

% depth
if ~isempty(depthloc)
    depthformat = cell(1,size(depthloc,2));
    for i=1:size(depthloc,2)
        depthformat{i} = annotsformatsplit(depthloc(i)+1: locs_all(find(locs_all==depthloc(i))+1)-1);
    end
    depthformatall = strcat([depthformat{:}],';');
    depthformatfin = ['depth;' depthformatall{:}];
else
    depthformatfin = [];
end

% seeg
if ~isempty(seegloc)
    seegformat = cell(1,size(seegloc,2));
    for i=1:size(seegloc,2)
        seegformat{i} = annotsformatsplit(seegloc(i)+1: locs_all(find(locs_all==seegloc(i))+1)-1);
    end
    seegformatall = strcat([seegformat{:}],';');
    seegformatfin = ['seeg;' seegformatall{:}];
else
    seegformatfin = [];
end

metadata.format_info=[ecogformatfin, stripformatfin, depthformatfin, seegformatfin];
end
