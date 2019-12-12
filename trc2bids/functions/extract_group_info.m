% extract group information
% assumption the included are only grid and strip
function ch_group = extract_group_info(metadata)

ch_label                                    = metadata.ch_label                    ;
ch2use_included = metadata.ch2use_included;

if contains(lower(metadata.format_info),'seeg')
    idx_depths = ch2use_included;
    idx_strips = false(size(ch_label));
    idx_grid = false(size(ch_label));
elseif contains(lower(metadata.format_info),'ecog')
    C = strsplit(metadata.format_info,{';','['});
    
    id_ngrid = regexpi(C,'1');
    ngridnum = find(cellfun(@isempty,id_ngrid)==0);
    
    id_strip = regexpi(C,'strip');
    stripnum = find(cellfun(@isempty,id_strip)==0);
    id_depth = regexpi(C,'depth');
    depthnum = find(cellfun(@isempty,id_depth)==0);
    
    idx_depths = zeros(size(ch_label));
    idx_strips = zeros(size(ch_label));
    for i=1:size(ngridnum,2)
        if ~isempty(stripnum) && isempty(depthnum)
            if any(ngridnum(i) > stripnum)
                idx_strip = regexpi(ch_label,C{ngridnum(i)-1});
                idx_strip = cellfun(@isempty,idx_strip);
                idx_strip = ~idx_strip;
                idx_strips = idx_strips + idx_strip;
            end
        elseif isempty(stripnum) && ~isempty(depthnum)
            if any(ngridnum(i) > depthnum)
                idx_depth = regexpi(ch_label,C{ngridnum(i)-1});
                idx_depth = cellfun(@isempty,idx_depth);
                idx_depth = ~idx_depth;
                idx_depths = idx_depths + idx_depth;
            end
        elseif ~isempty(stripnum) && ~isempty(depthnum)
            if depthnum < stripnum
                if any(ngridnum(i) > depthnum) && any(ngridnum(i) < stripnum)
                    idx_depth = regexpi(ch_label,C{ngridnum(i)-1});
                    idx_depth = cellfun(@isempty,idx_depth);
                    idx_depth = ~idx_depth;
                    idx_depths = idx_depths + idx_depth;
                elseif ngridnum(i) > stripnum
                    idx_strip = regexpi(ch_label,C{ngridnum(i)-1});
                    idx_strip = cellfun(@isempty,idx_strip);
                    idx_strip = ~idx_strip;
                    idx_strips = idx_strips + idx_strip;
                end
            elseif depthnum > stripnum
                if ngridnum(i) > stripnum && ngridnum(i) < depthnum
                    idx_strip = regexpi(ch_label,C{ngridnum(i)-1});
                    idx_strip = cellfun(@isempty,idx_strip);
                    idx_strip = ~idx_strip;
                    idx_strips = idx_strips + idx_strip;
                    
                elseif ngridnum(i) > depthnum
                    % this finds all electrodes with a specific name, so
                    % for example when depth electrode is called 'L', then
                    % the grid electrodes 'CL'  will also be noted as depth
                    % electrodes
                    idx_depth = regexpi(ch_label,C{ngridnum(i)-1});
                    % this searches for when the 'L' is mentioned as first,
                    % otherwise the 'L' is part of another electrodename
                    idx_depth = cellfun(@(x) find(x==1), idx_depth, 'UniformOutput', false);
                    idx_depth = cellfun(@isempty,idx_depth);
                    idx_depth = ~idx_depth;
                    idx_depths = idx_depths + idx_depth;
                end
            end
        end
    end
    idx_grid = ~idx_depths & ~idx_strips & ch2use_included;
end

idx_depths = logical(idx_depths);
idx_strips = logical(idx_strips);

ch_group                                    = cell(size(ch2use_included)) ;
if(any(idx_grid))
    [ch_group{idx_grid}]                    = deal('grid')                         ;
end
if(any(idx_strips))
    [ch_group{idx_strips}]                   = deal('strip')                        ;
end
if(any(idx_depths))
    [ch_group{idx_depths}]                   = deal('depth')                        ;
end
if(any(~ch2use_included))
    [ch_group{ ~ch2use_included }] = deal('other')                        ;
end
