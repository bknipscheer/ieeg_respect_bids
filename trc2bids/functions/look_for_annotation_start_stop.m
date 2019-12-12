function [artefacts]=look_for_annotation_start_stop(annots,str_start,str_stop,ch,header)

fs = header.Rate_Min;

start_art=find(contains(annots(:,2),str_start));
end_art=find(contains(annots(:,2),str_stop));

% in a seizure, it is possible that the involved electrodes cannot be
% described in one Sz_on or one Sz_off, this part deals with this
if strcmp(str_start,'Sz_on')
    
    sz_cont = find(contains(annots(:,2),'Sz_cont'), 1);
    if ~isempty(sz_cont)
        if annots{sz_cont,1} < 10*fs % Sz_cont is annotated at beginning of file, sz_cont should be added in start_art
                start_art = sort([start_art; sz_cont]);
        else
            end_art = sort([end_art; sz_cont]);
        end
    end
        
    
    if any(diff(end_art) == 1) % if two Sz_offs are after each other
        diffSzoffs = (diff(end_art));
        sampSzOffs = diff([annots{end_art,1}]);
        
        for i=1:size(sampSzOffs,2)
            if diffSzoffs(i) == 1 && sampSzOffs(i) < 1*fs % SzOffs belong to the same seizure
                if strcmp(annots{end_art(i),2}(end),';')
                    annot1 = annots{end_art(i),2};
                else
                    annot1 = [annots{end_art(i),2},';'];
                end
                annotsplit2 = strsplit(annots{end_art(i+1),2},'Sz_off;');
                annots{end_art(i),2} = '';
                annots{end_art(i+1),2} = [annot1, annotsplit2{2}];
            elseif diffSzoffs(i) == 1 && sampSzOffs(i) > 1*fs % two Sz_offs do not belong to the same seizure
                error('starts and ends did not match')
            end
        end
        
        end_art=find(contains(annots(:,2),str_stop));
    end
    
    if any(diff(start_art) == 1) % if two Sz_ons are after each other
        diffSzons = (diff(start_art));
        sampSzOns = diff([annots{start_art,1}]);
        
        for i=1:size(sampSzOns,2)
            if diffSzons(i) == 1 && sampSzOns(i) < 1*fs % SzOffs belong to the same seizure
                if strcmp(annots{start_art(i),2}(end),';')
                    annot1 = annots{start_art(i),2};
                else
                    annot1 = [annots{start_art(i),2},';'];
                end
                annotsplit2 = strsplit(annots{start_art(i+1),2},'Sz_on;');
                annots{start_art(i),2} = '';
                annots{start_art(i+1),2} = [annot1, annotsplit2{2}];
            elseif diffSzons(i) == 1 && sampSzOns(i) > 1*fs % two Sz_offs do not belong to the same seizure
                error('starts and ends did not match')
            end
        end
        start_art=find(contains(annots(:,2),str_start));
    end
end

if (length(start_art)~=length(end_art)) && ~strcmp(str_start,'Eyes_open')
    error('starts and ends did not match')
end

if strcmp(str_start,'Art_on')
    %this code works perfectly when an event is closed in which the same
    %electrodes are involved in 'on' and 'off', like artefact. This is not the situation
    %in case of a seizure. Different electrodes can be involved in the end
    %compared to the beginning. In case of a motortask/sleep, a subtype is
    %mentioned in '.._on' but not in '.._off'
    artefacts=cell(size(start_art));
    
    for i=1:numel(start_art)
        art=struct;
        matched_end=find(contains(annots(:,2),replace(annots{start_art(i),2},str_start,str_stop)));
        if(isempty(matched_end))
            error('start and stop %s does not match',annots{start_art(i),2});
        end
        if(length(matched_end)>1)
            matched_end=matched_end((matched_end-start_art(i))>0);
            [val,idx_closest]=min(matched_end);
            matched_end=matched_end(idx_closest);%take the closest in time
        end
        ch_art_idx=parse_annotation(annots{start_art(i),2},ch);
        
        if isempty({ch{logical(ch_art_idx)}})
            art.ch_names_on={'all'};
            art.ch_names_off={'all'};
        else
            art.ch_names_on={ch{logical(ch_art_idx)}};
            art.ch_names_off={ch{logical(ch_art_idx)}};
        end
        art.pos=[(annots{start_art(i),1}) annots{matched_end,1}];
        art.type = 'n/a';
        if any(contains(annots(:,2),'NB artefacts annotated in avg'))
            art.notes = 'Please note, artefacts annotated in avg!';
        else
            art.notes = 'n/a';
        end
        artefacts{i}=art;
    end
elseif  strcmp(str_start,'Sz_on') % in case of a seizure
    artefacts=cell(size(start_art));
    for i=1:numel(start_art)
        art=struct;
        
        % seizure onset
        annotsplit = strsplit(annots{start_art(i),2},';');
        annotsplit = annotsplit(~cellfun(@isempty,annotsplit));
        if size(annotsplit,2) >2 % both subtype of seizure and channelnames are mentioned 
            if strcmp(annotsplit{2},'clin') || strcmp(annotsplit{2},'cluster') || strcmp(annotsplit{2},'subclin') || strcmp(annotsplit{2},'aura')
                type = annotsplit{2};
            else
                type = 'unknown';
            end
            ch_art_idx = parse_annotation(annots{start_art(i),2},ch);
            ch_names_on = {ch{logical(ch_art_idx)}};
            
        elseif size(annotsplit,2) == 2 % either subtype of seizure or channelnames are mentiond
            if strcmp(annotsplit{2},'clin') || strcmp(annotsplit{2},'cluster') || strcmp(annotsplit{2},'subclin') || strcmp(annotsplit{2},'aura')
                type = annotsplit{2};
                
                if strcmp(annotsplit{1},'Sz_cont') % Sz_cont is annotated, so no channelnames
                    ch_names_on = {'continuation of seizure'};
                else
                    ch_names_on = {'diffuse'};
                end
            else
                type = 'unknown';
                ch_art_idx = parse_annotation(annots{start_art(i),2},ch);
                ch_names_on = {ch{logical(ch_art_idx)}};
            end
        elseif  size(annotsplit,2) ==1 && strcmp(annotsplit{1},'Sz_cont') % for example with Sz_cont;, it is possible that only Sz_cont is mentioned
            type = 'unknown';
            ch_names_on = {'continuation of seizure'};
        end
        
        % seizure offset
        annotsplit = strsplit(annots{end_art(i),2},';');
        annotsplit = annotsplit(~cellfun(@isempty,annotsplit));
        if size(annotsplit,2) > 1 % for example "Sz_off;C[13,14]
            ch_art_idx = parse_annotation(annots{end_art(i),2},ch);
            ch_names_off = {ch{logical(ch_art_idx)}};
        else
            if strcmp(annotsplit{1},'Sz_cont')
                ch_names_off = {'Seizure continues in next file'};
            else
                ch_names_off = {'diffuse'};
            end
        end
        
        art.type = type;
        art.ch_names_on = ch_names_on;
        art.ch_names_off = ch_names_off;
        art.pos = [(annots{start_art(i),1}) (annots{end_art(i),1})]; % each seizure is consecutive, so as long as a seizure's end has the correct annotation, this will work
        artefacts{i} = art;
    end
elseif  strcmp(str_start,'Eyes_open') % in case of eyes open
    all_art = sort([start_art;end_art]);
    artefacts=cell(size(all_art));
    
    for i=1:numel(all_art)
        art=struct;
       
        annotsplit = strsplit(annots{all_art(i),2},';');
        annotsplit = annotsplit(~cellfun(@isempty,annotsplit));
        
        if strcmp(annotsplit{1},'Eyes_open')
            type = 'opened';
        else % if annotsplit is Eyes_close
            type = 'closed';
        end        
        art.type = type;
        
        art.ch_names_on={'all'};
        art.ch_names_off={'all'};
        
        if i<numel(all_art)
            art.pos=[(annots{all_art(i),1}) annots{all_art(i+1),1}];
        elseif i==numel(all_art) % if last annotation is selected, this event takes till the end of the file
            art.pos=[(annots{all_art(i),1}) header.Num_Samples];
        end
        artefacts{i}=art;
    end
    
else 
    artefacts=cell(size(start_art));
    
    for i=1:numel(start_art)
        art=struct;
        
        %matched_end=find(contains(annots(:,2),replace(annots{start_art(i),2},str_start,str_stop)));
        %         if(isempty(matched_end))
        %             error('start and stop %s does not match',annots{start_art(i),2});
        %         end
        %         if(length(matched_end)>1)
        %             matched_end=matched_end((matched_end-start_art(i))>0);
        %             [val,idx_closest]=min(matched_end);
        %             matched_end=matched_end(idx_closest);%take the closest in time
        %         end
        
        annotsplit = strsplit(annots{start_art(i),2},';');
        annotsplit = annotsplit(~cellfun(@isempty,annotsplit));
        
        if size(annotsplit,2) == 1
            if strcmp(str_start,'Eyes_open')
                type = 'n/a';
            else % if str_start is motor/sens/language
                type = 'unknown';
            end
        elseif size(annotsplit,2) == 2
            type = annotsplit{2};
        end
        
        ch_art_idx=parse_annotation(annots{start_art(i),2},ch);
        
        art.type = type;
        if isempty({ch{logical(ch_art_idx)}})
            art.ch_names_on={'all'};
            art.ch_names_off={'all'};
        else
            art.ch_names_on={ch{logical(ch_art_idx)}};
            art.ch_names_off={ch{logical(ch_art_idx)}};
        end
        art.pos=[(annots{start_art(i),1}) annots{end_art(i),1}];
        artefacts{i}=art;
    end
end
