function [annots_new, eventsannots ] = add_event2annotation(event,evname,eventsannots, annots_new, header)

fs = header.Rate_Min;

if size(eventsannots.type,2) == 0
    cc = 1;
else
    cc = size(eventsannots.type,2)+1;
end

if(~isempty(event))
    for i=1:numel(event)
        
        eventsannots.type{cc}    = evname                           ;
        
        %         if strcmp(evname,'artefact')
        eventsannots.sub_type{cc} = event{i}.type;
        eventsannots.s_start{cc} = round(event{i}.pos(1)/fs,1); % time in seconds (1 decimal)
        eventsannots.samp_start{cc} = num2str(event{i}.pos(1))          ;
        eventsannots.s_end{cc}   = round(event{i}.pos(end)/fs,1); % time in seconds (1 decimal)
        eventsannots.samp_end{cc} = num2str(event{i}.pos(end))          ;
        eventsannots.duration{cc} = round(eventsannots.s_end{cc} - eventsannots.s_start{cc},3);
        eventsannots.stim_type{cc} = 'n/a';
        eventsannots.site_name{cc} = 'n/a';
        eventsannots.site_channum{cc} = 'n/a';
        eventsannots.stim_cur{cc} = 'n/a';
        eventsannots.freq{cc} = 'n/a';
       
        if isfield(event{i},'notes') && ~isempty(event{i}.notes)
            eventsannots.notes{cc} = event{i}.notes;
        else
            eventsannots.notes{cc} = 'n/a';
        end
        
        if size(event{i}.ch_names_on,2) == 1
            name = event{i}.ch_names_on{1}              ;
        else
            for j=1
                name =[event{i}.ch_names_on{1}];
            end
            
            for j=2:size(event{i}.ch_names_on,2)
                name = [name ,',', event{i}.ch_names_on{j}];
            end
        end
        eventsannots.ch_name_on{cc} = name;
        
        if size(event{i}.ch_names_off,2) == 1
            name = event{i}.ch_names_off{1}              ;
        else
            for j=1
                name =[event{i}.ch_names_off{1}];
            end
            
            for j=2:size(event{i}.ch_names_off,2)
                name = [name ,',', event{i}.ch_names_off{j}];
            end
        end
        
        eventsannots.ch_name_off{cc} = name;
        
        annots_new([annots_new{:,1}]==event{i}.pos(1),:)=[];
        annots_new([annots_new{:,1}]==event{i}.pos(end),:)=[];
        
        cc          = cc + 1                               ;
        
        %         end
    end
end
