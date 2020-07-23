function eventsannots = add_stimulation2annotation(cfg,metadata,header,annots_new,eventsannots,fevents_name)
% this function finds the stimulation annotations and the specific stimulus
% parameters with each stimulus.

trigger = metadata.trigger;
stimtriggers = trigger.pos(trigger.val>1000);

%% determine period of stimulation

for stim = 1:size(metadata.stimulation,1)
    stimperiod = metadata.stimulation{stim}.pos;
    evname = metadata.stimulation{stim}.type;
    
    spesend = stimperiod(2);
    
    %% find annotations with stimulation parameters
    numannots = find([annots_new{:,1}]>stimperiod(1) & [annots_new{:,1}]<stimperiod(2));
    div_annots = findStimAnnots(metadata,annots_new,numannots);
    
    annot_stim = div_annots.annot_stim;
    annot_neg = div_annots.annot_neg;
    annot_bi_mono = div_annots.annot_bi_mono;
    annot_pulsewidth = div_annots.annot_pulsdur;
    annot_freq = div_annots.annot_freq;
    annot_curr = div_annots.annot_curr;
    annot_note = div_annots.annot_note;
    
    %% default stimulation parameters
    % stimulus current
    if contains(lower(metadata.format_info),'seeg')
        default_stimcur = 2/1000;
    elseif contains(lower(metadata.format_info),'ecog')
        default_stimcur = 8/1000;
    end
    
    % stimulus type
    if contains(lower(evname),'spes') || contains(lower(evname),'rec2stim')
        default_bi_mono = 'monophasic';
    elseif contains(lower(evname),'slowesm')
        default_bi_mono = 'monophasic';
    elseif contains(lower(evname),'esm') || contains(lower(evname),'chocs') || contains(lower(evname),'treni')
        default_bi_mono = 'biphasic';
    end
    
    % stimulus frequency
    if contains(lower(evname),'spes')
        default_pulsewidth = 1/1000;
        default_freq = 0.2;
    elseif contains(lower(evname),'esm')
        default_freq = 50;
        default_pulsewidth = 1050/1000;
    elseif contains(lower(evname),'chocs1')
        default_freq = 1;
        default_pulsewidth = 1.025/1000; %1025 usec
    elseif contains(lower(evname),'chocs2')
        default_freq = 1;
        default_pulsewidth = 2.025/1000; %2025 usec
    end
    
    % note description: notification when stimcurr is unknown
    if any(contains(fieldnames(metadata),'stimcurr'))
        if contains(lower(metadata.stimcurr),'unknown')
            default_note = sprintf('Stimulation intensity is suggested to be %1.3f A but may differ when applied in eloquent tissue',default_stimcur);
        else
            default_note = 'n/a';
        end
    else
        default_note = 'n/a';
    end
    
    %% load ECoG if no triggers
    if sum(cellfun(@(x) contains(x,{'No trigger'}),annots_new(:,2)))>0
        dataName = [cfg(1).ieeg_dir{1},'/', replace(fevents_name,'_events.tsv','_ieeg.vhdr')];
        data = ft_read_data(dataName,'dataformat','brainvision_eeg');
    end
    
    %% for each stimulus pair annotation, determine stimulus parameters
    for i=1:size(annot_stim,1)
        
        if i<size(annot_stim,1)
            stim_stop = annot_stim{i+1,1};
        else
            stim_stop = spesend;
        end
        
        stim_start = annot_stim{i,1};
        stim_num = annot_stim{i,4};
        
        %% find which annotations are after the current stimulation pair annotation, and before the next stimulation pair annotation
        
        idx_neg = ismember(annot_neg,stim_start:stim_stop-1);
        if ~isempty(annot_curr)
            idx_curr = ismember(annot_curr(:,1),stim_start:stim_stop-1);
        else
            idx_curr = [];
        end
        idx_bi_mono = ismember(annot_bi_mono,stim_start:stim_stop-1);
        idx_pulsewidth = ismember(annot_pulsewidth,stim_start:stim_stop-1);
        idx_freq = ismember(annot_freq,stim_start:stim_stop-1);
        if ~isempty(annot_note)
            idx_note = ismember(vertcat(annot_note{:,1}),stim_start:stim_stop-1);
        else
            idx_note = [];
        end
        
        % an extra trigger should be entered when stimulation is not one
        % pulse with a certain pulse duration, but a pulse train (like with
        % esm/treni)
        if any(trigger.val == 0)
            idx_trigger = trigger.val == 0 & ismember(trigger.pos,stim_start:stim_stop-1);
        else
            idx_trigger = [];
        end
        
        %% find onset samples when stimulation starts
        
        % if no trigger is found in the entire data dile
        if sum(cellfun(@(x) contains(x,{'No trigger'}),annots_new(:,2)))>0 && isempty(find(ismember(stimtriggers,stim_start:stim_stop),1))
            locs = stim_start + findtrigger(data, stim_num,header.Rate_Min, stim_start,stim_stop);
            
        elseif any(ismember(vertcat(annots_new{contains(annots_new(:,2),'No trigger'),1}),stim_start:stim_stop)) % if No triggers is part of periof of stimulus pair
            locs = stim_start + findtrigger(data, stim_num,header.Rate_Min, stim_start,stim_stop);
            
        elseif ~isempty(find(ismember(stimtriggers,stim_start:stim_stop),1)) % when triggers are present
            locs = trigger.pos(ismember(trigger.pos,stim_start:stim_stop) & trigger.val>1000);
            
        else
            
            error('No triggers, but also no comment "No trigger" in trc-file')
        end
        
        % define onset of stimulation
        stim_locs = num2cell(locs');
        s_start = locs/header.Rate_Min;
        stim_s_start = num2cell(s_start');
        
        stim_s_end = cell(size(stim_locs,1),1);
        [stim_s_end{:}] = deal('n/a');
        
        stim_samp_end = cell(size(stim_locs,1),1);
        [stim_samp_end{:}] = deal('n/a');
        
        %% stim name (uses neg-annotation)
        
        samp_neg = [];
        if ~isempty(annot_neg)
            samp_neg = annot_neg(idx_neg);
            
            idx_stim_neg = false(1,size(stim_locs,1));
            for n=1:size(samp_neg,1)
                if n<size(samp_neg,1)
                    idx_stim_neg(n,:) =  locs>samp_neg(n) & locs<samp_neg(n+1);
                else
                    idx_stim_neg(n,:) = locs>samp_neg(n);
                end
            end
            
            idx_stim_neg = boolean(sum(idx_stim_neg,1));
            
            stim_site_name = cell(size(stim_locs,1),1);
            
            if any(idx_stim_neg)
                [stim_site_name{idx_stim_neg}] = deal([annot_stim{i,3} '-' annot_stim{i,2}]);
                [stim_site_name{~idx_stim_neg}] = deal([annot_stim{i,2} '-' annot_stim{i,3}]);
            else
                [stim_site_name{:}] =deal([annot_stim{i,2} '-' annot_stim{i,3}]);
            end
            
        else
            stim_site_name = cell(size(stim_locs,1),1);
            [stim_site_name{:}] =deal([annot_stim{i,2} '-' annot_stim{i,3}]);
            
        end
        
        %% stimulus current
        
        samp_curr = [];
        if ~isempty(annot_curr)
            
            samp_curr = annot_curr(idx_curr,1);
            val_curr = annot_curr(idx_curr,2);
            
            stim_cur = cell(size(stim_locs,1),1);
            [stim_cur{:}] = deal(default_stimcur);
            
            for n=1:size(samp_curr,1)
                if n<size(samp_curr,1)
                    idx_stim_curr = locs>samp_curr(n)& locs<samp_curr(n+1) ;
                    [stim_cur{idx_stim_curr}] = deal(val_curr(n)/1000);
                else
                    idx_stim_curr = locs>samp_curr(n);
                    [stim_cur{idx_stim_curr}] = deal(val_curr(n)/1000);
                end
            end
            
        else
            
            idx_stim_curr = true(1,size(stim_locs,1));
            stim_cur = cell(size(stim_locs,1),1);
            [stim_cur{idx_stim_curr}] = deal(default_stimcur);
        end
        
        %% stimulus frequency
        
        samp_freq = [];
        if ~isempty(annot_freq)
            samp_freq = annot_freq(idx_freq,1);
            val_freq = annot_freq(idx_freq,2);
            
            stim_freq = cell(size(stim_locs,1),1);
            [stim_freq{:}] = deal(default_freq);
            
            for n=1:size(samp_freq,1)
                if n<size(samp_freq,1)
                    idx_stim_freq = locs>samp_freq(n)& locs<samp_freq(n+1) ;
                    [stim_freq{idx_stim_freq}] = deal(val_freq(n));
                else
                    idx_stim_freq = locs>samp_freq(n);
                    [stim_freq{idx_stim_freq}] = deal(val_freq(n));
                end
            end
            
        else
            idx_stim_freq = true(1,size(stim_locs,1));
            stim_freq = cell(size(stim_locs,1),1);
            [stim_freq{idx_stim_freq}] = deal(default_freq);
        end
        
        %% pulse width
        
        samp_pulsewidth = [];
        if ~isempty(annot_pulsewidth)
            
            samp_pulsewidth = annot_pulsewidth(idx_pulsewidth,1);
            val_pulsewidth = annot_pulsewidth(idx_pulsewidth,2);
            
            stim_pulsewidth = cell(size(stim_locs,1),1);
            [stim_pulsewidth{:}] = deal(default_pulsewidth);
            
            for n=1:size(samp_pulsewidth,1)
                if n<size(samp_pulsewidth,1)
                    idx_stim_pulsewidth = locs>samp_pulsewidth(n)& locs<samp_pulsewidth(n+1) ;
                    [stim_pulsewidth{idx_stim_pulsewidth}] = deal(val_pulsewidth(n));
                else
                    idx_stim_pulsewidth = locs>samp_pulsewidth(n);
                    [stim_pulsewidth{idx_stim_pulsewidth}] = deal(val_pulsewidth(n));
                end
            end
            
        else
            idx_stim_pulsewidth = true(1,size(stim_locs,1));
            
            stim_pulsewidth = cell(size(stim_locs,1),1);
            [stim_pulsewidth{idx_stim_pulsewidth}] = deal(default_pulsewidth);
        end
        
        %% duration of stimulation (if spes than, duration is pulsdur, when it its ESM, than 50Hz can have a pulse duration of 200 usec, but a total duration of stimulation that is about 5sec)
        
        if any(idx_trigger)
            
            stim_duration = (trigger.pos(idx_trigger)-stim_start)/header.Rate_Min; % in seconds
        else
            stim_duration = cell(size(stim_locs,1),1);
            [stim_duration{:}] = deal(stim_pulsewidth{1});
        end
        
        
        %% stimulus type (monophasic/biphasic)
        
        if ~isempty(annot_bi_mono)
            warning('Bi/monophasic annotation is used, check code in add_stimulation2annotation')
            
            samp_note = annot_bi_mono(idx_bi_mono,1);
            
            stim_bi_mono = cell(size(stim_locs,1),1);
            [stim_bi_mono{:}] = deal(default_bi_mono);
            
            for n=1:size(samp_note,1)
                if n<size(samp_note,1)
                    idx_stim_bi_mono = locs>samp_note(n)& locs<samp_note(n+1) ;
                    [stim_bi_mono{idx_stim_bi_mono}] = deal('Biphasic');
                else
                    idx_stim_bi_mono = locs>samp_note(n);
                    [stim_bi_mono{idx_stim_bi_mono}] = deal('Biphasic');
                end
            end
            
        else
            idx_stim_bi_mono = true(1,size(stim_locs,1));
            
            stim_bi_mono = cell(size(stim_locs,1),1);
            [stim_bi_mono{idx_stim_bi_mono}] = deal(default_bi_mono);
        end
        
        
        %% type
        stim_type = cell(size(stim_locs,1),1);
        [stim_type{:}] = deal('electrical_stimulation');
        
        stim_sub_type = cell(size(stim_locs,1),1);
        [stim_sub_type{:}] = deal(evname);
        
        %% notes
        
        if ~isempty(annot_note)
            
            samp_note = vertcat(annot_note{idx_note,1});
            val_note = annot_note(idx_note,2);
            
            % find where stim pair is stimulated for second time but with
            % different settings (neg, bi, lower amplitude etc)
            allsamps = unique(vertcat(samp_note,samp_curr, samp_freq,samp_neg,samp_pulsewidth));
            sampstart = [stim_start; allsamps(allsamps>locs(1))];
            
            stim_note = cell(size(stim_locs,1),1);
            [stim_note{:}] = deal(default_note);
            
            for n=1:size(samp_note,1)
                
                % find start of stimulation with specific note
                dif_start_note = samp_note(n)-sampstart;
                idx_sampstart = dif_start_note == min(dif_start_note(dif_start_note>0));
                
                idx_stim_note = locs> sampstart(idx_sampstart) & locs<samp_note(n);
                [stim_note{idx_stim_note}] = deal(val_note{n});
            end
            
        else
            idx_stim_note = true(1,size(stim_locs,1));
            stim_note = cell(size(stim_locs,1),1);
            [stim_note{idx_stim_note}] = deal(default_note);
        end
        
        %% involved channels
        stim_ch_name_on = cell(size(stim_locs,1),1);
        [stim_ch_name_on{:}] = deal('n/a');
        
        stim_ch_name_off = cell(size(stim_locs,1),1);
        [stim_ch_name_off{:}] = deal('n/a');
        
        %% place all info in
        
        stimannots(stim).type{i} = stim_type;
        stimannots(stim).sub_type{i} = stim_sub_type;
        stimannots(stim).stim_type{i} = stim_bi_mono;
        stimannots(stim).samp_start{i} = stim_locs;
        stimannots(stim).s_start{i} = stim_s_start; % time in seconds
        stimannots(stim).site_name{i} = stim_site_name;
        stimannots(stim).stim_pulsewidth{i} = stim_pulsewidth;
        stimannots(stim).duration{i} = stim_duration;
        stimannots(stim).samp_end{i} = stim_samp_end;
        stimannots(stim).s_end{i} = stim_s_end;
        stimannots(stim).stim_cur{i} = stim_cur;
        stimannots(stim).notes{i} = stim_note;
        stimannots(stim).freq{i} = stim_freq;
        stimannots(stim).ch_name_on{i} = stim_ch_name_on;
        stimannots(stim).ch_name_off{i} = stim_ch_name_off;
        
    end
end

type = horzcat(stimannots(:).type);
sub_type = horzcat(stimannots(:).sub_type);
stim_type = horzcat(stimannots(:).stim_type);
samp_start = horzcat(stimannots(:).samp_start);
s_start = horzcat(stimannots(:).s_start);
site_name = horzcat(stimannots(:).site_name);
pulsewidth = horzcat(stimannots(:).stim_pulsewidth);
duration = horzcat(stimannots(:).duration);
samp_end = horzcat(stimannots(:).samp_end);
s_end = horzcat(stimannots(:).s_end);
stim_cur = horzcat(stimannots(:).stim_cur);
notes = horzcat(stimannots(:).notes);
freq = horzcat(stimannots(:).freq);
ch_name_on = horzcat(stimannots(:).ch_name_on);
ch_name_off = horzcat(stimannots(:).ch_name_off);

eventsannots.type = vertcat(vertcat(eventsannots.type(:)), vertcat(type{:}));
eventsannots.sub_type = vertcat(vertcat(eventsannots.sub_type(:)), vertcat(sub_type{:}));
eventsannots.stim_type = vertcat(vertcat(eventsannots.stim_type(:)), vertcat(stim_type{:}));
eventsannots.samp_start = vertcat(vertcat(eventsannots.samp_start(:)),vertcat(samp_start{:}));
eventsannots.s_start = vertcat(vertcat(eventsannots.s_start(:)), vertcat(s_start{:}));
eventsannots.site_name = vertcat(vertcat(eventsannots.site_name(:)), vertcat(site_name{:}));
eventsannots.duration = vertcat(vertcat(eventsannots.duration(:)), vertcat(duration{:}));
eventsannots.samp_end = vertcat(vertcat(eventsannots.samp_end(:)), vertcat(samp_end{:}));
eventsannots.s_end = vertcat(vertcat(eventsannots.s_end(:)), vertcat(s_end{:}));
eventsannots.stim_cur = vertcat(vertcat(eventsannots.stim_cur(:)), vertcat(stim_cur{:}));
eventsannots.pulsewidth = vertcat(vertcat(eventsannots.pulsewidth(:)), vertcat(pulsewidth{:}));
eventsannots.notes = vertcat(vertcat(eventsannots.notes(:)), vertcat(notes{:}));
eventsannots.freq = vertcat(vertcat(eventsannots.freq(:)), vertcat(freq{:}));
eventsannots.ch_name_on = vertcat(vertcat(eventsannots.ch_name_on(:)), vertcat(ch_name_on{:}));
eventsannots.ch_name_off = vertcat(vertcat(eventsannots.ch_name_off(:)), vertcat(ch_name_off{:}));


end


