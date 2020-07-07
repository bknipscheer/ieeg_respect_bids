function eventsannots = add_stimulation2annotation(cfg,metadata,header,annots_new,eventsannots,fevents_name)
% this function finds the stimulation annotations and the specific stimulus
% parameters with each stimulus. 

trigger = metadata.trigger;
stimtriggers = trigger.pos(trigger.pos>1000);

% there are stimtriggers and the annotation No trigger
if ~isempty(stimtriggers) && sum(cellfun(@(x) contains(x,{'No trigger'}),annots_new(:,2)))>0
   warning('Check add_stimulation2annotation.m to check what happens in this situation') 
end

%% determine period of stimulation
% FIX THIS: this code now assumes that only 1 type of stimulation is applied in
% one trc-file.
stimperiod = metadata.stimulation{1}.pos;
evname = metadata.stimulation{1}.type;

spesend = stimperiod(2);

%% find annotations with stimulation parameters
numannots = find([annots_new{:,1}]>stimperiod(1) & [annots_new{:,1}]<stimperiod(2));
div_annots = findStimAnnots(metadata,annots_new,numannots);

annot_stim = div_annots.annot_stim;
annot_neg = div_annots.annot_neg;
annot_bi_mono = div_annots.annot_bi_mono;
annot_pulsdur = div_annots.annot_pulsdur;
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
elseif contains(lower(evname),'esm') || contains(lower(evname),'chocs') || contains(lower(evname),'treni')
    default_bi_mono = 'biphasic';
end

% stimulus frequency
if contains(lower(evname),'spes')
    default_pulsdur = 1/1000;
    default_freq = 0.2;
end

% note description: notification when stimcurr is unknown
if contains(fieldnames(metadata),'stimcurr')
    if contains(lower(metadata.stimcurr),'unknown')
        default_note = sprintf('Stimulation intensity is suggested to be %i mA but may differ when applied in eloquent tissue',default_stimcur);
    else
        default_note = 'n/a';
    end
else
    default_note = 'n/a';
end

%% load ECoG if no triggers
if sum(cellfun(@(x) contains(x,{'No trigger'}),annots_new(:,2)))>0
    dataName = [cfg(1).ieeg_dir{1},'/', replace(fevents_name,'_events.tsv','_ieeg.vhdr')];
    data_raw = ft_read_data(dataName,'dataformat','brainvision_eeg');
    data = data_raw(metadata.ch2use_included,:);
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
    idx_curr = ismember(annot_curr(:,1),stim_start:stim_stop-1);
    idx_bi_mono = ismember(annot_bi_mono,stim_start:stim_stop-1);
    idx_pulsdur = ismember(annot_pulsdur,stim_start:stim_stop-1);
    idx_freq = ismember(annot_freq,stim_start:stim_stop-1);
    if ~isempty(annot_note)
        idx_note = ismember(vertcat(annot_note{:,1}),stim_start:stim_stop-1);
    else
        idx_note = [];
    end
    
    %% find onset samples when stimulation starts
    
    % if no trigger is found
    if sum(cellfun(@(x) contains(x,{'No trigger'}),annots_new(:,2)))>0
        locs = stim_start + findtrigger(data, stim_num,header.Rate_Min, stim_start,stim_stop);
        
    elseif ~isempty(find(ismember(stimtrigger,stim_start:stim_stop),1)) % when triggers are present
        locs = trigger.pos(ismember(stimtrigger,stim_start:stim_stop));
        
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
        
        idx_stim_neg = true(size(stim_locs,1),1);
        for n=1:size(samp_neg,1)
            if n<size(samp_neg,1)
                idx_stim_neg = locs>samp_neg(n)&& locs<samp_neg(n+1);
            else
                idx_stim_neg = locs>samp_neg(n);
            end
        end
        
        stim_site_name = cell(size(stim_locs,1),1);
        
        if any(idx_stim_neg)
            [stim_site_name{idx_stim_neg}] = deal([annot_stim{i,3} '-' annot_stim{i,2}]);
            [stim_site_name{~idx_stim_neg}] = deal([annot_stim{i,2} '-' annot_stim{i,3}]);
        else
            [stim_site_name{:}] =deal([annot_stim{i,2} '-' annot_stim{i,3}]);
        end
        
    else
        
        warning('No neg annotations, so check script add_stimulation2annotation.m')
        
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
        warning('Freq is used, check code in add_stimulation2annotation')
        
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
    
    %% pulse duration
    
    samp_pulsdur = [];
    if ~isempty(annot_pulsdur)
        warning('Pulsdur is used, check code in add_stimulation2annotation')
        
        samp_pulsdur = annot_pulsdur(idx_pulsdur,1);
        val_pulsdur = annot_pulsdur(idx_pulsdur,2);
        
        stim_pulsdur = cell(size(stim_locs,1),1);
        [stim_pulsdur{:}] = deal(stimpulsdurdefault);
        
        for n=1:size(samp_pulsdur,1)
            if n<size(samp_pulsdur,1)
                idx_stim_pulsdur = locs>samp_pulsdur(n)& locs<samp_pulsdur(n+1) ;
                [stim_pulsdur{idx_stim_pulsdur}] = deal(val_pulsdur(n));
            else
                idx_stim_pulsdur = locs>samp_pulsdur(n);
                [stim_pulsdur{idx_stim_pulsdur}] = deal(val_pulsdur(n));
            end
        end
        
    else
        idx_stim_pulsdur = true(1,size(stim_locs,1));
        
        stim_pulsdur = cell(size(stim_locs,1),1);
        [stim_pulsdur{idx_stim_pulsdur}] = deal(default_pulsdur);
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
        allsamps = unique(vertcat(samp_note,samp_curr, samp_freq,samp_neg,samp_pulsdur));
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
    
    stimannots.type{i} = stim_type;
    stimannots.sub_type{i} = stim_sub_type;
    stimannots.stim_type{i} = stim_bi_mono;
    stimannots.samp_start{i} = stim_locs;
    stimannots.s_start{i} = stim_s_start; % time in seconds (1 decimal)
    stimannots.site_name{i} = stim_site_name;
    stimannots.site_channum{i} = stim_site_channum;
    stimannots.duration{i} = stim_pulsdur;
    stimannots.samp_end{i} = stim_samp_end;
    stimannots.s_end{i} = stim_s_end;
    stimannots.stim_cur{i} = stim_cur;
    stimannots.notes{i} = stim_note;
    stimannots.freq{i} = stim_freq;
    stimannots.ch_name_on{i} = stim_ch_name_on;
    stimannots.ch_name_off{i} = stim_ch_name_off;
    
end

eventsannots.type = vertcat(vertcat(eventsannots.type(:)), vertcat(stimannots.type{:}));
eventsannots.sub_type = vertcat(vertcat(eventsannots.sub_type(:)), vertcat(stimannots.sub_type{:}));
eventsannots.stim_type = vertcat(vertcat(eventsannots.stim_type(:)), vertcat(stimannots.stim_type{:}));
eventsannots.samp_start = vertcat(vertcat(eventsannots.samp_start(:)),vertcat(stimannots.samp_start{:}));
eventsannots.s_start = vertcat(vertcat(eventsannots.s_start(:)), vertcat(stimannots.s_start{:}));
eventsannots.site_name = vertcat(vertcat(eventsannots.site_name(:)), vertcat(stimannots.site_name{:}));
eventsannots.site_channum = vertcat(vertcat(eventsannots.site_channum(:)), vertcat(stimannots.site_channum{:}));
eventsannots.duration = vertcat(vertcat(eventsannots.duration(:)), vertcat(stimannots.duration{:}));
eventsannots.samp_end = vertcat(vertcat(eventsannots.samp_end(:)), vertcat(stimannots.samp_end{:}));
eventsannots.s_end = vertcat(vertcat(eventsannots.s_end(:)), vertcat(stimannots.s_end{:}));
eventsannots.stim_cur = vertcat(vertcat(eventsannots.stim_cur(:)), vertcat(stimannots.stim_cur{:}));
eventsannots.notes = vertcat(vertcat(eventsannots.notes(:)), vertcat(stimannots.notes{:}));
eventsannots.freq = vertcat(vertcat(eventsannots.freq(:)), vertcat(stimannots.freq{:}));
eventsannots.ch_name_on = vertcat(vertcat(eventsannots.ch_name_on(:)), vertcat(stimannots.ch_name_on{:}));
eventsannots.ch_name_off = vertcat(vertcat(eventsannots.ch_name_off(:)), vertcat(stimannots.ch_name_off{:}));


end


