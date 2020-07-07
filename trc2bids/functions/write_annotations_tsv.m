%% write annotations to a tsv file _annotations
function annotation_tsv = write_annotations_tsv(cfg,metadata,header,annots,fevents_name)


%% type / sample start / sample end /  chname;
% ch_label  = metadata.ch_label;
% ch2use_included = metadata.ch2use_included;
% fs = header.Rate_Min;

eventsannots = struct();
eventsannots.type        = {};
eventsannots.sub_type    = {};
eventsannots.s_start     = {};
eventsannots.s_end       = {};
eventsannots.ch_name_on  = {};
eventsannots.ch_name_off = {};

annots_new = annots;

%% artefacts
[annots_new, eventsannots ] = add_event2annotation(metadata.artefacts,'artefact',eventsannots, annots_new,header);

%% seizures
[annots_new, eventsannots ] = add_event2annotation(metadata.seizure,'seizure',eventsannots, annots_new,header);

%% stimulation
[annots_new, eventsannots ] = add_event2annotation(metadata.stimulation,'stimulation',eventsannots, annots_new,header);

%% eyes open
[annots_new, eventsannots ] = add_event2annotation(metadata.eyes_open,'eyes',eventsannots, annots_new,header);

%% sleep wake transition
[annots_new, eventsannots ] = add_event2annotation(metadata.slaw_trans,'sleep-wake transition',eventsannots, annots_new,header);

%% sleep
[annots_new, eventsannots ] = add_event2annotation(metadata.sleep,'sleep',eventsannots, annots_new,header);

%% motortask
[annots_new, eventsannots ] = add_event2annotation(metadata.motortask,'motortask',eventsannots, annots_new,header);

%% language task
[annots_new, eventsannots ] = add_event2annotation(metadata.langtask,'languagetask',eventsannots, annots_new,header);

%% sensing task
[annots_new, eventsannots ] = add_event2annotation(metadata.senstask,'sensingtask',eventsannots, annots_new,header);

%% visual selection data segments - slow wave sleep (or ar least late NREM stage) - aim: 10 consecutive minutes
[annots_new, eventsannots ] = add_event2annotation(metadata.SWSselection,'sws selection',eventsannots, annots_new,header);

%% visual selection data segments - rapid eye-movement sleep - aim: 10 consecutive minutes
[annots_new, eventsannots ] = add_event2annotation(metadata.REMselection,'rem selection',eventsannots, annots_new,header);

%% visual selection data segments - interictal awake - aim: 10 consecutive minutes
[annots_new, eventsannots ] = add_event2annotation(metadata.IIAWselection,'iiaw selection',eventsannots, annots_new,header);

%% visual selection data segments - Epileptogenicity index - about 5s before visual SOZ to 10s after, but sometimes different due to big artefacts/other signal characteristics
[annots_new, eventsannots ] = add_event2annotation(metadata.EIselection,'EI selection',eventsannots, annots_new,header);

%% adding trigger data to events list
% skip following section if no spes/esm/stimulation has been annotated in
% the file
if ~isempty(metadata.stimulation)
    eventsannots = add_stimulation2annotation(cfg,metadata,header,annots_new,eventsannots,fevents_name);
end

%% strip eventsannots into separate variables
if isempty(eventsannots.s_start)
    s_start         = 'n/a';
    s_end           = 'n/a';
    duration        = 'n/a';
    type            = 'n/a';
    sub_type        = 'n/a';
    ch_name_on      = 'n/a';
    ch_name_off     = 'n/a';
    samp_start      = 'n/a';
    samp_end        = 'n/a';
    stim_type       = 'n/a';
    site_name       = 'n/a';
    site_channum    = 'n/a';
    stim_cur        = 'n/a';
    notes           = 'n/a';
    freq            = 'n/a';
else
    s_start         = eventsannots.s_start;
    s_end           = eventsannots.s_end;
    duration        = eventsannots.duration;
    type            = eventsannots.type;
    sub_type        = eventsannots.sub_type;
    ch_name_on      = eventsannots.ch_name_on;
    ch_name_off     = eventsannots.ch_name_off;
    samp_start      = eventsannots.samp_start;
    samp_end        = eventsannots.samp_end;
    stim_type       = eventsannots.stim_type;
    site_name       = eventsannots.site_name;
    site_channum    = eventsannots.site_channum;
    stim_cur        = eventsannots.stim_cur;
    notes           = eventsannots.notes;
    freq            = eventsannots.freq;
    
end

% make table
annotation_tsv  = table(s_start,duration, type, sub_type, ch_name_on, ch_name_off,s_end, samp_start, samp_end, stim_type, site_name, stim_cur, freq, notes,  ...
    'VariableNames',{'onset', 'duration','trial_type', 'sub_type','electrodes_involved_onset','electrodes_involved_offset','offset','sample_start','sample_end','electrical_stimulation_type','electrical_stimulation_site','electrical_stimulation_current','electrical_stimulation_frequency','notes' });

% write table
if ~isempty(annotation_tsv)
    for i=1:size(cfg(1).ieeg_dir,2)
        
        filename = fullfile(cfg(1).ieeg_dir{i},fevents_name);
        
        if isfile(filename)
            existing = read_tsv(filename);
        else
            existing = [];
        end % try
        if ~isempty(existing)
            ft_error('existing file is not empty');
        end
        write_tsv(filename, annotation_tsv);
    end
end

