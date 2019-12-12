%% create dataset descriptor
function create_eventDesc(proj_dir)

edesc_json.onset                                = 'onset of event in seconds' ;
edesc_json.offset                               = 'offset of event in seconds' ;
edesc_json.duration                             = 'duration of event in seconds' ;
edesc_json.trial_type                           = 'type of event (electrical stimulation/motor task/sensing task/artefact/sleep/sleep wake transition/eyes open)' ;
edesc_json.subtype                              = 'more description of event (sleep:nrem/rem, motor:Mario/hand/jump, sens:circle, electrical stimulation:SPES/ESM/REC2stim)' ;
edesc_json.electrodes_involved_onset            = 'electrodes involved in onset. For example: in seizure: electrodes involved in seizure onset or in artefact.' ;
edesc_json.electrodes_involved_offset           = 'electrodes involved in offset. For example: in seizure: electrodes involved in the end of a seizure or in an artefact.' ;
edesc_json.sample_start                         = 'onset of event in samples' ;
edesc_json.sample_end                           = 'offset of event in samples' ;
edesc_json.electrical_stimulation_type          = 'type of electrical stimulation [mono-/biphasic]';
edesc_json.electrical_stimulation_site          = 'electrode names of stimulus pair';
edesc_json.electrical_stimulation_site_num_1    = 'electrode one in stimulus pair' ;
edesc_json.electrical_stimulation_site_num_2    = 'electrode two in stimulus pair' ;
edesc_json.electrical_stimulation_current       = 'electrical stimulation current in Ampere';
edesc_json.frequency                            = 'electrical stimulation frequency in Hertz';
edesc_json.notes                                = 'notes about stimulation current';

if ~isempty(edesc_json)
    
    filename = fullfile(proj_dir,'event_description.json');
    if isfile(filename)
        existing = read_json(filename);
    else
        existing = [];
    end
    write_json(filename, mergeconfig(existing, edesc_json))
end
