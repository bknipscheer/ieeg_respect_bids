%% create scan descriptor
function create_scansDesc(proj_dir)

scansdesc_json.filename                 = 'name of the file' ;
scansdesc_json.format                   = '"included" means that Format annotation is included in the trc file, "not included" means that Format annotation was not included in the trc file. In at least 1 trc file, Format annotation must be present.' ;
scansdesc_json.artefact                 = 'number of artefacts annotated' ;
scansdesc_json.sleep_total              = 'period (s) of sleep. Sleep is not defined as rem or non-rem.' ;
scansdesc_json.sleep_rem                = 'period (s) of sleep, specifically defined as rem-sleep.' ;
scansdesc_json.sleep_nrem               = 'period (s) of sleep, specifically defined as non-rem-sleep.' ;
scansdesc_json.sleepwaketransition      = 'period (s) where the patient is waking up or falling asleep and it was difficult to see whether the patient was either sleeping or resting with eyes closed.' ;
scansdesc_json.seizure_total            = 'total number of seizures' ;
scansdesc_json.seizure_clinical         = 'total number of seizures with clinical symptoms visible in the video' ;
scansdesc_json.seizure_subclinical      = 'total number of seizures with no clinical symptoms visible in the video';
scansdesc_json.motor                    = 'period (s) with motor tasks';
scansdesc_json.spes                     = 'total number of electrical stimuli with Single Pulse Electrical Stimulation (0.2 Hz, 1000 us, 4-8mA, monophasic)';
scansdesc_json.rec2stim                 = 'total number of electrical stimuli for the clinical trial REC2Stim (clinicaltrials.gov: NCT04158531)';
scansdesc_json.esm                      = 'total period (s) with electrical stimulation mapping to delineat eloquent cortex';
scansdesc_json.chocs                    = 'total period (s) with CHOCs protocol (1Hz, biphasic: 1000 us, 15s, 1-3 mA and 2000 us, 30s 1-5mA) to evoke seizures and clinical symptoms';
scansdesc_json.language                 = 'period (s) with language tasks';
scansdesc_json.sens                     = 'period (s) with sensory tasks';
scansdesc_json.sws_se                   = 'bla';
scansdesc_json.rem_sel                  = 'bla';
scansdesc_json.iiaw_sel                 = 'bla';
scansdesc_json.EI_sel                   = 'bla';

if ~isempty(scansdesc_json)
    
    filename = fullfile(proj_dir,'scans.json');
    if isfile(filename)
        existing = read_json(filename);
    else
        existing = [];
    end
    write_json(filename, mergeconfig(existing, scansdesc_json))
end
