%% create dataset descriptor
function create_participantsDesc(proj_dir)

partdesc_json.participant_id    = 'name of the patient' ;
partdesc_json.session           = 'session number, session 1 is the first period of invasive recordings, 1a/b means that the electrodes have changes during the first session' ;
partdesc_json.age               = 'age of the patient during the session (years)' ;
partdesc_json.sex               = 'gender of the patient (male/female/unknown)' ;

if ~isempty(partdesc_json)
    
    filename = fullfile(proj_dir,'participants.json');
    if isfile(filename)
        existing = read_json(filename);
    else
        existing = [];
    end
    write_json(filename, mergeconfig(existing, partdesc_json))
end
