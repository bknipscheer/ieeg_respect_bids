%% create dataset descriptor
function create_datasetDesc(proj_dir)

ddesc_json.Name               = 'RESPect' ;
ddesc_json.BIDSVersion        = 'BEP010';
ddesc_json.License            = 'Not licenced yet';
ddesc_json.Authors            = {'van Blooijs D., Demuru M., Zweiphenning W.J.E.M, Leijten F.S.S., Zijlmans M.'};
ddesc_json.Acknowledgements   = 'Huiskamp G.J.M.';
ddesc_json.HowToAcknowledge   = 'possible paper to quote' ;
ddesc_json.Funding            = 'Epi-Sign Project and Epilepsiefonds #17-07' ;
ddesc_json.ReferencesAndLinks = {'articles and/or links'};
ddesc_json.DatasetDOI         = 'DOI of the dataset if online';


if ~isempty(ddesc_json)
    
    filename = fullfile(proj_dir,'dataset_description.json');
    if isfile(filename)
        existing = read_json(filename);
    else
        existing = [];
    end
    write_json(filename, mergeconfig(existing, ddesc_json))
    %     json_options.indent = ' ';
    %     jsonwrite(filename, mergeconfig(existing, ddesc_json), json_options)
end
