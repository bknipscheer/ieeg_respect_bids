function write_coordsystemJSON(cfg)

%%%% This script writes coordsystem JSON file
%%%% Make sure data matlab to JSON library is added
%%%% This can be found here: https://github.com/gllmflndn/JSONio

%%%% Dora Hermes, Jaap van der Aar, Giulio Castegnaro, Dorien van Blooijs 2019

electrodes_json_name = fullfile(cfg.ieeg_directory,...
    [cfg.sub_labels{:} '_' cfg.ses_label '_coordsystem.json']); 

% This line is to create a variable for the name of the file intendedFor
filename_T1w = fullfile(cfg.sub_labels{:},cfg.ses_label, 'anat', [cfg.sub_labels{:} '_' cfg.ses_label '_proc-deface_T1w.nii']);

% assign information and methodology
loc_json.iEEGCoordinateSystem  = 'Other';
loc_json.iEEGCoordinateUnits  = 'mm';
loc_json.iEEGCoordinateSystemDescription = 'The origin of the coordinate system is between the ears and the axis are in the RAS direction. The scaling is with respect to the individuals anatomical scan and no scaling or deformation have been applied to the individuals anatomical scan';
loc_json.IntendedFor = filename_T1w; 
loc_json.iEEGCoordinateProcessingDescription = 'Surface projection Hermes or Branco';
loc_json.iEEGCoordinateProcessingReference = 'Hermes et al., 2010 JNeuroMeth , Branco et al., 2018 JNeuroMeth';

jsonSaveDir = fileparts(electrodes_json_name);

% ensure there is a /ieeg/ folder 
if ~isfolder(jsonSaveDir)
    fprintf('Warning: directory to save json file does not exist, create: %s \n',jsonSaveDir)
end

% json_options.indent = '    '; % this just makes the json file look prettier 
% 
% % write JSON file
% jsonwrite(electrodes_json_name,loc_json,json_options)

write_json(electrodes_json_name, loc_json)

end


function write_json(filename, json)
json = remove_empty(json);
ft_info('writing %s\n', filename);
if ft_hastoolbox('jsonlab', 3)
    savejson('', json, filename);
else
    str = jsonencode(json);
    fid = fopen(filename, 'w');
    fwrite(fid, str);
    fclose(fid);
end
end

function s = remove_empty(s)
fn = fieldnames(s);
fn = fn(structfun(@isempty, s));
s = removefields(s, fn);
end



