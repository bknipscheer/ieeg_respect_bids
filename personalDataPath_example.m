function cfg = personalDataPath_example(varargin)

% function that contains local data path, is ignored in .gitignore

if ~isempty(varargin{1})
    if isstruct(varargin{1})
        if sum(contains(fieldnames(varargin{1}),'sub_labels'))
            if contains(varargin{1}.sub_labels,'studyA')
                % studyA
                foldername = input('Choose SystemPlus-folder: [name repositories where ieeg files are located]: ','s');
                if strcmp(foldername,'repositoryA')
                    cfg(1).proj_dirinput = '/location/to/ieeg-files/';
                elseif strcmp(foldername,'repositoryB')
                    cfg(1).proj_dirinput = '/location/to/ieeg-files';
                elseif strcmp(foldername,'repositoryC')
                    cfg(1).proj_dirinput = '/location/to/ieeg-files';
                else
                    error('Foldername is not recognized')
                end
                cfg(2).proj_dirinput = '/location/to/bidsfiles';
                cfg(1).proj_diroutput = '/location/to/bidsfiles';
                cfg(2).proj_diroutput = '/second/location/to/copy/bidsfiles'; % optional: this could remain empty
                
            elseif contains(varargin{1}.sub_labels,'studyB')
                % studyB
                cfg(1).proj_dirinput = '/location/to/ieeg-files';
                cfg(2).proj_dirinput = '/location/to/bidsfiles';
                cfg(1).proj_diroutput = '/location/to/bidsfiles';
            elseif contains(varargin{1}.sub_labels,'studyC')
                % studyC
                cfg(1).proj_dirinput = '/location/to/ieeg-files';
                cfg(2).proj_dirinput = '/location/to/bidsfiles';
                cfg(1).proj_diroutput = '/location/to/bidsfiles';
                
            end
            
            pat = [input('What is the PAT-folder in micromed database? [PAT_XXX] ','s'),'/'];
            cfg(1).pathname = fullfile(cfg(1).proj_dirinput,pat);

            
        end
    end
end

cfg(1).elec_input = '/location/to/electrodes.xlsx';

fieldtrip_folder  = '/location/to/fieldtrip/';
% copy the private folder in fieldtrip to somewhere else
fieldtrip_private = '/location/to/fieldtrip_private/';
jsonlab_folder    = '/location/to/jsonlab/';
addpath(fieldtrip_folder)
addpath(fieldtrip_private)
addpath(jsonlab_folder)
ft_defaults


end