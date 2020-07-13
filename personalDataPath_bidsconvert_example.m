function cfg = personalDataPath_bidsconvert_example(varargin)

% function that contains local data path, is ignored in .gitignore

if ~isempty(varargin{1})
    if isstruct(varargin{1})
        
        if sum(contains(fieldnames(varargin{1}),'sub_labels'))
            if contains(varargin{1}.sub_labels,'RESP')
                
                % SPES
                foldername = input('Choose SystemPlus-folder: testomgeving, RESPect_spes_scratch, RESPect_chronic_ECoG_trc: ','s');
                if strcmp(foldername,'testomgeving')
                    cfg(1).proj_dirinput = '/folder/to/ieeg-files/testomgeving/patients/';
                elseif strcmp(foldername,'RESPect_spes_scratch')
                    cfg(1).proj_dirinput = '/folder/to/ieeg-files/RESPect_spes_scratch/patients/';
                elseif strcmp(foldername,'RESPect_chronic_ECoG_trc')
                    cfg(1).proj_dirinput = '/folder/to/ieeg-files/RESPect_chronic_ECoG_trc/patients/';
                else
                    error('Foldername is not recognized')
                end
            end
            
            cfg(2).proj_dirinput = '/folder/to/ieeg-files/chronic_ECoG/';
            cfg(1).proj_diroutput = '/folder/to/BIDS-files/chronic_ECoG/';
            cfg(2).proj_diroutput = '/folder/to/BIDS-files//CCEP/'; % optional: this could remain empty
            
        elseif contains(varargin{1}.sub_labels,'REC2Stim')            % REC2Stim

            cfg(1).proj_dirinput = '/folder/to/ieeg-files/REC2Stim/patients/';
            cfg(2).proj_dirinput = '/folder/to/BIDS-files/REC2Stimstudy/';
            cfg(1).proj_diroutput = '/folder/to/BIDS-files/REC2Stimstudy/';
        
        elseif contains(varargin{1}.sub_labels,'PRIOS')            % prios study

            cfg(1).proj_dirinput = '/folder/to/ieeg-files/PRIOS_study/patients';
            cfg(2).proj_dirinput = '/folder/to/BIDS-files/PRIOS_study/';
            cfg(1).proj_diroutput = '/folder/to/BIDS-files/PRIOS_study/';
            
        end
        pat = [input('What is the PAT-folder in micromed database? [PAT_XXX] ','s'),'/'];
        cfg(1).pathname = fullfile(cfg(1).proj_dirinput,pat);        
        
    end
end


if contains(fieldnames(cfg),'no_fieldtrip')
    cfg.fieldtrip_folder  = '/folder/to/fieldtrip/';
    % copy the private folder in fieldtrip to somewhere else
    cfg.fieldtrip_private = '/folder/to/fieldtrip_private/';
    %%add those later to path to avoid errors with function 'dist'
    rmpath(cfg.fieldtrip_folder)
    rmpath(cfg.fieldtrip_private)
else
    fieldtrip_folder  = '/folder/to/fieldtrip/';
    % copy the private folder in fieldtrip to somewhere else
    fieldtrip_private = '/folder/to/fieldtrip_private/';
end

jsonlab_folder    = '/folder/to/jsonlab/';
addpath(fieldtrip_folder)
addpath(fieldtrip_private)
addpath(jsonlab_folder)
ft_defaults


end
