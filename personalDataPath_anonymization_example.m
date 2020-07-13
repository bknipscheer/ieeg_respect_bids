function cfg = personalDataPath_anonymization_example(varargin)

% function that contains local data path, is ignored in .gitignore

if ~isempty(varargin{1})
    if isstruct(varargin{1})
        if strcmp(varargin{1}.mode,'anonymization')
            
            cfg(1).proj_dirinput = '/home/dorien/Desktop/temp_ecog/';
            tempName = varargin{1}.sub_labels{:};
            
            % check whether RESP-number is entered correctly
            if strcmp(tempName,'') && ~isempty(respName)
                
            elseif contains(tempName,'RESP') || contains(tempName,'REC2Stim') || contains(tempName,'PRIOS')
                cfg.respName = tempName;
            else
                error('RESPect/REC2Stim name is not correct')
            end
            
        else
            if sum(contains(fieldnames(varargin{1}),'sub_labels'))
                if contains(varargin{1}.sub_labels,'RESP')
                    % for conversion trc-file to BIDS
                    if strcmp(varargin{1}.mode,'bidsconversion')
                        
                        % SPES
                        foldername = input('Choose SystemPlus-folder: testomgeving, RESPect_spes_scratch, RESPect_chronic_ECoG_trc: ','s');
                        if strcmp(foldername,'testomgeving')
                            cfg(1).proj_dirinput = '/home/dorien/Desktop/bulkstorage/db/respect-leijten/Dorien/testomgeving/patients/';
                        elseif strcmp(foldername,'RESPect_spes_scratch')
                            cfg(1).proj_dirinput = '/home/dorien/Desktop/bulkstorage/db/respect-leijten/Dorien/c_ecog/spes/RESPect_spes_scratch/patients/';
                        elseif strcmp(foldername,'RESPect_chronic_ECoG_trc')
                            cfg(1).proj_dirinput = '/home/dorien/Desktop/bulkstorage/db/respect-leijten/RESPect_chronic_ECoG_trc/patients/';
                        else
                            error('Foldername is not recognized')
                        end
                    end
                    
                    cfg(2).proj_dirinput = '/Fridge/KNF/chronic_ECoG/';
                    cfg(1).proj_diroutput = '/Fridge/KNF/chronic_ECoG/';
                    cfg(2).proj_diroutput = '/Fridge/KNF/CCEP/'; % optional: this could remain empty
                    
                elseif contains(varargin{1}.sub_labels,'REC2Stim')
                    % REC2Stim
                    cfg(1).proj_dirinput = '/home/dorien/Desktop/bulkstorage/db/respect-leijten/Dorien/REC2Stim/patients/';
                    cfg(2).proj_dirinput = '/Fridge/KNF/REC2Stimstudy/';
                    cfg(1).proj_diroutput = '/Fridge/KNF/REC2Stimstudy/';
                elseif contains(varargin{1}.sub_labels,'PRIOS')
                    % prios study
                    cfg(1).proj_dirinput = '/home/dorien/Desktop/bulkstorage/db/respect-leijten/PRIOS_study/patients';
                    cfg(2).proj_dirinput = '/Fridge/KNF/PRIOS_study/';
                    cfg(1).proj_diroutput = '/Fridge/KNF/PRIOS_study/';
                    
                end
                
                % for conversion trc-file to BIDS
                if strcmp(varargin{1}.mode,'bidsconversion')
                    
                    pat = [input('What is the PAT-folder in micromed database? [PAT_XXX] ','s'),'/'];
                    cfg(1).pathname = fullfile(cfg(1).proj_dirinput,pat);
                    
                elseif strcmp(varargin{1}.mode,'electrode_position')
                    % for electrode positions
                    cfg(1).sub_labels = varargin{1}.sub_labels;
                    cfg(1).ses_label = input('Session number (ses-X): ','s');
                    cfg(1).hemisphere = input('Hemisphere with implanted electrodes [l/r]: ','s');
                    cfg(1).freesurfer_directory = sprintf('%sderivatives/freesurfer/%s/%s/',cfg(1).proj_diroutput,cfg(1).sub_labels{:},cfg(1).ses_label);
                    cfg(1).anat_directory = sprintf('%s%s/%s/anat/',cfg(1).proj_diroutput,cfg(1).sub_labels{:},cfg(1).ses_label);
                    cfg(1).ieeg_directory = sprintf('%s%s/%s/ieeg/',cfg(1).proj_diroutput,cfg(1).sub_labels{:},cfg(1).ses_label);
                    cfg(1).surface_directory = sprintf('%sderivatives/surfaces/%s/%s/',cfg(1).proj_diroutput,cfg(1).sub_labels{:},cfg(1).ses_label);
                    cfg(1).elec_input = sprintf('%s%s/%s/ieeg/',cfg(1).proj_diroutput,cfg(1).sub_labels{:},cfg(1).ses_label);
                    cfg(1).path_talairach = '/Fridge/users/dorien/MRI_defaced/talairach_mixed_with_skull.gca';
                    cfg(1).path_face = '/Fridge/users/dorien/MRI_defaced/face.gca';
                end
                
            end
        end
    end
    
    cfg(1).elec_input = '/home/dorien/Desktop/bulkstorage/db/respect-leijten/Electrodes/';
    
    if contains(fieldnames(cfg),'no_fieldtrip')
        cfg.fieldtrip_folder  = '/home/dorien/git_rep/fieldtrip/';
        % copy the private folder in fieldtrip to somewhere else
        cfg.fieldtrip_private = '/home/dorien/git_rep/fieldtrip_private/';
        %%add those later to path to avoid errors with function 'dist'
        rmpath(cfg.fieldtrip_folder)
        rmpath(cfg.fieldtrip_private)
    else
        fieldtrip_folder  = '/home/dorien/git_rep/fieldtrip/';
        % copy the private folder in fieldtrip to somewhere else
        fieldtrip_private = '/home/dorien/git_rep/fieldtrip_private/';
    end
    
    jsonlab_folder    = '/home/dorien/git_rep/jsonlab/';
    addpath(fieldtrip_folder)
    addpath(fieldtrip_private)
    addpath(jsonlab_folder)
    ft_defaults
    
    
end
end