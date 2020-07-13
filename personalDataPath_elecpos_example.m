function cfg = personalDataPath_elecpos_example(varargin)

% function that contains local data path, is ignored in .gitignore

if ~isempty(varargin{1})
    if isstruct(varargin{1})
        
        if sum(contains(fieldnames(varargin{1}),'sub_labels'))
            if contains(varargin{1}.sub_labels,'RESP')
                cfg(1).proj_dirinput = '/folder/to/bids-files/CCEP/';
                cfg(2).proj_dirinput = '/folder/to/bids-files/chronic_ECoG/';
                cfg(1).proj_diroutput = '/folder/to/bids-files/CCEP/';
                cfg(2).proj_diroutput = '/folder/to/bids-files/chronic_ECoG/'; % optional: this could remain empty
                
            elseif contains(varargin{1}.sub_labels,'REC2Stim')
                % REC2Stim
                cfg(1).proj_dirinput = '/folder/to/ieeg-files/REC2Stim/patients/';
                cfg(2).proj_dirinput = '/folder/to/bids-files/REC2Stimstudy/';
                cfg(1).proj_diroutput = '/folder/to/bids-files/REC2Stimstudy/';
                
            elseif contains(varargin{1}.sub_labels,'PRIOS')
                % prios study
                cfg(1).proj_dirinput = '/folder/to/ieeg-files/PRIOS_study/patients';
                cfg(2).proj_dirinput = '/folder/to/bids-files/PRIOS_study/';
                cfg(1).proj_diroutput = '/folder/to/bids-files/PRIOS_study/';
                
            end
            
            % for electrode positions
            cfg(1).sub_labels = varargin{1}.sub_labels;
            cfg(1).ses_label = input('Session number (ses-X): ','s');
            cfg(1).hemisphere = input('Hemisphere with implanted electrodes [l/r]: ','s');
            cfg(1).freesurfer_directory = sprintf('%sderivatives/freesurfer/%s/%s/',cfg(1).proj_diroutput,cfg(1).sub_labels{:},cfg(1).ses_label);
            cfg(1).anat_directory = sprintf('%s%s/%s/anat/',cfg(1).proj_diroutput,cfg(1).sub_labels{:},cfg(1).ses_label);
            cfg(1).ieeg_directory = sprintf('%s%s/%s/ieeg/',cfg(1).proj_diroutput,cfg(1).sub_labels{:},cfg(1).ses_label);
            cfg(1).surface_directory = sprintf('%sderivatives/surfaces/%s/%s/',cfg(1).proj_diroutput,cfg(1).sub_labels{:},cfg(1).ses_label);
            cfg(1).elec_input = sprintf('%s%s/%s/ieeg/',cfg(1).proj_diroutput,cfg(1).sub_labels{:},cfg(1).ses_label);
            cfg(1).path_talairach = '/folder/to/talairach_mixed_with_skull.gca';
            cfg(1).path_face = '/folder/to/face.gca';
            
        end
    end
end

if any(contains(fieldnames(varargin{1}),'no_fieldtrip'))
    cfg(1).fieldtrip_folder  = '/folder/to/fieldtrip/';
    % copy the private folder in fieldtrip to somewhere else
    cfg(1).fieldtrip_private = '/folder/to/fieldtrip_private/';
    
    %add those later to path to avoid errors with function 'dist'
    rmpath(cfg(1).fieldtrip_folder)
    rmpath(cfg(1).fieldtrip_private)
else
    cfg(1).fieldtrip_folder  = '/folder/tofieldtrip/';
    % copy the private folder in fieldtrip to somewhere else
    cfg(1).fieldtrip_private = '/folder/to/fieldtrip_private/';
    
    addpath(cfg(1).fieldtrip_folder)
    addpath(cfg(1).fieldtrip_private)
    ft_defaults
    
end

jsonlab_folder    = '/folder/to/jsonlab/';
addpath(jsonlab_folder)


end
