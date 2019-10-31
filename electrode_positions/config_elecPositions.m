%% settings for electrode positions

% Make sure SPM functions are in your Matlab path

addpath(genpath('/home/dorien/git_rep/ieeg_respect_bids/electrode_positions/'))
cfg.fieldtrip_folder  = '/home/dorien/git_rep/fieldtrip/';
% copy the private folder in fieldtrip to somewhere else
cfg.fieldtrip_private = '/home/dorien/git_rep/fieldtrip_private/'; % CHECK: this might give difficulties with function 'dist' 
% addpath(cfg.fieldtrip_folder) 
% addpath(cfg.fieldtrip_private)
% ft_defaults

cfg.home_directory = '/Fridge/CCEP/';
cfg.path_talairach = '/Fridge/users/dorien/MRI_defaced/talairach_mixed_with_skull.gca';
cfg.path_face = '/Fridge/users/dorien/MRI_defaced/face.gca';
cfg.sub_labels = {['sub-' input('Patient number (RESPXXXX): ','s')]};
cfg.ses_label = input('Session number (ses-X): ','s');
cfg.hemisphere = input('Hemisphere with implanted electrodes [l/r]: ','s');
% cfg.freesurfer_directory = '/Fridge/users/dorien/dataBIDS/derivatives/freesurfer/';
cfg.freesurfer_directory = sprintf('%sderivatives/freesurfer/%s/%s/',cfg.home_directory,cfg.sub_labels{:},cfg.ses_label);
% cfg.anat_directory = sprintf('/Fridge/users/dorien/dataBIDS/%s/%s/anat/',cfg.sub_labels{:},cfg.ses_label);
cfg.anat_directory = sprintf('%s%s/%s/anat/',cfg.home_directory,cfg.sub_labels{:},cfg.ses_label);
cfg.ieeg_directory = sprintf('%s%s/%s/ieeg/',cfg.home_directory,cfg.sub_labels{:},cfg.ses_label);
cfg.surface_directory = sprintf('%sderivatives/surfaces/%s/%s/',cfg.home_directory,cfg.sub_labels{:},cfg.ses_label);
cfg.elec_input = sprintf('%s%s/%s/ieeg/',cfg.home_directory,cfg.sub_labels{:},cfg.ses_label);
