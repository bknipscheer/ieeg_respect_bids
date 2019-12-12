%% config_TRC2BIDS
% make sure you have the toolboxes jsonlab and fieldtrip available

%% set paths
addpath(genpath('git_rep/ieeg_respect_bids/trc2bids/'))
addpath('git_rep/ieeg_respect_bids/micromed_utils/')
addpath('git_rep/ieeg_respect_bids/external/')

fieldtrip_folder  = '/home/dorien/git_rep/fieldtrip/';
% copy the private folder in fieldtrip to somewhere else
fieldtrip_private = '/home/dorien/git_rep/fieldtrip_private/';
jsonlab_folder    = '/home/dorien/git_rep/jsonlab/';
addpath(fieldtrip_folder) 
addpath(fieldtrip_private)
addpath(jsonlab_folder)
ft_defaults

clear

%% folder settings and selecting the correct patient

% CLE
% cfg(1).proj_dirinput = '/home/dorien/Desktop/bulk/smb-share:server=smb-ds.bsc01.gd.umcutrecht.nl,share=ds_her_respect-leijten/Dorien/c_ecog/sz_cle/RESPect_sz_scratch/patients';
% cfg(1).proj_diroutput = '/Fridge/chronic_ECoG';%'/Fridge/users/Dorien/sz_cle';

% SPES
% cfg(1).proj_dirinput = '/home/dorien/Desktop/bulk/smb-share:server=smb-ds.bsc01.gd.umcutrecht.nl,share=ds_her_respect-leijten/Dorien/c_ecog/spes/RESPect_spes_scratch/patients';
cfg(1).proj_dirinput = '/home/dorien/Desktop/bulk/smb-share:server=smb-ds.bsc01.gd.umcutrecht.nl,share=ds_her_respect-leijten/RESPect_chronic_ECoG_trc/patients';
cfg(2).proj_dirinput = '/Fridge/chronic_ECoG'; 
cfg(1).proj_diroutput = '/Fridge/chronic_ECoG'; 
cfg(2).proj_diroutput = '/Fridge/CCEP'; % optional: this could remain empty

% [~, pathname] = uigetfile('*.TRC;*.trc','Select *.TRC file',[cfg(1).proj_dirinput]);
pat = [input('What is the PAT-folder in micromed database? [PAT_XXX] ','s'),'/'];
pathname = fullfile(cfg(1).proj_dirinput,pat);

