%% check electrode positions and move to /Fridge/CCEP and /Fridge/chronic_ECoG
% author: Dorien van Blooijs
% date: September 2019

% addpath(genpath('/home/dorien/git_rep/Paper_Hermes_2010_JNeuroMeth/'))
% addpath(genpath('/home/dorien/git_rep/jsonlab/'))
% addpath(genpath('/home/dorien/git_rep/JSONio/'))
addpath(('/home/dorien/git_rep/ieeg_respect_bids/electrode_positions/'))
addpath(genpath('/home/dorien/git_rep/BasicCode_ECoG_DvB/'))

dataBase.sub_label = 'sub-RESP0295';
dataBase.ses_label = 'ses-1';
cfg.hemi_cap = 'L';
cfg.show_labels = 'yes';
cfg.change_size = 'no';
cfg.change_color = 'no';

visualizeMRI_ECoG(dataBase,cfg)

