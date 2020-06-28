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

cfg(1).sub_labels = {['sub-' input('Patient number (RESPXXXX)/(REC2StimXX)/(PRIOSXX): ','s')]};

% CLE
% cfg(1).proj_dirinput = '/home/dorien/Desktop/bulk/smb-share:server=smb-ds.bsc01.gd.umcutrecht.nl,share=ds_her_respect-leijten/Dorien/c_ecog/sz_cle/RESPect_sz_scratch/patients';
% cfg(1).proj_diroutput = '/Fridge/chronic_ECoG';%'/Fridge/users/Dorien/sz_cle';

if contains(cfg(1).sub_labels,'RESP')
    % SPES
    foldername = input('Choose SystemPlus-folder: testomgeving, RESPect_spes_scratch, RESPect_chronic_ECoG_trc: ','s');
    if strcmp(foldername,'testomgeving')
        cfg(1).proj_dirinput = '/home/dorien/Desktop/db_respect/Dorien/testomgeving/patients';
    elseif strcmp(foldername,'RESPect_spes_scratch')
        cfg(1).proj_dirinput = '/home/dorien/Desktop/db_respect/Dorien/c_ecog/spes/RESPect_spes_scratch/patients';
    elseif strcmp(foldername,'RESPect_chronic_ECoG_trc')
        cfg(1).proj_dirinput = '/home/dorien/Desktop/bulkstorage/remote_conn/smb-share:server=arch11-smb-ds.arch11.gd.umcutrecht.nl,share=her$/snap/Respect-Leijten/RESPect_chronic_ECoG_trc/patients';
    else
       error('Foldername is not recognized') 
    end
    cfg(2).proj_dirinput = '/Fridge/KNF/chronic_ECoG';
    cfg(1).proj_diroutput = '/Fridge/KNF/chronic_ECoG';
    cfg(2).proj_diroutput = '/Fridge/KNF/CCEP'; % optional: this could remain empty
    
elseif contains(cfg(1).sub_labels,'REC2Stim')
    % REC2Stim
    cfg(1).proj_dirinput = '/home/dorien/Desktop/db/Dorien/REC2Stim/patients';
    cfg(2).proj_dirinput = '/Fridge/KNF/REC2Stimstudy';
    cfg(1).proj_diroutput = '/Fridge/KNF/REC2Stimstudy';
elseif contains(cfg(1).sub_labels,'PRIOS')
    % prios study
    cfg(1).proj_dirinput = '/home/dorien/Desktop/bulkstorage/db/respect-leijten/PRIOS_study/patients';
    cfg(2).proj_dirinput = '/Fridge/KNF/REC2Stimstudy/PRIOS_study';
    cfg(1).proj_diroutput = '/Fridge/KNF/REC2Stimstudy/PRIOS_study';
    
end

% [~, pathname] = uigetfile('*.TRC;*.trc','Select *.TRC file',[cfg(1).proj_dirinput]);
pat = [input('What is the PAT-folder in micromed database? [PAT_XXX] ','s'),'/'];
pathname = fullfile(cfg(1).proj_dirinput,pat);

