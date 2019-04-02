%% example cecog TRC to bids
% make sure you have the toolboxes jsonlab and fieldtrip available

addpath('Desktop/git_rep/ieeg_respect_bids/trc2bids/')
addpath('Desktop/git_rep/ieeg_respect_bids/micromed_utils/')
addpath('Desktop/git_rep/ieeg_respect_bids/external/')

fieldtrip_folder  = '/home/dorien/Desktop/git_rep/fieldtrip/';
% copy the private folder in fieldtrip to somewhere else
fieldtrip_private = '/home/dorien/Desktop/git_rep/fieldtrip_private/';
jsonlab_folder    = '/home/dorien/Desktop/git_rep/jsonlab/';
addpath(fieldtrip_folder) 
addpath(fieldtrip_private)
addpath(jsonlab_folder)

%% TRC to bids

cfg.proj_dirinput = '/home/dorien/Desktop/bulk/smb-share:server=smb-ds.bsc01.gd.umcutrecht.nl,share=ds_her_respect-leijten/Dorien/c_ecog/spes/RESPect_spes_scratch/patients';
cfg.proj_diroutput = '/Fridge/CCEP';

[filename, pathname] = uigetfile('*.TRC;*.trc','Select *.TRC file',[cfg.proj_dirinput]);
cfg.filename = [pathname,filename];

pathsplit = strsplit(pathname,{'/'});
patient = pathsplit{end-1};
filesplit = strsplit(filename,{'_','.TRC'});
file = filesplit{end-1};

fprintf('Running %s, writing EEG: %s to BIDS \n', patient,file)

[status,msg,output,metadata,annots] = annotatedTRC2bids(cfg)

