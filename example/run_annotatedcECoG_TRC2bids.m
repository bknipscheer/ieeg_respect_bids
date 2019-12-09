%% example cecog TRC to bids
% make sure you have the toolboxes jsonlab and fieldtrip available

clear all
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
ft_defaults
clear all

%% TRC to bids - run all files in patient-folder
% CLE
% cfg(1).proj_dirinput = '/home/dorien/Desktop/bulk/smb-share:server=smb-ds.bsc01.gd.umcutrecht.nl,share=ds_her_respect-leijten/Dorien/c_ecog/sz_cle/RESPect_sz_scratch/patients';
% cfg(1).proj_diroutput = '/Fridge/chronic_ECoG';%'/Fridge/users/Dorien/sz_cle';

% SPES
cfg(1).proj_dirinput = '/home/dorien/Desktop/bulk/smb-share:server=smb-ds.bsc01.gd.umcutrecht.nl,share=ds_her_respect-leijten/Dorien/c_ecog/spes/RESPect_spes_scratch/patients';
cfg(2).proj_dirinput = '/Fridge/chronic_ECoG'; 
cfg(1).proj_diroutput = '/Fridge/chronic_ECoG'; 
cfg(2).proj_diroutput = '/Fridge/CCEP'; % optional: this could remain empty

[~, pathname] = uigetfile('*.TRC;*.trc','Select *.TRC file',[cfg(1).proj_dirinput]);
files = dir(pathname);
runall = struct;

% run all files within your input directory
for i=1:size(files,1) 
    runall(i).file = files(i).name;
    if contains(files(i).name,'EEG_')
        
        cfg(1).filename = [pathname,files(i).name];
        
        pathsplit = strsplit(pathname,{'/'});
        patient = pathsplit{end-1};
        filesplit = strsplit(files(i).name,{'_','.TRC'});
        file = filesplit{end-1};
        
        fprintf('Running %s, writing EEG: %s to BIDS \n', patient,file)
        [runall(i).status,runall(i).msg,runall(i).output,runall(i).metadata,runall(i).annots] = annotatedTRC2bids(cfg);
    end
end

if any([runall(:).status])
    disp('All runs are done, but some still have errors. Fix them manually!')
else
    disp('All runs are completed')
end

%% run files which gave errors again 

for i=1:size(runall,2)
    
    if runall(i).status ==1 
        cfg(1).filename = [pathname,runall(i).file];
        
        pathsplit = strsplit(pathname,{'/'});
        patient = pathsplit{end-1};
        filesplit = strsplit(runall(i).file,{'_','.TRC'});
        file = filesplit{end-1};
        
        fprintf('Running %s, writing EEG: %s to BIDS \n', patient,file)
        [runall(i).status,runall(i).msg,runall(i).output,runall(i).metadata,runall(i).annots] = annotatedTRC2bids(cfg);

    end
end

if any([runall(:).status])
    disp('All runs are done, but some still have errors')
else
    disp('All runs are completed')
end

%% run one single file instead of all files within the input directory

% CLE
% cfg(1).proj_dirinput = '/home/dorien/Desktop/bulk/smb-share:server=smb-ds.bsc01.gd.umcutrecht.nl,share=ds_her_respect-leijten/Dorien/c_ecog/sz_cle/RESPect_sz_scratch/patients';
% cfg(1).proj_diroutput = '/Fridge/chronic_ECoG';%'/Fridge/users/Dorien/sz_cle';

% SPES
cfg(1).proj_dirinput = '/home/dorien/Desktop/bulk/smb-share:server=smb-ds.bsc01.gd.umcutrecht.nl,share=ds_her_respect-leijten/Dorien/c_ecog/spes/RESPect_spes_scratch/patients';
cfg(2).proj_dirinput = '/Fridge/chronic_ECoG'; 
cfg(1).proj_diroutput = '/Fridge/chronic_ECoG';
cfg(2).proj_diroutput = '/Fridge/CCEP'; % optional: this could remain empty

[filename, pathname] = uigetfile('*.TRC;*.trc','Select *.TRC file',[cfg(1).proj_dirinput]);
        
cfg(1).filename = [pathname,filename];

pathsplit = strsplit(pathname,{'/'});
patient = pathsplit{end-1};
filesplit = strsplit(filename,{'_','.TRC'});
file = filesplit{end-1};

fprintf('Running %s, writing EEG: %s to BIDS \n', patient,file)
[status,msg,output,metadata,annots] = annotatedTRC2bids(cfg);


