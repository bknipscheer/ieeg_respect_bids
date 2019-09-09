% code to localize subdural eeg and stereo eeg electrodes
%   Created by:
%
%     Copyright (C) 2009  D. Hermes, Dept of Neurology and Neurosurgery, University Medical Center Utrecht
%                   2019  D. van Blooijs, Dept of Neurology and Neurosurgery, University Medical Center Utrecht

%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
%
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
%
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.

% DvB - made it usable with BIDS electrodes.tsv September 2019

%% 1) run segment anatomical MR in SPM5 <-- dit stuk kan denk ik gewoon weg? zie uitgeschreven stuk hieronder
% startup SPM (spm functions are used)
% coregister + reslice CT to anatomical MR
% in SPM5 reference image = MR, source image = CT
% segment MR in freesurfer


%% 0) preparations
% This script has a part that should be run in a linux terminal, and part
% that can be run in matlab. The parts that should be run in a linux
% terminal have "run in linux terminal" in the section title. 

% Make sure SPM functions are in your path

addpath(genpath('/home/dorien/git_rep/Paper_Hermes_2010_JNeuroMeth/'))
addpath(('/home/dorien/git_rep/ieeg_respect_bids/electrode_positions/'))

% cfg.dataPath = '/Fridge/chronic_ECoG';
cfg.path_talairach = '/Fridge/users/dorien/MRI_defaced/talairach_mixed_with_skull.gca ';
cfg.path_face = '/Fridge/users/dorien/MRI_defaced/face.gca';
cfg.freesurfer_directory = '/Fridge/users/dorien/ccep/dataBIDS/derivatives/freesurfer/';
cfg.sub_labels = {['sub-' input('Patient number (RESPXXXX): ','s')]};
cfg.ses_label = input('Session number (ses-X): ','s');
cfg.anat_directory = sprintf('/Fridge/users/dorien/ccep/dataBIDS/%s/%s/anat/',cfg.sub_labels{:},cfg.ses_label);
cfg.ieeg_directory = sprintf('/Fridge/users/dorien/ccep/dataBIDS/%s/%s/ieeg/',cfg.sub_labels{:},cfg.ses_label);

%% defacing MRI - RUN IN LINUX TERMINAL!

% Find the original MRI and copy this to the folder you're using.
% Rename the original MRI to (run the line below and copy the printed line in the command window):
fprintf('%s_%s_proc_deface_T1w.nii\n',...
    cfg.sub_labels{:},...
    cfg.ses_label);

% Right click in the folder with the original MRI and start Linux terminal.
% Copy the printed lines in the command window to deface the MRI in the linux terminal:
fprintf('mri_deface %s_%s_T1w.nii %s  %s %s_%s_proc_deface_T1w.nii\n',...
    cfg.sub_labels{:},...
    cfg.ses_label,...
    cfg.path_talairach,...
    cfg.path_face,...
    cfg.sub_labels{:},...
    cfg.ses_label);

%% run freesurfer to add Destrieux atlases - RUN IN LINUX TERMINAL!

% Make a freesurfer folder
if exist(cfg.freesurfer_directory, 'dir')
    fprintf('%s exists already\n',cfg.freesurfer_directory)
else
    mkdir(cfg.freesurfer_directory)
end

% Right click in the folder with the original MRI and start Linux terminal.
% Copy the printed lines in the command window into the linux terminal:
fprintf('export SUBJECTS_DIR=%s\n',cfg.freesurfer_directory)
% Copy the printed lines in the command window to run Freesurfer in the linux terminal:
fprintf('recon-all -autorecon-all -s %s -i %s%s_%s_proc-deface_T1w.nii -cw256\n',...
    cfg.sub_labels{:},...
    cfg.anat_directory,...
    cfg.sub_labels{:},...
    cfg.ses_label)

% This takes up to 12 hours to run! In the end, you will see a subject
% folder in the freesurfer folder.


%% 2) generate surface (The Hull) to project electrodes to
% only for ECoG

% get_mask(subject,gray,white,outputdir,degree of smoothing,threshold) 
% e.g. get_mask(6,0.1) or get_mask(16,0.3)

% if using freesurfer: 
get_mask_V3(cfg.sub_labels{:},... % subject name
    './data/t1_class.nii',... % freesurfer class file
    './',... % where you want to safe the file
    'r',... % 'l' for left 'r' for right
    13,0.3); % settings for smoothing and threshold


%% 3) select electrodes from ct
ctmr
% view result
% save image: saves as nifti hdr and img files

%% 4) sort unprojected electrodes
% open electrodes.tsv
elec_input = sprintf('/Fridge/CCEP/%s/%s/ieeg/',cfg.sub_labels{:},cfg.ses_label);
[filename, pathname] = uigetfile('*.tsv;*.tsv','Select electroces.tsv file',elec_input);
tb_elecs = readtable(fullfile(pathname,filename),'FileType','text','Delimiter','\t');

% Make ieeg folder
if exist(cfg.ieeg_directory, 'dir')
    fprintf('%s exists already\n',cfg.ieeg_directory)
else
    mkdir(cfg.ieeg_directory)
end

% sort unprojected electrodes
cfg.saveFile = sprintf('%s%s_%s_electrodes_temp.mat',cfg.ieeg_directory,cfg.sub_labels{:},cfg.ses_label);
sortElectrodes(tb_elecs,cfg); % [electrode labels, folder to save]
% loads img file with electrodes from previous step
% saves in electrodes_temp.tsv;

% TODO: change NaN to n/a!!

%% 5) plot electrodes 2 surface
% electrodes2surf(subject,localnorm index,do not project electrodes closer than 3 mm to surface)

%% 6) combine electrode files into one and make an image
% not necessary since we use electrodes.tsv directly
% elecmatrix=nan(126,3);
% 
% a = load('/home/dorien/Desktopelectrodes_loc1.mat'); %
% elecmatrix(1:64,:) = a.elecmatrix;
% a = load('/home/dorien/Desktopelectrodes_loc2.mat'); %
% elecmatrix(65:90,:) = a.elecmatrix;
% a = load('/home/dorien/Desktopelectrodes_loc3.mat'); %
% elecmatrix(97:126,:) = a.elecmatrix;
% 
% save('/home/dorien/Desktop/Fransen_elecmatrix','elecmatrix')

%% 7) ?
[output,els,els_ind,outputStruct]=position2reslicedImage(elecmatrix,'./data_freesurfer/name_t1.nii');

for filenummer=1:100
    save(['./data/' subject '_electrodes_surface_loc_all' int2str(filenummer) '.mat'],'elecmatrix');
    outputStruct.fname=['./data/electrodes_surface_all' int2str(filenummer) '.img' ];
    if ~exist(outputStruct.fname,'file')>0
        disp(['saving ' outputStruct.fname]);
        % save the data
        spm_write_vol(outputStruct,output);
        break
    end
end




