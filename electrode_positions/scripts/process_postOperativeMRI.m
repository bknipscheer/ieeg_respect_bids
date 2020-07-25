% code to display electrodes coregistered in pre-operative T1 weighted MRI with
% post-operative T1 weighted MRI

%   Created by:
%     Copyright (C) 2019  D. van Blooijs, Dept of Neurology and Neurosurgery, University Medical Center Utrecht

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


% This script has a part that should be run in a linux terminal, and part
% that can be run in matlab. The parts that should be run in a linux
% terminal have "run in linux terminal" in the section title.

%% patient characteristics - matlab

clear 
cfg(1).sub_labels = {['sub-' input('Patient number (RESPXXXX)/(REC2StimXX)/(PRIOSXX): ','s')]};
cfg(1).no_fieldtrip = 'yes';
cfg(1).mode = 'electrodeposition_postMRI';

% set paths
cfg = setLocalDataPath(cfg);

%% STEP 3: defacing MRI - RUN IN LINUX TERMINAL!

clc

% Find the original MRI and copy this to the folder you're using.
% Rename the original MRI to (run the line below and copy the printed line in the command window):
fprintf('\n ----- RENAME T1WEIGHTED MRI TO: ----- \n %s_%s_T1w_post.nii\n',...
    cfg(1).sub_labels{:},...
    cfg(1).ses_label);

% Right click in the folder with the original MRI and start Linux terminal.
% Copy the printed lines in the command window to deface the MRI in the linux terminal:
fprintf('\n ----- OPEN %ssourcedata/%s/%s/anat/  ----- \n ----- CLICK WITH RIGHT MOUSE AND OPEN LINUX TERMINAL -----\n ----- RUN LINE BELOW IN LINUX TERMINAL ----- \n mri_deface %s_%s_T1w_post.nii %s  %s %s_%s_proc-deface_T1w_post.nii\n',...
    cfg(1).proj_diroutput,...
    cfg(1).sub_labels{:},...
    cfg(1).ses_label,...
    cfg(1).sub_labels{:},...
    cfg(1).ses_label,...
    cfg(1).path_talairach,...
    cfg(1).path_face,...
    cfg(1).sub_labels{:},...
    cfg(1).ses_label);
% this takes around 5 minutes

fprintf('\n ----- RUN LINE BELOW IN LINUX TERMINAL, OPEN DEFACED MRI TO CHECK DEFACING ----- \n mricron \n')

%% STEP 4: run freesurfer to segment brain add Destrieux atlases - RUN IN LINUX TERMINAL!
clc

% Copy defaced .nii to correct folder
if exist(cfg(1).anat_directory,'dir')
    copyfile(fullfile(cfg(1).proj_diroutput,'sourcedata',cfg(1).sub_labels{:},cfg(1).ses_label,'anat',...
        [cfg(1).sub_labels{:},'_',cfg(1).ses_label,'_proc-deface_T1w_post.nii']),[cfg(1).anat_directory,cfg(1).sub_labels{:},'_',cfg(1).ses_label,'_proc-deface_T1w_post.nii'])
else
    mkdir(cfg(1).anat_directory)
        copyfile(fullfile(cfg(1).proj_diroutput,'sourcedata',cfg(1).sub_labels{:},cfg(1).ses_label,'anat',...
        [cfg(1).sub_labels{:},'_',cfg(1).ses_label,'_proc-deface_T1w_post.nii']),[cfg(1).anat_directory,cfg(1).sub_labels{:},'_',cfg(1).ses_label,'_proc-deface_T1w_post.nii'])
end

% Make a freesurfer folder
if exist(cfg(1).freesurfer_directory, 'dir')
%     fprintf('\n%s exists already\n',cfg(1).freesurfer_directory)
else
    mkdir(cfg(1).freesurfer_directory)
end

% Right click in the folder with the original MRI and start Linux terminal.
% Copy the printed lines in the command window into the linux terminal:
fprintf('\n ----- OPEN %ssourcedata/%s/%s/anat/  ----- \n ----- CLICK WITH RIGHT MOUSE AND OPEN LINUX TERMINAL ----- \n ----- RUN LINE BELOW IN LINUX TERMINAL ----- \nexport SUBJECTS_DIR=%s\n',...
    cfg(1).proj_diroutput,...
    cfg(1).sub_labels{:},...
    cfg(1).ses_label,...
    cfg(1).freesurfer_directory)
% Copy the printed lines in the command window to run Freesurfer in the linux terminal:
fprintf('\n ----- RUN LINE BELOW IN LINUX TERMINAL ----- \nrecon-all -autorecon-all -s %s -i %s%s_%s_proc-deface_T1w_post.nii -cw256\n',...
    cfg(1).sub_labels{:},...
    cfg(1).anat_directory,...
    cfg(1).sub_labels{:},...
    cfg(1).ses_label)

% This takes up to 12 hours to run! In the end, you will see a subject
% folder in the freesurfer folder.


%% STEP 13: convert freesurfer file to .gii - RUN IN LINUX TERMINAL
clc

% Make surface folder
if exist(cfg(1).surface_directory, 'dir')
%     fprintf('%s exists already\n',cfg(1).surface_directory)
else
    mkdir(cfg(1).surface_directory)
end

% Right click in the dataBIDS/derivatives/freesurfer/sub-,./surf-folder and
% start Linux terminal.
% Copy the printed lines in the command window into the linux terminal:
fprintf('\n ----- OPEN %ssurf/ ---- \n ---- CLICK WITH YOUR RIGHT MOUSE AND OPEN LINUX TERMINAL ----- \n ----- RUN LINE BELOW IN LINUX TERMINAL ----- \nmris_convert %sh.pial %sh.pial.gii\n',cfg(1).freesurfer_directory,cfg(1).hemisphere,cfg(1).hemisphere)

%% STEP 14: Convert the .gii coordinates to the MRI native space - matlab
% load the Freesurfer gifti (freesurfer coordinates)
g = gifti(fullfile(cfg(1).freesurfer_directory,'surf',[cfg(1).hemisphere,'h.pial.gii']));

% convert from freesurfer space to original space
% the transformation matrix is in the /freesurfer/sub/mri/orig.mgz file:
mri_orig = fullfile(cfg(1).freesurfer_directory,'mri','orig.mgz');
orig = MRIread(mri_orig); % MRIread is a function from vistasoft
Torig = orig.tkrvox2ras;
Norig = orig.vox2ras;
freeSurfer2T1 = Norig*inv(Torig);

% convert freesurfer vertices to original T1 space
vert_mat = double(([g.vertices ones(size(g.vertices,1),1)])');
vert_mat = freeSurfer2T1*vert_mat;
vert_mat(4,:) = [];
vert_mat = vert_mat';
g.vertices = vert_mat; clear vert_mat

% save correct coordinates back as a gifti
gifti_name = fullfile(cfg(1).surface_directory, ...
[cfg(1).sub_labels{:},'_',cfg(1).ses_label,'_T1w_pial.' cfg(1).hemisphere '.surf.gii']);

save(g,gifti_name,'Base64Binary')

disp('gifti converted to original space')

%% STEP 16: CHECK ATLAS WITH ELECTRODE POSITIONS - matlab

fprintf('------ OPEN THE ELECTRODES.TSV \n-----')
% open electrodes.tsv
[filename, pathname] = uigetfile('*.tsv;*.tsv','Select electroces.tsv file',cfg(1).elec_input);
tb_elecs = readtable(fullfile(pathname,filename),'FileType','text','Delimiter','\t');

cfg(1).show_labels = 'yes';
cfg(1).change_size = 'no';
cfg(1).change_color = 'no';
cfg(1).view_atlas ='no';
% cfg(1).atlas = 'Destrieux'; % [DKT/Destrieux]
% cfg(1).atlas = 'DKT'; % [DKT/Destrieux]
cfg(1).view_elec ='yes';
cfg(1).elec_offset = 'no';
cfg(1).transparency = 0.6;

check_atlas_elec_MRI_pre_post(cfg(1),tb_elecs)








