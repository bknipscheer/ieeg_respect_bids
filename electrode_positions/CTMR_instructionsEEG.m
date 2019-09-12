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

%% 0) preparations
% This script has a part that should be run in a linux terminal, and part
% that can be run in matlab. The parts that should be run in a linux
% terminal have "run in linux terminal" in the section title.

% Make sure SPM functions are in your Matlab path

addpath(genpath('/home/dorien/git_rep/Paper_Hermes_2010_JNeuroMeth/'))
addpath(('/home/dorien/git_rep/ieeg_respect_bids/electrode_positions/'))
addpath(genpath('/home/dorien/git_rep/jsonlab/'))

% cfg.dataPath = '/Fridge/chronic_ECoG';
cfg.path_talairach = '/Fridge/users/dorien/MRI_defaced/talairach_mixed_with_skull.gca';
cfg.path_face = '/Fridge/users/dorien/MRI_defaced/face.gca';
cfg.freesurfer_directory = '/Fridge/users/dorien/dataBIDS/derivatives/freesurfer/';
cfg.sub_labels = {['sub-' input('Patient number (RESPXXXX): ','s')]};
cfg.ses_label = input('Session number (ses-X): ','s');
cfg.hemisphere = input('Hemisphere with implanted electrodes [l/r]: ','s');
cfg.anat_directory = sprintf('/Fridge/users/dorien/dataBIDS/%s/%s/anat/',cfg.sub_labels{:},cfg.ses_label);
cfg.ieeg_directory = sprintf('/Fridge/users/dorien/dataBIDS/%s/%s/ieeg/',cfg.sub_labels{:},cfg.ses_label);
cfg.elec_input = sprintf('/Fridge/CCEP/%s/%s/ieeg/',cfg.sub_labels{:},cfg.ses_label);

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

%% run freesurfer to segment brain add Destrieux atlases - RUN IN LINUX TERMINAL!

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
% only for ECoG, because this is necessary to correct for brain-shift.
% ECoG electrodes are projected to the hull.

% Freesurfer creates a file called /mri/ribbon.mgz and we want to convert
% this file to a .nii file so we can read it and use it to create the hull
% (a tight balloon of where the electrodes should be on the pre-op MRI)

% RUN IN LINUX TERMINAL:
% Right click in the freesurfer/mri-folder and start Linux terminal.
% Copy the printed lines in the command window into the linux terminal:
fprintf('mri_convert ribbon.mgz t1_class.nii\n')

settings_hull = [13,0.3];

% Go back to Matlab and create the hull
get_mask_V3(cfg.sub_labels{:},... % subject name
    [cfg.freesurfer_directory,cfg.sub_labels{:},'/mri/t1_class.nii'],... % freesurfer class file
    cfg.anat_directory,... % where you want to safe the file
    cfg.hemisphere,... % 'l' for left 'r' for right --> only use the hemisphere where electrodes are located
    settings_hull(1),...% setting for smoothing
    settings_hull(2)); % settings for  threshold
% the hull is saved as sub-RESPXXXX_surface1_13_03.img

%% check hull - RUN IN Linux TERMINAL
% type 'mricron'
% load the MRI
% put the hull as overlay on top of the mri
% check whether the hull looks like it matches the dura (should be a tight
% baloon around the grey matter)

%% 3) select electrodes from ct
% the order in which you click electrodes does not matter. Just make sure
% you click all electrodes implanted!

ctmr
% view result
% save image: saves as nifti hdr and img files
% this is saved as electrodes1.hdr and electrodes1.img

%% 4) sort unprojected electrodes
% open electrodes.tsv
[filename, pathname] = uigetfile('*.tsv;*.tsv','Select electroces.tsv file',cfg.elec_input);
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
fprintf('Matched electrodes are saved in %s\n',cfg.saveFile)
% loads img file with electrodes from previous step
% saves in electrodes_temp.mat;

%% 5) plot electrodes 2 surface
% corrects for the brain shift - ONLY for ECoG
% 1xN STRIP: do not project any electrodes that are already close to the
%               surfaces, such as subtemporal and interhemispheric
% DEPTH: do not project

% open electrodes.tsv
[filename, pathname] = uigetfile('*.tsv;*.tsv','Select electroces.tsv file',cfg.elec_input);
tb_elecs = readtable(fullfile(pathname,filename),'FileType','text','Delimiter','\t');

% only select letters from channelname
letters = regexp(tb_elecs.name,'[a-z_A-Z]');
channame = cell(size(tb_elecs,1),1);
for chan = 1:size(tb_elecs,1)
    channame{chan,1} = tb_elecs.name{chan}(letters{chan});
end

% load json-file with formats of specific electrode strips/grids
files = dir(cfg.elec_input);
jsonfile = find(contains({files(:).name},'_ieeg.json')==1,1);
json = loadjson(fullfile(cfg.elec_input, files(jsonfile).name));
chansplit = strsplit(json.iEEGElectrodeGroups,{';','['});
changroups = chansplit(diff(contains(chansplit,']'))==1);

% electrodes on strip subtemporal or interhemispheric should be skipped
% because they are already located close to brain surface.
fprintf('Which electrodes should be skipped since they are already close to surfaces? \n')
fprintf('(such as subtemporal, interhemispheric)? \n')
fprintf('choose from %s %s %s %s %s %s', changroups{:});
skip_elec = input(': ','s');
skip_elec = strsplit(skip_elec,{', ',',',' '});
skip_elec = skip_elec(~cellfun(@isempty,skip_elec));

% determine format of each specific electrode strip/grid and correct for
% brainshift
if contains(json.iEEGElectrodeGroups,'seeg','IgnoreCase',1)
    disp('This session is seeg, so no correction for brain shift is needed')
    
elseif contains(json.iEEGElectrodeGroups,'ecog','IgnoreCase',1)
    
    format = strsplit(json.iEEGElectrodeGroups,';');
    elecmatrix_shift = NaN(size(elecmatrix));
    
    for i=1:size(format,2)
        if contains(format{i},'[')
            % find channame and format of this group of electrodes
            formatelec = strsplit(format{i},{'[','x',']'});
            formatelec = formatelec(~cellfun(@isempty,formatelec));
            
            % find which electrodes belong to this specific group
            num_elecs = find(strcmp(formatelec{1},channame)==1);
            
            % check whether all electrodes in group are included
            if str2double(formatelec{2}) * str2double(formatelec{3}) == numel(num_elecs)
            else
                disp('ERROR: mismatch between found electrodes and expected number of electrodes in format!')
            end
            
            % skip electrodes mentioned above
            if ~contains(format{i},skip_elec)
                % format of specific grid/strip
                if any(contains(formatelec,'1')) % it is a 1xN strip/depth electrode
                    settings = [0,2];
                elseif any(contains(formatelec,'2')) % it is a 2xN strip
                    settings = [4,1];
                else % it is a grid (larger than 2xN)
                    settings = [5,1];
                end
                % correct location of electrodes for brain shift
                [out_els,out_els_ind] = electrodes2surf(...
                    cfg.sub_labels{:},... % subject name
                    settings(1),... % 5 for grid, 4 for 2xN strip, 0 for 1xN strip
                    settings(2),... % 1 for grid or 2xN strip, 2 for 1xN strip
                    num_elecs,... % matrix indices (rows) in elecmatrix, e.g. 1:64 is for grid C1-C64
                    [cfg.ieeg_directory, cfg.sub_labels{:},'_',cfg.ses_label,'_electrodes_temp.mat'],... % file that contains elecmatrix
                    [cfg.anat_directory, cfg.sub_labels{:},'_surface1_',num2str(settings_hull(1)),'_0',num2str(settings_hull(2)*10),'.img'],... % hull we just created
                    [cfg.anat_directory, cfg.sub_labels{:},'_',cfg.ses_label,'_proc-deface_T1w.nii'],... % T1 file (mr.img for same image space with electrode positions)
                    cfg.anat_directory);
                
                % saves automatically a matrix with projected electrode positions and an image
                % with projected electrodes
                % saves as electrodes_onsurface_filenumber_inputnr2
                elecmatrix_shift(num_elecs,:) = out_els;
                
            elseif contains(format{i},skip_elec)
                elecmatrix_shift(num_elecs,:) = elecmatrix(num_elecs,:);
            end

        end
    end
end

tb_elecs.x = elecmatrix_shift(:,1);
tb_elecs.y = elecmatrix_shift(:,2);
tb_elecs.z = elecmatrix_shift(:,3);

%% save electrode positions, corrected for brain shift to electrodes.tsv

saveFile = replace(cfg.saveFile,'_temp.mat','.tsv');
writetable(tb_elecs, saveFile, 'Delimiter', 'tab', 'FileType', 'text');
fprintf('Electrode positions, corrected for brainshift, are saved in %s\n',saveFile)
% TODO: change NaN to n/a!!

%% 6) combine electrode files into one and make an image
% Only necessary if you did step 5), because each grid/strip is projected
% separately...
% The separate projection can maybe be fixed by using code from iElvis or
% the Dykstra method...

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

%% 7) Write electrode positions as numbers in a nifti
% This is not necessary, only if you want to do some extra checks or so.
% e.g. it can be nice to visualize the projected electrodes in MRIcron.
[output,els,els_ind,outputStruct] = position2reslicedImage(elecmatrix_shift,[cfg.anat_directory, cfg.sub_labels{:},'_',cfg.ses_label,'_proc-deface_T1w.nii']);

for filenummer=1:100
    save([cfg.ieeg_directory cfg.sub_labels{:} '_' cfg.ses_label,'_electrodes_surface_loc_all' int2str(filenummer) '.mat'],'elecmatrix_shift');
    outputStruct.fname=[cfg.anat_directory,cfg.sub_labels{:},'_',cfg.ses_label,'_electrodes_surface_all' int2str(filenummer) '.img' ];
    if ~exist(outputStruct.fname,'file')>0
        disp(['saving ' outputStruct.fname]);
        % save the data
        spm_write_vol(outputStruct,output);
        break
    end
end




