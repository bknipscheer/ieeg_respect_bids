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

%% 0) preparations - matlab
% This script has a part that should be run in a linux terminal, and part
% that can be run in matlab. The parts that should be run in a linux
% terminal have "run in linux terminal" in the section title.
clear
close all
clc

config_elecPositions

%% STEP 3: defacing MRI - RUN IN LINUX TERMINAL!

clc

% Find the original MRI and copy this to the folder you're using.
% Rename the original MRI to (run the line below and copy the printed line in the command window):
fprintf('\n ----- RENAME T1WEIGHTED MRI TO: ----- \n %s_%s_T1w.nii\n',...
    cfg.sub_labels{:},...
    cfg.ses_label);

% Right click in the folder with the original MRI and start Linux terminal.
% Copy the printed lines in the command window to deface the MRI in the linux terminal:
fprintf('\n ----- OPEN %ssourcedata/%s/%s/anat/  ----- \n ----- CLICK WITH RIGHT MOUSE AND OPEN LINUX TERMINAL -----\n ----- RUN LINE BELOW IN LINUX TERMINAL ----- \n mri_deface %s_%s_T1w.nii %s  %s %s_%s_proc-deface_T1w.nii\n',...
    cfg.home_directory,...
    cfg.sub_labels{:},...
    cfg.ses_label,...
    cfg.sub_labels{:},...
    cfg.ses_label,...
    cfg.path_talairach,...
    cfg.path_face,...
    cfg.sub_labels{:},...
    cfg.ses_label);
% this takes around 5 minutes

fprintf('\n ----- RUN LINE BELOW IN LINUX TERMINAL, OPEN DEFACED MRI TO CHECK DEFACING ----- \n mricron \n')

% Copy the defaced MRI to the/sub-RESPXXXX/ses-X/anat-folder
fprintf('\n ----- IF DEFACING WAS CORRECT, COPY DEFACED MRI TO %s -----\n', cfg.anat_directory)

%% STEP 4: run freesurfer to segment brain add Destrieux atlases - RUN IN LINUX TERMINAL!
clc

% Make a freesurfer folder
if exist(cfg.freesurfer_directory, 'dir')
%     fprintf('\n%s exists already\n',cfg.freesurfer_directory)
else
    mkdir(cfg.freesurfer_directory)
end

% Right click in the folder with the original MRI and start Linux terminal.
% Copy the printed lines in the command window into the linux terminal:
fprintf('\n ----- OPEN %ssourcedata/%s/%s/anat/  ----- \n ----- CLICK WITH RIGHT MOUSE AND OPEN LINUX TERMINAL ----- \n ----- RUN LINE BELOW IN LINUX TERMINAL ----- \nexport SUBJECTS_DIR=%s\n',...
    cfg.home_directory,...
    cfg.sub_labels{:},...
    cfg.ses_label,...
    cfg.freesurfer_directory)
% Copy the printed lines in the command window to run Freesurfer in the linux terminal:
fprintf('\n ----- RUN LINE BELOW IN LINUX TERMINAL ----- \nrecon-all -autorecon-all -s %s -i %s%s_%s_proc-deface_T1w.nii -cw256\n',...
    cfg.sub_labels{:},...
    cfg.anat_directory,...
    cfg.sub_labels{:},...
    cfg.ses_label)

% This takes up to 12 hours to run! In the end, you will see a subject
% folder in the freesurfer folder.

%% STEP5: generate surface (The Hull) to project electrodes to - RUN IN LINUX TERMINAL
clc
% only for ECoG, because this is necessary to correct for brain-shift.
% ECoG electrodes are projected to the hull.

% Freesurfer creates a file called /mri/ribbon.mgz and we want to convert
% this file to a .nii file so we can read it and use it to create the hull
% (a tight balloon of where the electrodes should be on the pre-op MRI)

% Right click in the freesurfer/mri-folder and start Linux terminal.
% Copy the printed lines in the command window into the linux terminal:
fprintf('\n ----- OPEN %smri ----- \n ----- CLICK WITH RIGHT MOUSE AND OPEN LINUX TERMINAL ----- \n ----- RUN LINE BELOW IN LINUX TERMINAL ----- \nmri_convert ribbon.mgz t1_class.nii\n',cfg.freesurfer_directory)

%% STEP 6: Create the hull - matlab
settings_hull = [13,... % setting for smoothing: default 13
                0.2]; % setting for threshold: default 0.3

get_mask_V3(cfg.sub_labels{:},... % subject name
    [cfg.freesurfer_directory,'mri/t1_class.nii'],... % freesurfer class file
    cfg.anat_directory,... % where you want to safe the file
    cfg.hemisphere,... % 'l' for left 'r' for right --> only use the hemisphere where electrodes are located
    settings_hull(1),...% setting for smoothing
    settings_hull(2)); % settings for  threshold
% the hull is saved as sub-RESPXXXX_surface1_13_03.img

%% STEP 7: check hull - RUN IN Linux TERMINAL
% type 'mricron'
% load the MRI
% put the hull as overlay on top of the mri
% check whether the hull looks like it matches the dura (should be a tight
% baloon around the grey matter)

fprintf('\n ----- RUN LINE BELOW IN LINUX TERMINAL, OPEN DEFACED MRI AND PUT HULL AS OVERLAY ON TOP ----- \n mricron \n \n ----- CHECK WHETHER THE HULL IS A TIGHT BALLOON AROUND THE CORTEX ----- \n')

%% STEP 8: select electrodes from ct - matlab
% the order in which you click electrodes does not matter. Just make sure
% you click all electrodes implanted!
clc
fprintf(' ----- OPEN THE CT-SCAN AND CLICK ON ALL ELECTRODES. YOU CAN CHECK WHETHER YOU HAVE ALL ELECTRODES BY CLICKING ON VIEW RESULT ----- \n')

ctmr
% view result
% save image: saves as nifti hdr and img files
% this is saved as electrodes1.hdr and electrodes1.img

%% STEP 9: sort unprojected electrodes - matlab

fprintf('------ OPEN THE ELECTRODES.TSV \n-----')
% open electrodes.tsv
[filename, pathname] = uigetfile('*.tsv;*.tsv','Select electroces.tsv file',cfg.elec_input);
tb_elecs = readtable(fullfile(pathname,filename),'FileType','text','Delimiter','\t');

% Make ieeg folder
if exist(cfg.ieeg_directory, 'dir')
    fprintf('%s exists already\n',cfg.ieeg_directory)
else
    mkdir(cfg.ieeg_directory)
end

fprintf('------ OPEN THE CLICKED ELECTRODES YOU SAVED IN THE PREVIOUS STEP \n-----')

% sort unprojected electrodes
cfg.saveFile = sprintf('%s%s_%s_electrodes_temp.mat',cfg.ieeg_directory,cfg.sub_labels{:},cfg.ses_label);
sortElectrodes(tb_elecs,cfg); % [electrode labels, folder to save]
fprintf('Matched electrodes are saved in %s\n',cfg.saveFile)
% loads img file with electrodes from previous step
% saves in electrodes_temp.mat;

%% STEP 10: plot electrodes 2 surface - matlab
% corrects for the brain shift - ONLY for ECoG
% 1xN STRIP: do not project any electrodes that are already close to the
%               surfaces, such as subtemporal and interhemispheric
% DEPTH: do not project

% open electrodes.tsv
[filename, pathname] = uigetfile('*.tsv;*.tsv','Select electroces.tsv file',cfg.elec_input);
tb_elecs = readtable(fullfile(pathname,filename),'FileType','text','Delimiter','\t');

% log_elec_incl = ~strcmp(tb_elecs.group,'other');
% tb_elecs = tb_elecs(log_elec_incl,:);

[filename, pathname] = uigetfile('*.mat','Select electrodes_temp.mat',cfg.elec_input);
load(fullfile(pathname,filename));

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
fprintf('(such as subtemporal, interhemispheric, depth)? \n')
fprintf('choose from %s %s %s %s %s %s', changroups{:});
skip_elec = input(': ','s');
skip_elec = strsplit(skip_elec,{', ',',',' '});
skip_elec = skip_elec(~cellfun(@isempty,skip_elec));

% determine format of each specific electrode strip/grid and correct for
% brainshift
if contains(json.iEEGElectrodeGroups,'seeg','IgnoreCase',1)
    disp('This session is seeg, so no correction for brain shift is needed')
    
    elecmatrix_shift = elecmatrix;
    
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
            num_elecs = num_elecs(~isnan(elecmatrix(num_elecs,1)) );
            
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
                     [cfg.anat_directory, cfg.sub_labels{:},'_surface1_',num2str(settings_hull(1)),'_0',num2str(settings_hull(2)*10),'.img'],... % hull we just created
                   cfg.anat_directory);
%                      [cfg.anat_directory, cfg.sub_labels{:},'_',cfg.ses_label,'_proc-deface_T1w.nii'],... % T1 file (mr.img for same image space with electrode positions)
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


%% STEP 11: save electrode positions, corrected for brain shift to electrodes.tsv - matlab

saveFile = replace(cfg.saveFile,'_temp.mat','.tsv');
writetable(tb_elecs, saveFile, 'Delimiter', 'tab', 'FileType', 'text');
fprintf('Electrode positions, corrected for brainshift, are saved in %s\n',saveFile)
% TODO: change NaN to n/a!!

%% STEP 12: Write electrode positions as numbers in a nifti - matlab
% This is not necessary, only if you want to do some extra checks or so.
% e.g. it can be nice to visualize the projected electrodes in MRIcron.
% to visualize in MRIcron: open the defaced MRI and add the surface_all.img
% as overlay

[output,els,els_ind,outputStruct] = position2reslicedImage(elecmatrix_shift,[cfg.anat_directory, cfg.sub_labels{:},'_',cfg.ses_label,'_proc-deface_T1w.nii']);

for filenummer=1:100
    save([cfg.ieeg_directory cfg.sub_labels{:} '_' cfg.ses_label,'_electrodes_surface_loc_all' int2str(filenummer) '.mat'],'elecmatrix_shift');
    outputStruct.fname=[cfg.anat_directory,cfg.sub_labels{:},'_',cfg.ses_label,'_electrodes_surface_all' int2str(filenummer) '.img' ];
    if ~exist(outputStruct.fname,'file')>0
        fprintf('----- SAVING %s ------ \n', outputStruct.fname);
        % save the data
        spm_write_vol(outputStruct,output);
        break
    end
end


%% STEP 13: convert freesurfer file to .gii - RUN IN LINUX TERMINAL
clc

% Make surface folder
if exist(cfg.surface_directory, 'dir')
%     fprintf('%s exists already\n',cfg.surface_directory)
else
    mkdir(cfg.surface_directory)
end

% Right click in the dataBIDS/derivatives/freesurfer/sub-,./surf-folder and
% start Linux terminal.
% Copy the printed lines in the command window into the linux terminal:
fprintf('\n ----- OPEN %s%s/surf/ ---- \n ---- CLICK WITH YOUR RIGHT MOUSE AND OPEN LINUX TERMINAL ----- \n ----- RUN LINE BELOW IN LINUX TERMINAL ----- \nmris_convert %sh.pial %sh.pial.gii\n',cfg.freesurfer_directory,cfg.sub_labels{1},cfg.hemisphere,cfg.hemisphere)

%% STEP 14: Convert the .gii coordinates to the MRI native space - matlab
% load the Freesurfer gifti (freesurfer coordinates)
g = gifti(fullfile(cfg.freesurfer_directory,'surf',[cfg.hemisphere,'h.pial.gii']));

% convert from freesurfer space to original space
% the transformation matrix is in the /freesurfer/sub/mri/orig.mgz file:
mri_orig = fullfile(cfg.freesurfer_directory,'mri','orig.mgz');
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
gifti_name = fullfile(cfg.surface_directory, ...
[cfg.sub_labels{:},'_',cfg.ses_label,'_T1w_pial.' cfg.hemisphere '.surf.gii']);

save(g,gifti_name,'Base64Binary')

disp('gifti converted to original space')

%% STEP 15: add labels of atlases to tsv-file - matlab

[tb_elecs_atlases, cfg.destrieux_labels, cfg.DKT_labels] = lookupAtlases(cfg,tb_elecs);

disp('Atlases added')

%% STEP 16: CHECK ATLAS WITH ELECTRODE POSITIONS - matlab

cfg.show_labels = 'yes';
cfg.change_size = 'no';
cfg.change_color = 'no';
cfg.view_atlas ='yes';
% cfg.atlas = 'Destrieux'; % [DKT/Destrieux]
cfg.atlas = 'DKT'; % [DKT/Destrieux]
cfg.view_elec ='yes';

check_atlas_elec_MRI(cfg,tb_elecs_atlases)

%% STEP 17: save electrodes.tsv, make electrodes descriptor, write coordsystem and add hemisphere to existing ieeg_json files

addpath(cfg.fieldtrip_folder) 
addpath(cfg.fieldtrip_private)
ft_defaults

writetable(tb_elecs_atlases, ...
    fullfile(cfg.ieeg_directory, ...
    [cfg.sub_labels{:} '_' cfg.ses_label '_electrodes.tsv' ]),...
    'Filetype','text','Delimiter','\t');

fprintf('Saved %s\n',fullfile(cfg.ieeg_directory, ...
    [cfg.sub_labels{:} '_' cfg.ses_label '_electrodes.tsv' ]))

% 6. create electrodes descriptor

create_elecDesc(cfg.home_directory,cfg)

% 7. write coordsystem

write_coordsystemJSON(cfg)

% 8. add hemisphere to ieeg_json files 

D = dir(cfg.ieeg_directory);
ieeg_json_filenums = contains({D(:).name},'_ieeg.json');

for i=1:size(D,1)
    if ieeg_json_filenums(i) == 1
        ieeg_json = read_json([D(i).folder '/' D(i).name] );
        if strcmpi(cfg.hemisphere,'r')
            ieeg_json.iEEGPlacementScheme = 'right';
        elseif strcmpi(cfg.hemisphere,'l')
            ieeg_json.iEEGPlacementScheme = 'left';
        elseif strcmpi(cfg.hemisphere,'l,r') ||strcmpi(cfg.hemisphere,'r,l')
            ieeg_json.iEEGPlacementScheme = 'left,right';
        end
        write_json([D(i).folder '/' D(i).name], ieeg_json)
    end
end

%% STEP 17: save everything to /Fridge/chronic_ECoG

cfg.home_directory = '/Fridge/chronic_ECoG/';
cfg.freesurfer_directory = sprintf('%sderivatives/freesurfer/%s/%s/',cfg.home_directory,cfg.sub_labels{:},cfg.ses_label);
cfg.anat_directory = sprintf('%s%s/%s/anat/',cfg.home_directory,cfg.sub_labels{:},cfg.ses_label);
cfg.ieeg_directory = sprintf('%s%s/%s/ieeg/',cfg.home_directory,cfg.sub_labels{:},cfg.ses_label);
cfg.surface_directory = sprintf('%sderivatives/surfaces/%s/%s/',cfg.home_directory,cfg.sub_labels{:},cfg.ses_label);
cfg.elec_input = sprintf('%s%s/%s/ieeg/',cfg.home_directory,cfg.sub_labels{:},cfg.ses_label);

% 5B. save electrodes.tsv 
writetable(tb_elecs_atlases, ...
    fullfile(cfg.ieeg_directory, ...
    [cfg.sub_labels{:} '_' cfg.ses_label '_electrodes.tsv' ]),...
    'Filetype','text','Delimiter','\t');

disp('Saved electrodes.tsv')

% 6. create electrodes descriptor

create_elecDesc(cfg.home_directory,cfg)

% 7. write coordsystem

write_coordsystemJSON(cfg)

% 8. add hemisphere to ieeg_json files 

D = dir(cfg.ieeg_directory);
ieeg_json_filenums = contains({D(:).name},'_ieeg.json');

for i=1:size(D,1)
    if ieeg_json_filenums(i) == 1
        ieeg_json = read_json([D(i).folder '/' D(i).name] );
        if strcmpi(cfg.hemisphere,'r')
            ieeg_json.iEEGPlacementScheme = 'right';
        elseif strcmpi(cfg.hemisphere,'l')
            ieeg_json.iEEGPlacementScheme = 'left';
        elseif strcmpi(cfg.hemisphere,'l,r') ||strcmpi(cfg.hemisphere,'r,l')
            ieeg_json.iEEGPlacementScheme = 'left,right';
        end
        write_json([D(i).folder '/' D(i).name], ieeg_json)
    end
end



%% FUNCTIONS
            
function json = read_json(filename)
ft_info('reading %s\n', filename);
if ft_hastoolbox('jsonlab', 3)
    json = loadjson(filename);
else
    fid = fopen(filename, 'r');
    str = fread(fid, [1 inf], 'char=>char');
    fclose(fid);
    json = jsondecode(str);
end
end

function write_json(filename, json)
json = remove_empty(json);
ft_info('writing %s\n', filename);
if ft_hastoolbox('jsonlab', 3)
    savejson('', json, filename);
else
    str = jsonencode(json);
    fid = fopen(filename, 'w');
    fwrite(fid, str);
    fclose(fid);
end
end

function s = remove_empty(s)
fn = fieldnames(s);
fn = fn(structfun(@isempty, s));
s = removefields(s, fn);
end





