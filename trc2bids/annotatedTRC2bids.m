%  Convert annotated (see annotation scheme in docs) micromed file (.TRC) to Brain Imaging Data Structure (BIDS)
%  it generate all the required directory structure and files
%
%  cfg.proj_dir - directory name where to store the files
%  cfg.filename - name of the micromed file to convert
%
%  output structure with some information about the subject
%  output.subjName - name of the subject
%
%


% The following functions rely and take inspiration on fieltrip data2bids.m  function
% (https://github.com/fieldtrip/fieldtrip.git)
%
% fieldtrip toolbox should be on the path (plus the fieldtrip_private folder)
% (see http://www.fieldtriptoolbox.org/faq/matlab_does_not_see_the_functions_in_the_private_directory/)
%
% jsonlab toolbox
% https://github.com/fangq/jsonlab.git

% some external function to read micromed TRC files is used
% https://github.com/fieldtrip/fieldtrip/blob/master/fileio/private/read_micromed_trc.m
% copied in the external folder


%
%
%     Copyright (C) 2019 Matteo Demuru
%	  Copyright (C) 2019 Dorien van Blooijs
%     Copyright (C) 2019 Willemiek Zweiphenning
%
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
%     along with this program.  If not, see <https://www.gnu.org/licenses/>.

function [status,msg,output,metadata,annots] = annotatedTRC2bids(cfg)

try
    output.subjName = '';
    output.sitName  = '';
    msg = '';
    
    check_input(cfg,'proj_dirinput');
    check_input(cfg,'proj_diroutput');
    check_input(cfg,'filename');
    
    proj_dirinput = cell(1);
    proj_diroutput = cell(1);
    proj_diroutputcopy = cell(1);
    
    for n=1:size(cfg,2)
        if ~isempty(cfg(n).proj_dirinput)
            proj_dirinput{n}  = cfg(n).proj_dirinput;
        end
        if ~isempty(cfg(n).proj_diroutput)
            if n==1
                proj_diroutput{n} = cfg(n).proj_diroutput;
            elseif n>1
                proj_diroutputcopy{n-1} = cfg(n).proj_diroutput;
            end
        end
    end
    filename  = cfg(1).filename;
    
    [indir,fname,exte] = fileparts(filename);
    %create the subject level dir if not exist
    
    [header,data,data_time,trigger,annots] = read_TRC_HDR_DATA_TRIGS_ANNOTS(filename);
    
    if(isempty(header) || isempty(data) || isempty(data_time) || isempty(annots)) %|| isempty(trigger)
        error('TRC reading failed')  ;
    end
    ch_label = deblank({header.elec.Name}');
    sub_label = strcat('sub-',upper(deblank(header.name)));
    
    output.subjName = sub_label;
    
    [status,msg,metadata] = extract_metadata_from_annotations(header,annots,ch_label,trigger,sub_label,cfg);
    
    output.sitName = replace(strcat('ses-',deblank(metadata.ses_name)),' ','');
    output.runName = replace(strcat('run-',deblank(metadata.run_name),header.hour,header.min),' ','');
    output.taskName = replace(deblank(metadata.task_name),' ','');
    
    if(status==0)
        %% move trc with the proper naming and start to create the folder structure
        % for now for simplicity using the .trc even though is not one of the
        % allowed format (later it should be moved to source)
        
        %proj-dir/
        %   sub-<label>/
        %       ses-<label>/
        %           ieeg/
        %               sub-<label>_ses-<label>_task-<task_label>_ieeg.<allowed_extension>
        %               sub-<label>_ses-<label>_task-<task_label>_ieeg.json
        %               sub-<label>_ses-<label>_task-<task_label>_channels.tsv
        %               sub-<label>_ses-<label>_task-<task_label>_events.tsv
        %               sub-<label>_ses-<label>_electrodes.tsv
        %               sub-<label>_ses-<label>_coordsystem.json
        %
        %               sub-<label>_ses-<label>_photo.jpg
        
        task_label    = strcat('task-',replace(deblank(metadata.task_name),' ',''));
        if strfind(lower(task_label),'spes')~=0
            task_desc = 'No task, electrical stimulation is performed. Patient is resting with eyes open/closed. The latter is not always specified.';
        elseif strfind(lower(task_label),'rest') ~=0
            task_desc = 'Patient is resting with eyes open/closed. The latter is not always specified.';
        elseif strfind(lower(task_label),'sleep') ~=0
            task_desc = 'Patient is sleeping.';
        elseif strfind(lower(task_label),'slawtrans') ~=0
            task_desc = 'Patient is trying to fall asleep or is waking up.';
        elseif strfind(lower(task_label),'motor') ~=0
            task_desc = 'Patient is doing a motor task.';
        elseif strfind(lower(task_label),'esm') ~=0
            task_desc = 'Electrical stimulation mapping is performed to delineate functional areas.';
        elseif strfind(lower(task_label),'sens') ~=0
            task_desc = 'Patient is doing a sensing task.';
        elseif strfind(lower(task_label),'Language') ~= 0
            task_desc = 'Patient is doing a language task.';
        else
            task_desc = 'Not specified'
        end
        ses_label     = strcat('ses-',deblank(metadata.ses_name),' ','');
        run_label     = strcat('run-',deblank(metadata.run_name),header.hour,header.min,' ','');
        
        % make directories
        sub_dir       = fullfile([proj_diroutput,proj_diroutputcopy],sub_label);
        ses_dir       = fullfile([proj_diroutput,proj_diroutputcopy],sub_label,ses_label);
        ieeg_dir      = fullfile([proj_diroutput,proj_diroutputcopy],sub_label,ses_label,'ieeg');
        ieeg_file     = strcat(sub_label,'_',ses_label,'_',task_label,'_',run_label);
        
        for i=1:size(sub_dir,2)
            mydirMaker(sub_dir{i});
            mydirMaker(ses_dir{i});
            mydirMaker(ieeg_dir{i});
        end
        
        %check if it is empty, otherwise remove tsv,json,eeg,vhdr,trc,vmrk
        for i=1:size(ieeg_dir,2)
            ieeg_files = dir(ieeg_dir{i});
            
            if contains([ieeg_files(:).name],ieeg_file)
                
                delete(fullfile(ieeg_dir{i},[ieeg_file '*.tsv']))  ; 
                delete(fullfile(ieeg_dir{i},[ieeg_file '*.json'])) ;
                delete(fullfile(ieeg_dir{i},[ieeg_file,'*.eeg']))  ;
                delete(fullfile(ieeg_dir{i},[ieeg_file,'*.vhdr'])) ;
                delete(fullfile(ieeg_dir{i},[ieeg_file,'*.vmrk'])) ;
                delete(fullfile(ieeg_dir{i},[ieeg_file,'*.TRC']))  ;
                
            end
        end
        
        fieeg_name = strcat(sub_label,'_',ses_label,'_',task_label,'_',run_label,'_','ieeg',exte);
        fieeg_json_name = strcat(sub_label,'_',ses_label,'_',task_label,'_',run_label,'_','ieeg','.json');
        fchs_name = strcat(sub_label,'_',ses_label,'_',task_label,'_',run_label,'_','channels','.tsv');
        fevents_name = strcat(sub_label,'_',ses_label,'_',task_label,'_',run_label,'_','events','.tsv');
        felec_name = strcat(sub_label,'_',ses_label,'_','electrodes','.tsv');
        fcoords_name = strcat(sub_label,'_',ses_label,'_','coordsystem','.json');
        fpic_name = strcat(sub_label,'_',ses_label,'_','photo','.jpg');
        
        % file ieeg of the recording to .vhdr extension
        fileTRC = cell(1);
        fileVHDR = cell(1);
        fileVHDRcopy = cell(1);
        
        for i=1:size(ieeg_dir,2)
            fileTRC{i}  = fullfile(ieeg_dir{i},fieeg_name);
            if i==1
                fileVHDR{i} = replace(fileTRC{i},'.TRC','.vhdr');
            elseif i>1
                fileVHDRcopy{i-1} = replace(fileTRC{i},'.TRC','.vhdr');
            end
        end
        
        %% create Brainvision format from TRC
        
        cfg = [];
        cfg.dataset                     = filename; 
        cfg.continuous = 'yes';
        data2write = ft_preprocessing(cfg);
        
        cfg = [];
        cfg.outputfile                  = fileVHDR{1};
        
        cfg.mri.writesidecar       = 'no';
        cfg.meg.writesidecar        = 'no';
        cfg.eeg.writesidecar        = 'no';
        cfg.ieeg.writesidecar       = 'no';
        cfg.channels.writesidecar   = 'no';
        cfg.events.writesidecar     = 'no';
        
        % write .vhdr, .eeg, .vmrk
        data2bids(cfg, data2write)
                
        % to each output-file in fileVHDRcopy, the data should be copied
        if ~isempty(fileVHDRcopy{1})
            for i=1:size(fileVHDRcopy,2)
                copyfile(cfg.outputfile,fileVHDRcopy{i});
                fprintf('Copy to %s\n', fileVHDRcopy{i})
                file_eeginp = replace(cfg.outputfile,'.vhdr','.eeg');
                file_eegoutp = replace(fileVHDRcopy{i},'.vhdr','.eeg');
                copyfile(file_eeginp,file_eegoutp);
                fprintf('Copy to %s\n',file_eegoutp)
                file_vmrkinp = replace(cfg.outputfile,'.vhdr','.vmrk');
                file_vmrkoutp = replace(fileVHDRcopy{i},'.vhdr','.vmrk');
                copyfile(file_vmrkinp,file_vmrkoutp);
                fprintf('Copy to %s\n',file_vmrkoutp)
            end
        end
        
        %% create json sidecar for ieeg file
        cfg                             = [];
        cfg.ieeg                        = struct;
        cfg.channels                    = struct;
        cfg.electrodes                  = struct;
        cfg.coordsystem                 = struct;
        
        cfg.outputfile                  = fileVHDR{1};

        cfg.TaskName                    = task_label;
        cfg.TaskDescription             = task_desc;
        cfg.InstitutionName             = 'University Medical Center Utrecht';
        cfg.InstitutionalDepartmentName = 'Clinical Neurophysiology Department';
        cfg.InstitutionAddress          = 'Heidelberglaan 100, 3584 CX Utrecht';
        cfg.Manufacturer                = 'Micromed';
        cfg.ManufacturersModelName      = header.acquisition_eq;%sprintf('Acqui.eq:%i  File_type:%i',header.acquisition_eq,header.file_type);
        cfg.DeviceSerialNumber          = '';
        cfg.SoftwareVersions            = num2str(header.Header_Type);
        cfg.SoftwareFilters             = 'n/a';
        if strfind(cfg.ManufacturersModelName,'LTM') ~=0
            cfg.HardwareFilters.HighpassFilter.CutoffFrequency =             0.15;
            if header.Rate_Min/2.21 < 468
                cfg.HardwareFilters.LowpassFilter.CutoffFrequency = round(header.Rate_Min/2.21);
            else
                cfg.HardwareFilters.LowpassFilter.CutoffFrequency  =             468;
            end
        elseif strcmp(cfg.ManufacturersModelName,'SD128')
            cfg.HardwareFilters.HighpassFilter.CutoffFrequency =             0.15;
            cfg.HardwareFilters.LowpassFilter.CutoffFrequency  =             round(header.Rate_Min/3.81);
        elseif strcmp(cfg.ManufacturersModelName,'SD64')
            cfg.HardwareFilters.HighpassFilter.CutoffFrequency =             0.15;
            cfg.HardwareFilters.LowpassFilter.CutoffFrequency  =             round(header.Rate_Min/3.81);
        end
        
        %% create _channels.tsv
        % run json, channels.tsv
        json_sidecar_and_ch_and_ele_tsv(header,metadata,cfg);
        
        % to each output-file in fileVHDRcopy, the jsondata and channels.tsv should be copied
        if ~isempty(fileVHDRcopy{1})
            for i=1:size(fileVHDRcopy,2)
                cfg.outputfile = fileVHDRcopy{i};
                json_sidecar_and_ch_and_ele_tsv(header,metadata,cfg);
            end
        end
        
        %% create coordsystem.json
        cfg.outputfile                  = fileVHDR{1};
        cfg.coordsystem.iEEGCoordinateSystem                = []  ;
        cfg.coordsystem.iEEGCoordinateUnits                 = []      ;
        cfg.coordsystem.iEEGCoordinateProcessingDescription = []    ;
        cfg.coordsystem.IntendedFor                         =  fpic_name;
        
        json_coordsystem(cfg)
        
         % to each output-file in fileVHDRcopy, the json coordsystem should be copied
        if ~isempty(fileVHDRcopy{1})
            for i=1:size(fileVHDRcopy,2)
                cfg.outputfile = fileVHDRcopy{i};
                json_coordsystem(cfg);
            end
        end
        
        %% write annotations of the TRC
        cfg.outputfile                  = fileVHDR{1};
        annotations_tsv = write_annotations_tsv(header,metadata,annots,cfg);
        
        % to each output-file in fileVHDRcopy, the events.tsv should be copied
        if ~isempty(fileVHDRcopy{1})
            for i=1:size(fileVHDRcopy,2)
                file_eventsinp = replace(cfg.outputfile,'_ieeg.vhdr','_events.tsv');
                file_eventsoutp = replace(fileVHDRcopy{i},'_ieeg.vhdr','_events.tsv');
                copyfile(file_eventsinp,file_eventsoutp);
                fprintf('copying to %s \n',file_eventsoutp)
            end
        end
        
        %% write scans-file
        cfg.outputfile                  = fileVHDR{1};
        write_scans_tsv(cfg,metadata,annotations_tsv)
        
        % to each output-file in fileVHDRcopy, the scans.tsv should be adapted to the scans in that directory
        if ~isempty(fileVHDRcopy{1})
            for i=1:size(fileVHDRcopy,2)
                cfg.outputfile = fileVHDRcopy{i};
                write_scans_tsv(cfg,metadata,annotations_tsv)
            end
        end
                
        %% write participants-file
        cfg.outputfile                  = fileVHDR{1};
        write_participants_tsv(cfg,header)

        % to each output-file in fileVHDRcopy, the participants.tsv should be adapted to the scans in that directory
        if ~isempty(fileVHDRcopy{1})
            for i=1:size(fileVHDRcopy,2)
                cfg.outputfile = fileVHDRcopy{i};
                write_participants_tsv(cfg,header)
            end
        end
        
        %% write dataset descriptor
        create_datasetDesc(proj_diroutput{1})
        if ~isempty(proj_diroutputcopy{1})
            for i=1:size(proj_diroutputcopy,2)
                create_datasetDesc(proj_diroutputcopy{i})
            end
        end
        
        %% write event descriptor
        create_eventDesc(proj_diroutput{1})
        if ~isempty(proj_diroutputcopy{1})
            for i=1:size(proj_diroutputcopy,2)
                create_eventDesc(proj_diroutputcopy{i})
            end
        end
        
    else
        %% errors in parsing the data
        error(msg)
    end
    
catch ME
    
    status = 1;
    if (isempty(msg))
        msg = sprintf('%s err:%s --func:%s',filename,ME.message,ME.stack(1).name);
    else
        msg = sprintf('%s err:%s %s --func:%s',filename,msg,ME.message,ME.stack(1).name);
    end
    
end


%% create dataset descriptor
function create_datasetDesc(proj_dir)

ddesc_json.Name               = 'RESPect' ;
ddesc_json.BIDSVersion        = 'BEP010';
ddesc_json.License            = 'Not licenced yet';
ddesc_json.Authors            = {'van Blooijs D., Demuru M., Zweiphenning W.J.E.M, Leijten F.S.S., Zijlmans M.'};
ddesc_json.Acknowledgements   = 'Huiskamp G.J.M.';
ddesc_json.HowToAcknowledge   = 'possible paper to quote' ;
ddesc_json.Funding            = 'Epi-Sign Project and Epilepsiefonds #17-07' ;
ddesc_json.ReferencesAndLinks = {'articles and/or links'};
ddesc_json.DatasetDOI         = 'DOI of the dataset if online';


if ~isempty(ddesc_json)
    
    filename = fullfile(proj_dir,'dataset_description.json');
    if isfile(filename)
        existing = read_json(filename);
    else
        existing = [];
    end
    write_json(filename, mergeconfig(existing, ddesc_json))
    %     json_options.indent = ' ';
    %     jsonwrite(filename, mergeconfig(existing, ddesc_json), json_options)
end

%% create dataset descriptor
function create_eventDesc(proj_dir)

edesc_json.onset                                = 'onset of event in seconds' ;
edesc_json.offset                               = 'offset of event in seconds' ;
edesc_json.duration                             = 'duration of event in seconds' ;
edesc_json.trial_type                           = 'type of event (electrical stimulation/motor task/sensing task/artefact/sleep/sleep wake transition/eyes open)' ;
edesc_json.subtype                              = 'more description of event (sleep:nrem/rem, motor:Mario/hand/jump, sens:circle, electrical stimulation:SPES/ESM)' ;
edesc_json.electrodes_involved_onset            = 'electrodes involved in onset. For example: in seizure: electrodes involved in seizure onset or in artefact.' ;
edesc_json.electrodes_involved_offset           = 'electrodes involved in offset. For example: in seizure: electrodes involved in the end of a seizure or in an artefact.' ;
edesc_json.sample_start                         = 'onset of event in samples' ;
edesc_json.sample_end                           = 'offset of event in samples' ;
edesc_json.electrical_stimulation_type          = 'type of electrical stimulation [mono-/biphasic]';
edesc_json.electrical_stimulation_site          = 'electrode names of stimulus pair';
edesc_json.electrical_stimulation_site_num_1    = 'electrode one in stimulus pair' ;
edesc_json.electrical_stimulation_site_num_2    = 'electrode two in stimulus pair' ;
edesc_json.electrical_stimulation_current       = 'electrical stimulation current in A';
edesc_json.notes                                = 'notes about stimulation current';

if ~isempty(edesc_json)
    
    filename = fullfile(proj_dir,'event_description.json');
    if isfile(filename)
        existing = read_json(filename);
    else
        existing = [];
    end
    write_json(filename, mergeconfig(existing, edesc_json))
    %     json_options.indent = ' ';
    %     jsonwrite(filename, mergeconfig(existing, edesc_json), json_options)
end

%% function for json and tsv ieeg following fieldtrip style
function json_sidecar_and_ch_and_ele_tsv(header,metadata,cfg)
ch2use_included = metadata.ch2use_included;


%% Generic fields for all data types
cfg.TaskName                          = ft_getopt(cfg, 'TaskName'                    ); % REQUIRED. Name of the task (for resting state use the “rest” prefix). Different Tasks SHOULD NOT have the same name. The Task label is derived from this field by removing all non alphanumeric ([a-zA-Z0-9]) characters.
cfg.TaskDescription                   = ft_getopt(cfg, 'TaskDescription',''          ); % OPTIONAL. Description of the task.
cfg.Manufacturer                      = ft_getopt(cfg, 'Manufacturer'                ); % OPTIONAL. Manufacturer of the MEG system ("CTF", "​Elekta/Neuromag​", "​4D/BTi​", "​KIT/Yokogawa​", "​ITAB​", "KRISS", "Other")
cfg.ManufacturersModelName            = ft_getopt(cfg, 'ManufacturersModelName'      ); % OPTIONAL. Manufacturer’s designation of the MEG scanner model (e.g. "CTF-275"). See ​Appendix VII​ with preferred names
cfg.DeviceSerialNumber                = ft_getopt(cfg, 'DeviceSerialNumber',''       ); % OPTIONAL. The serial number of the equipment that produced the composite instances. A pseudonym can also be used to prevent the equipment from being identifiable, as long as each pseudonym is unique within the dataset.
cfg.SoftwareVersions                  = ft_getopt(cfg, 'SoftwareVersions',''         ); % OPTIONAL. Manufacturer’s designation of the acquisition software.
cfg.InstitutionName                   = ft_getopt(cfg, 'InstitutionName'             ); % OPTIONAL. The name of the institution in charge of the equipment that produced the composite instances.
cfg.InstitutionAddress                = ft_getopt(cfg, 'InstitutionAddress'          ); % OPTIONAL. The address of the institution in charge of the equipment that produced the composite instances.
cfg.InstitutionalDepartmentName       = ft_getopt(cfg, 'InstitutionalDepartmentName' ); % The department in the institution in charge of the equipment that produced the composite instances. Corresponds to DICOM Tag 0008, 1040 ”Institutional Department Name”.
cfg.SoftwareFilters                   = ft_getopt(cfg, 'SoftwareFilters');

%% IEEG inherited fields used

cfg.ieeg.ECOGChannelCount             = ft_getopt(cfg.ieeg, 'ECOGChannelCount'  ); %RECOMMENDED
cfg.ieeg.SEEGChannelCount             = ft_getopt(cfg.ieeg, 'SEEGChannelCount'  ); %RECOMMENDED
cfg.ieeg.EEGChannelCount              = ft_getopt(cfg.ieeg, 'EEGChannelCount'   ); %RECOMMENDED
cfg.ieeg.EOGChannelCount              = ft_getopt(cfg.ieeg, 'EOGChannelCount'   ); %RECOMMENDED
cfg.ieeg.ECGChannelCount              = ft_getopt(cfg.ieeg, 'ECGChannelCount'   ); %RECOMMENDED
cfg.ieeg.EMGChannelCount              = ft_getopt(cfg.ieeg, 'EMGChannelCount'   ); %RECOMMENDED
cfg.ieeg.RecordingDuration            = ft_getopt(cfg.ieeg, 'RecordingDuration' ); %RECOMMENDED
cfg.ieeg.RecordingType                = ft_getopt(cfg.ieeg, 'RecordingType'     ); %RECOMMENDED
cfg.ieeg.EpochLength                  = ft_getopt(cfg.ieeg, 'EpochLength'       ); %RECOMMENDED


%% IEEG specific fields
cfg.ieeg.SamplingFrequency            = ft_getopt(cfg.ieeg, 'SamplingFrequency'          ); % REQUIRED.
cfg.ieeg.PowerLineFrequency           = ft_getopt(cfg.ieeg, 'PowerLineFrequency'         ); % REQUIRED.
cfg.ieeg.iEEGReference                = ft_getopt(cfg.ieeg, 'iEEGReference'              ); % REQUIRED.
cfg.ieeg.ElectrodeManufacturer        = ft_getopt(cfg.ieeg, 'ElectrodeManufacturer'      ); %RECOMMENDED
cfg.ieeg.iEEGElectrodeGroups          = ft_getopt(cfg.ieeg, 'iEEGElectrodeGroups'        ); %RECOMMENDED


ft_warning('iEEG metadata fields need to be updated with the draft specification at http://bit.ly/bids_ieeg');


%% columns in the channels.tsv
cfg.channels.name               = ft_getopt(cfg.channels, 'name'               , nan);  % REQUIRED. Channel name (e.g., MRT012, MEG023)
cfg.channels.type               = ft_getopt(cfg.channels, 'type'               , nan);  % REQUIRED. Type of channel; MUST use the channel types listed below.
cfg.channels.units              = ft_getopt(cfg.channels, 'units'              , nan);  % REQUIRED. Physical unit of the data values recorded by this channel in SI (see Appendix V: Units for allowed symbols).
cfg.channels.description        = ft_getopt(cfg.channels, 'description'        , nan);  % OPTIONAL. Brief free-text description of the channel, or other information of interest. See examples below.
cfg.channels.sampling_frequency = ft_getopt(cfg.channels, 'sampling_frequency' , nan);  % OPTIONAL. Sampling rate of the channel in Hz.
cfg.channels.low_cutoff         = ft_getopt(cfg.channels, 'low_cutoff'         , nan);  % OPTIONAL. Frequencies used for the high-pass filter applied to the channel in Hz. If no high-pass filter applied, use n/a.
cfg.channels.high_cutoff        = ft_getopt(cfg.channels, 'high_cutoff'        , nan);  % OPTIONAL. Frequencies used for the low-pass filter applied to the channel in Hz. If no low-pass filter applied, use n/a. Note that hardware anti-aliasing in A/D conversion of all MEG/EEG electronics applies a low-pass filter; specify its frequency here if applicable.
cfg.channels.reference          = ft_getopt(cfg.channels, 'reference'          , nan);  % OPTIONAL.
cfg.channels.group              = ft_getopt(cfg.channels, 'group'              , nan);  % OPTIONAL.
cfg.channels.notch              = ft_getopt(cfg.channels, 'notch'              , nan);  % OPTIONAL. Frequencies used for the notch filter applied to the channel, in Hz. If no notch filter applied, use n/a.
cfg.channels.software_filters   = ft_getopt(cfg.channels, 'software_filters'   , nan);  % OPTIONAL. List of temporal and/or spatial software filters applied (e.g. "SSS", "SpatialCompensation"). Note that parameters should be defined in the general MEG sidecar .json file. Indicate n/a in the absence of software filters applied.
cfg.channels.status             = ft_getopt(cfg.channels, 'status'             , nan);  % OPTIONAL. Data quality observed on the channel (good/bad). A channel is considered bad if its data quality is compromised by excessive noise. Description of noise type SHOULD be provided in [status_description].
cfg.channels.status_description = ft_getopt(cfg.channels, 'status_description' , nan);  % OPTIONAL. Freeform text description of noise or artifact affecting data quality on the channel. It is meant to explain why the channel was declared bad in [status].


%% *_electrodes.tsv
cfg.electrodes.name             = ft_getopt(cfg.electrodes, 'name'               , nan);
cfg.electrodes.x                = ft_getopt(cfg.electrodes, 'x'                  , nan);
cfg.electrodes.y                = ft_getopt(cfg.electrodes, 'y'                  , nan);
cfg.electrodes.z                = ft_getopt(cfg.electrodes, 'z'                  , nan);
cfg.electrodes.size             = ft_getopt(cfg.electrodes, 'size'               , nan);
cfg.electrodes.group            = ft_getopt(cfg.electrodes, 'group'              , nan);
cfg.electrodes.material         = ft_getopt(cfg.electrodes, 'material'           , nan);
cfg.electrodes.manufacturer     = ft_getopt(cfg.electrodes, 'manufacturer'       , nan);
cfg.electrodes.silicon          = ft_getopt(cfg.electrodes, 'silicon'            , nan);
cfg.electrodes.soz              = ft_getopt(cfg.electrodes, 'soz'                , nan);
cfg.electrodes.ra               = ft_getopt(cfg.electrodes, 'ra'                 , nan);
cfg.electrodes.edge             = ft_getopt(cfg.electrodes, 'edge'               , nan);


%% start with empty  descriptions
ieeg_json    = [];
channels_tsv = [];

ieeg_json.TaskName                          = cfg.TaskName;
ieeg_json.TaskDescription                   = cfg.TaskDescription;
ieeg_json.Manufacturer                      = cfg.Manufacturer;
ieeg_json.ManufacturersModelName            = cfg.ManufacturersModelName;
ieeg_json.DeviceSerialNumber                = cfg.DeviceSerialNumber;
ieeg_json.SoftwareVersions                  = cfg.SoftwareVersions;
ieeg_json.SoftwareFilters                   = cfg.SoftwareFilters;
ieeg_json.InstitutionName                   = cfg.InstitutionName;
ieeg_json.InstitutionAddress                = cfg.InstitutionAddress;
ieeg_json.InstitutionalDepartmentName       = cfg.InstitutionalDepartmentName;
ieeg_json.HardwareFilters                   = cfg.HardwareFilters;
ch_label                                    = metadata.ch_label;

%% IEEG inherited fields used
if strcmpi(metadata.elec_info,'ecog')
    ieeg_json.ECOGChannelCount             = sum(ch2use_included);
    ieeg_json.SEEGChannelCount             = 0;
elseif strcmpi(metadata.elec_info,'seeg')
    ieeg_json.ECOGChannelCount             = 0;
    ieeg_json.SEEGChannelCount             = sum(ch2use_included);
end

%ieeg_json.EEGChannelCount              =
%ieeg_json.EOGChannelCount              =
ieeg_json.ECGChannelCount              = sum(~cellfun(@isempty,regexpi(ch_label,'ECG')));
%ieeg_json.EMGChannelCount              =
ieeg_json.RecordingDuration            = header.Num_Samples/header.Rate_Min;
ieeg_json.RecordingType                = 'continuous';
ieeg_json.EpochLength                  = 0;


%% IEEG specific fields
ieeg_json.SamplingFrequency            = header.Rate_Min;
ieeg_json.PowerLineFrequency           = 50;
ieeg_json.iEEGReference                = 'probably mastoid';
ieeg_json.ElectrodeManufacturer        = 'AdTech';
ieeg_json.iEEGElectrodeGroups          = metadata.format_info;
if strfind(cfg.TaskName,'SPES') ~=0
    ieeg_json.ElectricalStimulation        = 'true';
end

fn = {'name' 'type' 'units' 'low_cutoff' 'high_cutoff' 'reference' 'group' 'sampling_frequency'...
    'description' 'notch' 'status' 'status_description'};
for i=1:numel(fn)
    if numel(cfg.channels.(fn{i}))==1
        cfg.channels.(fn{i}) = repmat(cfg.channels.(fn{i}), header.Num_Chan, 1);
    end
end




%% iEEG  channels.tsv file
name                                = mergevector({header.elec(:).Name}', cfg.channels.name)                                   ;

type                                = cell(size(name))                                                                         ;
if(any(metadata.ch2use_included))
    if strcmpi(metadata.elec_info,'ECoG')
        [type{metadata.ch2use_included}]    = deal('ECOG');
    elseif strcmpi(metadata.elec_info,'SEEG')
        [type{metadata.ch2use_included}]    = deal('SEEG');
    end
end

if(any(~metadata.ch2use_included))
    [type{~metadata.ch2use_included}]   = deal('OTHER');
end
idx_ecg                             = ~cellfun(@isempty,regexpi(ch_label,'ECG'))                                               ;
idx_ecg                             = idx_ecg'                                                                                 ;

if(any(idx_ecg))
    [type{idx_ecg}]                     = deal('ECG')                                                                              ;
end

units                               = mergevector({header.elec(:).Unit}', cfg.channels.units)                                  ;
% units                               = replace(units(:),'u',char(181))                                                          ;
sampling_frequency                  = mergevector(repmat(header.Rate_Min, header.Num_Chan, 1), cfg.channels.sampling_frequency);
low_cutoff                          = cell(size(name)) ;
[low_cutoff{:}]                     = deal(cfg.HardwareFilters.LowpassFilter.CutoffFrequency);%{header.elec(:).Prefiltering_LowPass_Limit}'                                             ;
high_cutoff                         = cell(size(name)) ;
[high_cutoff{:}]                    = deal(cfg.HardwareFilters.HighpassFilter.CutoffFrequency);%{header.elec(:).Prefiltering_HiPass_Limit}'                                              ;
reference                           = {header.elec(:).Ref}'                                                                    ;

group                               = extract_group_info(metadata)                                                             ;

notch                               = repmat('n/a',header.Num_Chan, 1)                                                         ;
% software_filters                    = repmat('n/a',header.Num_Chan, 1)                                                         ;

[ch_status,ch_status_desc]          = status_and_description(metadata)                                                         ;
status                              = ch_status                                                                                ;
status_description                  = ch_status_desc                                                                           ;

channels_tsv                        = table(name, type, units,  low_cutoff,    ...
    high_cutoff, reference, group, sampling_frequency,   ...
    notch, status, status_description                                                  );

%% electrode table
fn = {'name' 'x' 'y' 'z' 'size' 'group' 'material' 'manufacturer' 'silicon' 'soz' 'ra' 'edge'};
for i=1:numel(fn)
    if numel(cfg.electrodes.(fn{i}))==1
        cfg.electrodes.(fn{i}) = repmat(cfg.electrodes.(fn{i}), header.Num_Chan, 1);
    end
end

%name                                = mergevector({header.elec(:).Name}', cfg.electrodes.name)                                   ;
x                                         = repmat({0},header.Num_Chan,1)                                                            ;
y                                         = repmat({0},header.Num_Chan,1)                                                            ;
z                                         = repmat({0},header.Num_Chan,1)                                                            ;
e_size                                    = repmat({'n/a'},header.Num_Chan,1)                                                          ; %TODO ask
material                                  = repmat({'n/a'},header.Num_Chan,1)                                                          ; %TODO ask
manufacturer                              = repmat({'n/a'},header.Num_Chan,1)                                                          ; %TODO ask
silicon                                   = repmat({'no'},header.Num_Chan,1)                                                          ; %TODO ask
soz                                       = repmat({'no'},header.Num_Chan,1)                                                          ; %TODO ask
resected                                  = repmat({'no'},header.Num_Chan,1)                                                          ; %TODO ask
edge                                      = repmat({'no'},header.Num_Chan,1)                                                          ; %TODO ask
if strcmpi(metadata.elec_info,'SEEG')
    screw                                     = repmat({'no'},header.Num_Chan,1)                                                          ; %TODO ask
    csf                                       = repmat({'no'},header.Num_Chan,1)                                                          ; %TODO ask
    lesion                                    = repmat({'no'},header.Num_Chan,1)                                                          ; %TODO ask
    gm                                        = repmat({'no'},header.Num_Chan,1)                                                          ; %TODO ask
    wm                                        = repmat({'no'},header.Num_Chan,1)                                                          ; %TODO ask
    hipp                                      = repmat({'no'},header.Num_Chan,1)                                                          ; %TODO ask
    amyg                                      = repmat({'no'},header.Num_Chan,1)                                                          ; %TODO ask
    gliosis                                   = repmat({'no'},header.Num_Chan,1)                                                          ; %TODO ask
end

if(any(metadata.ch2use_included))
    [e_size{metadata.ch2use_included}]        = deal('4.2')                                                                                ;
end

if(any(metadata.ch2use_included))
    [material{metadata.ch2use_included}]      = deal('Platinum')                                                                                ;
end

if(any(metadata.ch2use_included))
    [manufacturer{metadata.ch2use_included}]  = deal('AdTech')                                                                                ;
end

if(any(metadata.ch2use_silicon))
    [silicon{metadata.ch2use_silicon}]  = deal('yes')                                                                                ;
end

if(any(metadata.ch2use_soz))
    [soz{metadata.ch2use_soz}]  = deal('yes')                                                                                ;
end

if(any(metadata.ch2use_resected))
    [resected{metadata.ch2use_resected}]  = deal('yes')                                                                                ;
end

if(any(metadata.ch2use_edge))
    [edge{metadata.ch2use_edge}]  = deal('yes')                                                                                ;
end

if strcmpi(metadata.elec_info,'SEEG')
    if(any(metadata.ch2use_lesion))
        [lesion{metadata.ch2use_lesion}]  = deal('yes')                                                                                ;
    end
    if(any(metadata.ch2use_wm))
        [wm{metadata.ch2use_wm}]  = deal('yes')                                                                                ;
    end
    if(any(metadata.ch2use_gm))
        [gm{metadata.ch2use_gm}]  = deal('yes')                                                                                ;
    end
    if(any(metadata.ch2use_csf))
        [csf{metadata.ch2use_csf}]  = deal('yes')                                                                                ;
    end
     if(any(metadata.ch2use_screw))
        [screw{metadata.ch2use_screw}]  = deal('yes')                                                                                ;
    end
    if(any(metadata.ch2use_hipp))
        [hipp{metadata.ch2use_hipp}]  = deal('yes')                                                                                ;
    end
 
    if(any(metadata.ch2use_amyg))
        [amyg{metadata.ch2use_amyg}]  = deal('yes')                                                                                ;
    end
    if(any(metadata.ch2use_gliosis))
        [gliosis{metadata.ch2use_gliosis}]  = deal('yes')                                                                                ;
    end

    electrodes_tsv                            = table(name, x , y, z, e_size, ...
        group, material, manufacturer, silicon, soz, resected, edge ,...
        screw,csf,wm,gm,hipp,amyg,lesion,gliosis,...
    'VariableNames',{'name', 'x', 'y', 'z', 'size', ...
    'group', 'material', 'manufacturer', 'silicon' 'soz','resected','edge',...
    'screw','csf','whitematter','graymatter','hippocampus','amygdala','lesion','gliosis'})     ;
    
else
    electrodes_tsv                            = table(name, x , y, z, e_size, group, material, manufacturer, silicon, soz, resected, edge ,...
    'VariableNames',{'name', 'x', 'y', 'z', 'size', 'group', 'material', 'manufacturer', 'silicon' 'soz','resected','edge'})     ;
end


if ~isempty(ieeg_json)
    [p, f, x] = fileparts(cfg.outputfile);
    filename = fullfile(p, [f '.json']);
    if isfile(filename)
        existing = read_json(filename);
    else
        existing = [];
    end
    write_json(filename, mergeconfig(existing, ieeg_json))
    %     json_options.indent = ' ';
    %     jsonwrite(filename, mergeconfig(existing, ieeg_json), json_options)
end

if ~isempty(channels_tsv)
    [p, f, x] = fileparts(cfg.outputfile);
    g = strsplit(f,'_ieeg');
    filename = fullfile(p, [g{1} '_channels.tsv']);
    if isfile(filename)
        existing = read_tsv(filename);
    else
        existing = [];
    end % try
    if ~isempty(existing)
        ft_error('existing file is not empty');
    end
    write_tsv(filename, channels_tsv);
end

if metadata.incl_exist == 1 % format and included were annotated in the file, file should be saved in /Fridge/chronic_ECoG and in current dataset
    [p] = fileparts(cfg.outputfile);
    %g = strsplit(f,'_ieeg');
    q = strsplit(p,'/');
    qsub = contains(q,'sub');
    qses = contains(q,'ses');
    filename_dataset = fullfile(p, [q{qsub} '_' q{qses} '_electrodes.tsv']);
    
    if isfile(filename_dataset)
        
        cc_elec_old = readtable(filename_dataset,'FileType','text','Delimiter','\t');
        
        struct1 = table2struct(cc_elec_old);
        struct2 = table2struct(electrodes_tsv);
        if ~isequal(struct1,struct2)
            fprintf('%s exists!\n',filename_dataset)
            n=1;
            while isfile(filename_dataset)
                nameminelec = strsplit(filename_dataset,'electrodes');
                filename_dataset = [nameminelec{1} 'electrodes_' num2str(n) '.tsv'];
                n=n+1;
            end
        end
    end
   write_tsv(filename_dataset, electrodes_tsv);
end


%% write json coordsystem
function json_coordsystem(cfg)

cfg.coordsystem.iEEGCoordinateSystem                = ft_getopt(cfg.coordsystem, 'iEEGCoordinateSystem'               , nan);
cfg.coordsystem.iEEGCoordinateUnits                 = ft_getopt(cfg.coordsystem, 'iEEGCoordinateUnits'                , nan);
cfg.coordsystem.iEEGCoordinateProcessingDescription = ft_getopt(cfg.coordsystem, 'iEEGCoordinateProcessingDescription', nan);
cfg.coordsystem.IntendedFor                         = ft_getopt(cfg.coordsystem, 'IntendedFor'                         ,nan);

coordsystem_json=[];
coordsystem_json.iEEGCoordinateSystem                    = cfg.coordsystem.iEEGCoordinateSystem                                  ;
coordsystem_json.iEEGCoordinateUnits                     = cfg.coordsystem.iEEGCoordinateUnits                                   ;
coordsystem_json.iEEGCoordinateProcessingDescription     = cfg.coordsystem.iEEGCoordinateProcessingDescription                   ;
coordsystem_json.IntendedFor                             = cfg.coordsystem.IntendedFor                                           ;

if ~isempty(coordsystem_json)
    [p, f, x] = fileparts(cfg.outputfile);
    g = strsplit(f,'_task');
    filename = fullfile(p, [g{1} '_coordsystem.json']);
    filename = replace(filename,'_task-acute','');
    if isfile(filename)
        existing = read_json(filename);
    else
        existing = [];
    end
    write_json(filename, mergeconfig(existing, coordsystem_json))
    %     json_options.indent = ' ';
    %     jsonwrite(filename, mergeconfig(existing, coordsystem_json), json_options)
    
end


%% write annotations to a tsv file _annotations
function annotation_tsv = write_annotations_tsv(header,metadata,annots,cfg)

%% type / sample start / sample end /  chname;
ch_label  = metadata.ch_label;
ch2use_included = metadata.ch2use_included;
fs = header.Rate_Min;

eventsannots.type        = {};
eventsannots.sub_type    = {};
eventsannots.s_start     = {};
eventsannots.s_end       = {};
% ch_name     = {};
eventsannots.ch_name_on  = {};
eventsannots.ch_name_off = {};

annots_new = annots;

%% artefacts

[annots_new, eventsannots ] = add_event2annotation(metadata.artefacts,'artefact',eventsannots, annots_new,header);

% if(~isempty(artefact))
%     for i=1:numel(artefact)
%
%         type{cc}    = 'artefact'                           ;
%         sub_type{cc} = 'n/a' ;
%         ch_name_on{cc} = 'n/a';
%         ch_name_off{cc} = 'n/a';
%         s_start{cc} = round(artefact{i}.pos(1)/fs,1); % time in seconds (1 decimal)
%         samp_start{cc} = num2str(artefact{i}.pos(1))          ;
%         s_end{cc}   = round(artefact{i}.pos(end)/fs,1); % time in seconds (1 decimal)
%         samp_end{cc} = num2str(artefact{i}.pos(end))          ;
%
%         if(isempty(artefact{i}.ch_names))
%            name = 'all';%metadata.ch_label{metadata.ch2use_included} ;
%             % error('artefact channel name wrong')
%         else
%             if size(artefact{i}.ch_names,2) == 1
%                 name = artefact{i}.ch_names{1}              ;
%             else
%                 for j=1
%                     name =[artefact{i}.ch_names{1}];
%                 end
%
%                 for j=2:size(artefact{i}.ch_names,2)
%                     name = [name ,',', artefact{i}.ch_names{j}];
%                 end
%             end
%         end
%         ch_name_on{cc} = name;
%         ch_name_off{cc} = ch_name_on{cc};
%
%         annots_new([annots_new{:,1}]==artefact{i}.pos(1),:)=[];
%         annots_new([annots_new{:,1}]==artefact{i}.pos(end),:)=[];
%
%         duration{cc} = round(s_end{cc} - s_start{cc},3);
%         stim_type{cc} = 'n/a';
%         site_name{cc} = 'n/a';
%         site_channum{cc} = 'n/a';
%         stim_cur{cc} = 'n/a';
%         notes{cc} = 'n/a';
%
%         cc          = cc + 1                               ;
%
%     end
% end

%% seizures

[annots_new, eventsannots ] = add_event2annotation(metadata.seizure,'seizure',eventsannots, annots_new,header);


% seizure = metadata.seizure;
%
% if(~isempty(seizure))
%     for i=1:numel(seizure)
%
%         type{cc}    = 'seizure'                           ;
%         sub_type{cc} = seizure{i}.type;
%         s_start{cc} = round(seizure{i}.pos(1)/fs,1); % time in seconds (1 decimal)
%         samp_start{cc} = num2str(seizure{i}.pos(1))          ;
%         s_end{cc}   = round(seizure{i}.pos(end)/fs,1); % time in seconds (1 decimal)
%         samp_end{cc} = num2str(seizure{i}.pos(end))          ;
%        ch_name_on{cc} = 'n/a';
%        ch_name_off{cc} = 'n/a';
%
%         if size(seizure{i}.ch_names_on,2) ==1
%             name = seizure{i}.ch_names_on{1}              ;
%         else
%             for j=1
%                 name =[seizure{i}.ch_names_on{1}];
%             end
%
%             for j=2:size(seizure{i}.ch_names_on,2)
%                 name = [name ,',', seizure{i}.ch_names_on{j}];
%             end
%
%         end
%         ch_name_on{cc} = name;
%
%         if size(seizure{i}.ch_names_off,2) ==1
%             name = seizure{i}.ch_names_off{1}              ;
%         else
%             for j=1
%                 name =[seizure{i}.ch_names_off{1}];
%             end
%
%             for j=2:size(seizure{i}.ch_names_off,2)
%                 name = [name ,',', seizure{i}.ch_names_off{j}];
%             end
%         end
%         ch_name_off{cc} = name;
%
%         annots_new([annots_new{:,1}]==seizure{i}.pos(1),:)=[];
%         annots_new([annots_new{:,1}]==seizure{i}.pos(end),:)=[];
%         duration{cc} = round(s_end{cc} - s_start{cc},3);
%         stim_type{cc} = 'n/a';
%         site_name{cc} = 'n/a';
%         site_channum{cc} = 'n/a';
%         stim_cur{cc} = 'n/a';
%         notes{cc} = 'n/a';
%
%         cc          = cc + 1                               ;
%
%     end
% end

%% stimulation

[annots_new, eventsannots ] = add_event2annotation(metadata.stimulation,'stimulation',eventsannots, annots_new,header);

% stimulation = metadata.stimulation;
%
% if(~isempty(stimulation))
%     for i=1:numel(stimulation)
%
%         type{cc}    = 'stimulation'                           ;
%         sub_type{cc} = 'n/a' ;
%         s_start{cc} = round(stimulation{i}.pos(1)/fs,1); % time in seconds (1 decimal)
%         samp_start{cc} = num2str(stimulation{i}.pos(1))          ;
%         s_end{cc}   = round(stimulation{i}.pos(end)/fs,1); % time in seconds (1 decimal)
%         samp_end{cc} = num2str(stimulation{i}.pos(end))          ;
%
%         if(isempty(stimulation{i}.ch_names))
%             ch_name_on{cc} = 'all';%metadata.ch_label{metadata.ch2use_included} ;
%             ch_name_off{cc} = ch_name_on{cc}; % error('artefact channel name wrong')
%         else
%             ch_name_on{cc} = stimulation{i}.ch_names{1}              ;
%             ch_name_off{cc} = ch_name_on{cc};
%         end
%
%         annots_new([annots_new{:,1}]==stimulation{i}.pos(1),:)=[];
%         annots_new([annots_new{:,1}]==stimulation{i}.pos(end),:)=[];
%         duration{cc} = round(s_end{cc} - s_start{cc},3);
%         stim_type{cc} = 'n/a';
%         site_name{cc} = 'n/a';
%         site_channum{cc} = 'n/a';
%         stim_cur{cc} = 'n/a';
%         notes{cc} = 'n/a';
%
%         cc          = cc + 1                               ;
%
%     end
% end

%% eyes open

[annots_new, eventsannots ] = add_event2annotation(metadata.eyes_open,'eyes',eventsannots, annots_new,header);

% eyes_open = metadata.eyes_open;
%
% if(~isempty(eyes_open))
%     for i=1:numel(eyes_open)
%
%         type{cc}    = 'eyes open'                           ;
%         sub_type{cc} = 'n/a' ;
%         s_start{cc} = round(eyes_open{i}.pos(1)/fs,1); % time in seconds (1 decimal)
%         samp_start{cc} = num2str(eyes_open{i}.pos(1))          ;
%         s_end{cc}   = round(eyes_open{i}.pos(end)/fs,1); % time in seconds (1 decimal)
%         samp_end{cc} = num2str(eyes_open{i}.pos(end))          ;
%
%         if(isempty(eyes_open{i}.ch_names))
%             ch_name_on{cc} = 'all';%metadata.ch_label{metadata.ch2use_included} ;
%             error('artefact channel name wrong')
%             ch_name_off{cc} = ch_name_on{cc} ;
%         else
%             ch_name_on{cc} = eyes_open{i}.ch_names{1}              ;
%             ch_name_off{cc} = ch_name_on{cc} ;
%         end
%
%         annots_new([annots_new{:,1}]==eyes_open{i}.pos(1),:)=[];
%         annots_new([annots_new{:,1}]==eyes_open{i}.pos(end),:)=[];
%         duration{cc} = round(s_end{cc} - s_start{cc},3);
%         stim_type{cc} = 'n/a';
%         site_name{cc} = 'n/a';
%         site_channum{cc} = 'n/a';
%         stim_cur{cc} = 'n/a';
%         notes{cc} = 'n/a';
%
%         cc          = cc + 1                               ;
%
%     end
% end




%% sleep wake transition

[annots_new, eventsannots ] = add_event2annotation(metadata.slaw_trans,'sleep-wake transition',eventsannots, annots_new,header);

% slaw_trans = metadata.slaw_trans;
%
% if(~isempty(slaw_trans))
%     for i=1:numel(slaw_trans)
%
%         type{cc}    = 'sleep-wake transition'                           ;
%         sub_type{cc} = 'n/a' ;
%         s_start{cc} = round(slaw_trans{i}.pos(1)/fs,1); % time in seconds (1 decimal)
%         samp_start{cc} = num2str(slaw_trans{i}.pos(1))          ;
%         s_end{cc}   = round(slaw_trans{i}.pos(end)/fs,1); % time in seconds (1 decimal)
%         samp_end{cc} = num2str(slaw_trans{i}.pos(end))          ;
%
%         if(isempty(slaw_trans{i}.ch_names))
%             ch_name_on{cc} = 'all';%metadata.ch_label{metadata.ch2use_included} ;
%             error('artefact channel name wrong')
%             ch_name_off{cc} = ch_name_on{cc} ;
%         else
%             ch_name_on{cc} = slaw_trans{i}.ch_names{1}              ;
%             ch_name_off{cc} = ch_name_on{cc} ;
%         end
%
%         annots_new([annots_new{:,1}]==slaw_trans{i}.pos(1),:)=[];
%         annots_new([annots_new{:,1}]==slaw_trans{i}.pos(end),:)=[];
%         duration{cc} = round(s_end{cc} - s_start{cc},3);
%         stim_type{cc} = 'n/a';
%         site_name{cc} = 'n/a';
%         site_channum{cc} = 'n/a';
%         stim_cur{cc} = 'n/a';
%         notes{cc} = 'n/a';
%
%         cc          = cc + 1                               ;
%
%     end
% end


%% sleep

[annots_new, eventsannots ] = add_event2annotation(metadata.sleep,'sleep',eventsannots, annots_new,header);

% sleep = metadata.sleep;
%
% if(~isempty(sleep))
%     for i=1:numel(sleep)
%
%         type{cc}    = 'sleep'                           ;
%         s_start{cc} = round(sleep{i}.pos(1)/fs,1); % time in seconds (1 decimal)
%         samp_start{cc} = num2str(sleep{i}.pos(1))          ;
%         s_end{cc}   = round(sleep{i}.pos(end)/fs,1); % time in seconds (1 decimal)
%         samp_end{cc} = num2str(sleep{i}.pos(end))          ;
%
%         if(isempty(sleep{i}.ch_names))
%             ch_name_on{cc} = 'all';%metadata.ch_label{metadata.ch2use_included} ;
%             % error('artefact channel name wrong')
%             ch_name_off{cc} = ch_name_on{cc} ;
%         else
%             ch_name_on{cc} = slaw_trans{i}.ch_names{1}              ;
%             ch_name_off{cc} = ch_name_on{cc} ;
%         end
%
%         if (isempty(sleep{i}.type))
%             sub_type{cc} = 'unknown' ;
%         else
%             sub_type{cc} = sleep{i}.type ;
%         end
%
%         annots_new([annots_new{:,1}]==sleep{i}.pos(1),:)=[];
%         annots_new([annots_new{:,1}]==sleep{i}.pos(end),:)=[];
%         duration{cc} = round(s_end{cc} - s_start{cc},3);
%         stim_type{cc} = 'n/a';
%         site_name{cc} = 'n/a';
%         site_channum{cc} = 'n/a';
%         stim_cur{cc} = 'n/a';
%         notes{cc} = 'n/a';
%
%         cc          = cc + 1                               ;
%
%     end
% end

%% motortask

[annots_new, eventsannots ] = add_event2annotation(metadata.motortask,'motortask',eventsannots, annots_new,header);

% motortask = metadata.motortask;
%
% if(~isempty(motortask))
%     for i=1:numel(motortask)
%
%         type{cc}    = 'motortask'                           ;
%         sub_type{cc} = 'n/a' ;
%         s_start{cc} = round(motortask{i}.pos(1)/fs,1); % time in seconds (1 decimal)
%         samp_start{cc} = num2str(motortask{i}.pos(1))          ;
%         s_end{cc}   = round(motortask{i}.pos(end)/fs,1); % time in seconds (1 decimal)
%         samp_end{cc} = num2str(motortask{i}.pos(end))          ;
%
%         if(isempty(motortask{i}.ch_names))
%             ch_name_on{cc} = 'all';%metadata.ch_label{metadata.ch2use_included} ;
%             % error('artefact channel name wrong')
%             ch_name_off{cc} = ch_name_on{cc} ;
%         else
%             ch_name_on{cc} = motortask{i}.ch_names{1}              ;
%             ch_name_off{cc} = ch_name_on{cc} ;
%         end
%
%         if (isempty(motortask{i}.type))
%             sub_type{cc} = 'unknown' ;
%         else
%             sub_type{cc} = motortask{i}.type ;
%         end
%
%         annots_new([annots_new{:,1}]==motortask{i}.pos(1),:)=[];
%         annots_new([annots_new{:,1}]==motortask{i}.pos(end),:)=[];
%         duration{cc} = round(s_end{cc} - s_start{cc},3);
%         stim_type{cc} = 'n/a';
%         site_name{cc} = 'n/a';
%         site_channum{cc} = 'n/a';
%         stim_cur{cc} = 'n/a';
%         notes{cc} = 'n/a';
%
%
%         cc          = cc + 1                               ;
%
%     end
% end

%% language task

[annots_new, eventsannots ] = add_event2annotation(metadata.langtask,'languagetask',eventsannots, annots_new,header);

% langtask = metadata.langtask;
%
% if(~isempty(langtask))
%     for i=1:numel(langtask)
%
%         type{cc}    = 'languagetask'                           ;
%         sub_type{cc} = 'n/a' ;
%         s_start{cc} = round(langtask{i}.pos(1)/fs,1); % time in seconds (1 decimal)
%         samp_start{cc} = num2str(langtask{i}.pos(1))          ;
%         s_end{cc}   = round(langtask{i}.pos(end)/fs,1); % time in seconds (1 decimal)
%         samp_end{cc} = num2str(langtask{i}.pos(end))          ;
%
%         if(isempty(langtask{i}.ch_names))
%             ch_name_on{cc} = 'all';%metadata.ch_label{metadata.ch2use_included} ;
%             % error('artefact channel name wrong')
%             ch_name_off{cc} = ch_name_on{cc} ;
%         else
%             ch_name_on{cc} = langtask{i}.ch_names{1}              ;
%             ch_name_off{cc} = ch_name_on{cc} ;
%         end
%
%         if (isempty(langtask{i}.type))
%             sub_type{cc} = 'unknown' ;
%         else
%             sub_type{cc} = langtask{i}.type ;
%         end
%
%         annots_new([annots_new{:,1}]==langtask{i}.pos(1),:)=[];
%         annots_new([annots_new{:,1}]==langtask{i}.pos(end),:)=[];
%         duration{cc} = round(s_end{cc} - s_start{cc},3);
%         stim_type{cc} = 'n/a';
%         site_name{cc} = 'n/a';
%         site_channum{cc} = 'n/a';
%         stim_cur{cc} = 'n/a';
%         notes{cc} = 'n/a';
%
%
%         cc          = cc + 1                               ;
%
%     end
% end

%% sensing task

[annots_new, eventsannots ] = add_event2annotation(metadata.senstask,'sensingtask',eventsannots, annots_new,header);

% senstask = metadata.senstask;
%
% if(~isempty(senstask))
%     for i=1:numel(senstask)
%
%         type{cc}    = 'sensing task'                           ;
%         sub_type{cc} = 'n/a' ;
%         s_start{cc} = round(senstask{i}.pos(1)/fs,1); % time in seconds (1 decimal)
%         samp_start{cc} = num2str(senstask{i}.pos(1))          ;
%         s_end{cc}   = round(senstask{i}.pos(end)/fs,1); % time in seconds (1 decimal)
%         samp_end{cc} = num2str(senstask{i}.pos(end))          ;
%
%         if(isempty(senstask{i}.ch_names))
%             ch_name_on{cc} = 'all';%metadata.ch_label{metadata.ch2use_included} ;
%             % error('artefact channel name wrong')
%             ch_name_off{cc} = ch_name_on{cc} ;
%         else
%             ch_name_on{cc} = senstask{i}.ch_names{1}              ;
%             ch_name_off{cc} = ch_name_on{cc} ;
%         end
%
%         if (isempty(senstask{i}.type))
%             sub_type{cc} = 'unknown' ;
%         else
%             sub_type{cc} = senstask{i}.type ;
%         end
%
%         annots_new([annots_new{:,1}]==senstask{i}.pos(1),:)=[];
%         annots_new([annots_new{:,1}]==senstask{i}.pos(end),:)=[];
%         duration{cc} = round(s_end{cc} - s_start{cc},3);
%         stim_type{cc} = 'n/a';
%         site_name{cc} = 'n/a';
%         site_channum{cc} = 'n/a';
%         stim_cur{cc} = 'n/a';
%         notes{cc} = 'n/a';
%
%
%         cc          = cc + 1                               ;
%
%     end
% end

%% bsuppression
% bsuppression = metadata.bsuppression;
%
% if(~isempty(bsuppression))
%     for i=1:numel(bsuppression)
%
%         type{cc}    = 'bsuppression'                       ;
%         sub_type{cc} = 'n/a' ;
%         s_start{cc} = num2str(bsuppression{i}.pos(1))      ;
%         s_end{cc}   = num2str(bsuppression{i}.pos(end))    ;
%         ch_name_on{cc} = 'all'                                ;
%         ch_name_off{cc} = 'all' ;
%
%         cc          = cc + 1                               ;
%
%     end
% end

%% addnotes
% addnotes = metadata.add_notes;
%
% if(~isempty(addnotes))
%     for i=1:numel(addnotes)
%
%         type{cc}    = 'oddbehaviour'                       ;
%         sub_type{cc} = 'n/a' ;
%         s_start{cc} = num2str(addnotes{i}.pos(1))          ;
%         s_end{cc}   = num2str(addnotes{i}.pos(end))        ;
%
%         if(isempty(addnotes{i}.ch_names))
%             error('artefact channel name wrong')
%         end
%
%         ch_name_on{cc} = num2str(addnotes{i}.ch_names{1})     ;
%         ch_name_off{cc} = 'n/a' ;
%         cc          = cc + 1                               ;
%
%     end
% end

%% resected channels --> not in events but in electrodes.tsv
%
% resected = metadata.ch2use_resected;
%
% if(sum(resected))
%     idx_res  = find(resected);
%
%     for i=1:numel(idx_res)
%
%         type{cc}    = 'resected'                            ;
%         s_start{cc} = '1'                                   ;
%         s_end{cc}   = 'Inf'                                 ;
%
%         ch_name_on{cc} = ch_label(idx_res(i))                  ;
%         ch_name_off{cc} = 'n/a' ;
%         cc          = cc + 1                                ;
%
%     end
% end


%% resected channels --> not in events but in electrodes.tsv
%
% edge = metadata.ch2use_edge;
%
% if(sum(edge))
%     idx_edge  = find(edge);
%
%     for i=1:numel(idx_edge)
%
%         type{cc}    = 'edge'                                ;
%         s_start{cc} = '1'                                   ;
%         s_end{cc}   = 'Inf'                                 ;
%         ch_name_on{cc} = ch_label(idx_edge(i))                 ;
%         ch_name_off{cc} = 'n/a' ;
%         cc          = cc + 1                                ;
%
%     end
% end

%% triggers for good epochs -- in cECoG we do not annotate good epochs
% trigger = metadata.trigger;
% BEG_GS  = 222             ;
% END_GS  = 223             ;
%
% if(~isempty(trigger))
%     idx_begins  = find(trigger.val==BEG_GS);
%     idx_ends    = find(trigger.val==END_GS);
%
%     for i=1:numel(idx_begins)
%
%         type{cc}    = 'trial'                               ;
%         s_start{cc} = trigger.pos(idx_begins(i))            ;
%         s_end{cc}   = trigger.pos(idx_ends(i))              ;
%         ch_name{cc} = 'ALL'                                 ;
%         cc          = cc + 1                                ;
%
%     end
% end

%% electrical stimulation mapping
[annots_new, eventsannots ] = add_event2annotation(metadata.esm,'esm',eventsannots, annots_new,header);

% esm = metadata.esm;
%
% if(~isempty(esm))
%     for i=1:numel(esm)
%
%         type{cc}    = 'esm'                           ;
%         sub_type{cc} = 'n/a' ;
%         s_start{cc} = round(esm{i}.pos(1)/fs,1); % time in seconds (1 decimal)
%         samp_start{cc} = num2str(esm{i}.pos(1))          ;
%         s_end{cc}   = round(esm{i}.pos(end)/fs,1); % time in seconds (1 decimal)
%         samp_end{cc} = num2str(esm{i}.pos(end))          ;
%
%         if(isempty(esm{i}.ch_names))
%             ch_name_on{cc} = 'all';%metadata.ch_label{metadata.ch2use_included} ;
%             % error('artefact channel name wrong')
%             ch_name_off{cc} = ch_name_on{cc} ;
%         else
%             ch_name_on{cc} = langtask{i}.ch_names{1}              ;
%             ch_name_off{cc} = ch_name_on{cc} ;
%         end
%
%         if (isempty(esm{i}.type))
%             sub_type{cc} = 'n/a' ;
%         else
%             sub_type{cc} = esm{i}.type ;
%         end
%
%         annots_new([annots_new{:,1}]==esm{i}.pos(1),:)=[];
%         annots_new([annots_new{:,1}]==esm{i}.pos(end),:)=[];
%         duration{cc} = round(s_end{cc} - s_start{cc},3);
%         stim_type{cc} = 'n/a';
%         site_name{cc} = 'n/a';
%         site_channum{cc} = 'n/a';
%         stim_cur{cc} = 'n/a';
%         notes{cc} = 'n/a';
%
%         cc          = cc + 1                               ;
%
%     end
% end

%% visual selection data segments - slow wave sleep (or ar least late NREM stage) - aim: 10 consecutive minutes
[annots_new, eventsannots ] = add_event2annotation(metadata.SWSselection,'sws selection',eventsannots, annots_new,header);

%% visual selection data segments - rapid eye-movement sleep - aim: 10 consecutive minutes
[annots_new, eventsannots ] = add_event2annotation(metadata.REMselection,'rem selection',eventsannots, annots_new,header);

%% visual selection data segments - interictal awake - aim: 10 consecutive minutes
[annots_new, eventsannots ] = add_event2annotation(metadata.IIAWselection,'iiaw selection',eventsannots, annots_new,header);

%% visual selection data segments - Epileptogenicity index - about 5s before visual SOZ to 10s after, but sometimes different due to big artefacts/other signal characteristics
[annots_new, eventsannots ] = add_event2annotation(metadata.EIselection,'EI selection',eventsannots, annots_new,header);

%% adding trigger data to events list
% skip following section if no spes/esm/stimulation has been annotated in
% the file
if ~isempty(metadata.spes) || ~isempty(metadata.esm) || ~isempty(metadata.stimulation)
    trigger = metadata.trigger;
    % determine period of esm
    if any(cellfun(@(x) contains(x,{'esm'}),lower(eventsannots.sub_type)) & cellfun(@(x) ~contains(x,{'slow'}),lower(eventsannots.sub_type)))
        esmfind =    cellfun(@(x) contains(x,{'esm'}),lower(eventsannots.sub_type))==1;
        esmstart = str2double(eventsannots.samp_start{esmfind});
        esmend = str2double(eventsannots.samp_end{esmfind});
        
    elseif any(cellfun(@(x) contains(x,{'esm'}),lower(eventsannots.type)) & cellfun(@(x) ~contains(x,{'slow'}),lower(eventsannots.type)))
        esmfind =    cellfun(@(x) contains(x,{'esm'}),lower(eventsannots.type))==1;
        esmstart = str2double(eventsannots.samp_start{esmfind});
        esmend = str2double(eventsannots.samp_end{esmfind});
    else
        esmstart = 0;
        esmend = 0;
    end
    
    % determine period of spes
    if any(cellfun(@(x) contains(x,{'spes'}),lower(eventsannots.sub_type)))
        spesfind =    cellfun(@(x) contains(x,{'spes'}),lower(eventsannots.sub_type))==1;
        spesstart = str2double(eventsannots.samp_start{spesfind});
        spesend = str2double(eventsannots.samp_end{spesfind});
        
    elseif any(cellfun(@(x) contains(x,{'spes'}),lower(eventsannots.type)))
        spesfind =    cellfun(@(x) contains(x,{'spes'}),lower(eventsannots.type))==1;
        spesstart = str2double(eventsannots.samp_start{spesfind});
        spesend = str2double(eventsannots.samp_end{spesfind});
    else
        spesstart = 0;
        spesend = 0;
    end
    
    % determine period of slowesm (1hz cortical stimulation)
    if any(cellfun(@(x) contains(x,{'esm'}),lower(eventsannots.sub_type)) & cellfun(@(x) contains(x,{'slow'}),lower(eventsannots.sub_type)))
        slowesmfind =    cellfun(@(x) contains(x,{'esm'}),lower(eventsannots.sub_type))==1;
        slowesmstart = str2double(eventsannots.samp_start{slowesmfind});
        slowesmend = str2double(eventsannots.samp_end{slowesmfind});
        
    elseif any(cellfun(@(x) contains(x,{'esm'}),lower(eventsannots.type)) & cellfun(@(x) contains(x,{'slow'}),lower(eventsannots.type)))
        slowesmfind =    cellfun(@(x) contains(x,{'esm'}),lower(eventsannots.type))==1;
        slowesmstart = str2double(eventsannots.samp_start{slowesmfind});
        slowesmend = str2double(eventsannots.samp_end{slowesmfind});
    else
        slowesmstart = 0;
        slowesmend = 0;
    end
    
    if ~isempty(trigger.pos)
        idx_start = find(trigger.val >1000); % with cortical stimulation, triggers are added automatically with a number >1000
        for i=1:numel(idx_start)
            
            if trigger.pos(idx_start(i)) > spesstart && trigger.pos(idx_start(i)) < spesend % if trigger belongs to spes period
                eventsannots = add_spestrigger2annotation(idx_start(i),'SPES',eventsannots,annots_new, header,metadata,trigger)  ;
            elseif trigger.pos(idx_start(i)) > esmstart && trigger.pos(idx_start(i)) < esmend % if trigger belongs to esm period
                eventsannots = add_esmtrigger2annotation(idx_start(i),eventsannots,annots_new, header,metadata,trigger);
            elseif trigger.pos(idx_start(i)) > slowesmstart && trigger.pos(idx_start(i)) < slowesmend % if trigger belongs to slowesm period
                eventsannots = add_spestrigger2annotation(idx_start(i),'slowESM',eventsannots,annots_new, header,metadata,trigger)  ;
            else
                error('Trigger does not belong to SPES period or ESM period')
            end
            
        end
    end
    % in older ECoGs, there are no triggers, but stimulation NEED TO BE FIXED
    if ~isempty(sum(cellfun(@(x) contains(x,{'No trigger'}),annots_new(:,2)))) 
        if strcmpi(metadata.elec_info,'SEEG')
            stimcurdefault = 2;
        elseif strcmpi(metadata.elec_info,'ECoG')
            stimcurdefault = 8;
        end
        
        if contains(lower(metadata.stimcurr),'unknown')
            note = sprintf('Stimulation intensity is suggested to be %i mA but may differ when applied in eloquent tissue, triggers were added automatically in Matlab',stimcurdefault);
        else
            note = 'Triggers were added automatically in Matlab';
        end
        
        % load ECoG
        dataName = cfg.outputfile;
        data_raw = ft_read_data(dataName,'dataformat','brainvision_eeg');
        data = data_raw(ch2use_included,:);
        ch_label = {ch_label{ch2use_included}};
        %data(metadata.ch2use_bad,:) = NaN;
        
        % annotation sample
        numnotrigger = find(cellfun(@(x) contains(x,{'No trigger'}),annots_new(:,2))==1);
        for i=1:size(numnotrigger,1) % for each 'No trigger'-annotation
            if annots_new{numnotrigger(i),1} > spesstart && annots_new{numnotrigger(i),1} < spesend % if 'No trigger' is within SPESperiod
                if ~isempty(trigger.pos) % if triggers are present, then period ends with next trigger
                    periodend = trigger.pos(find(annots_new{numnotrigger(i),1}<trigger.pos,1,'first'))-round(0.25*fs); %-0.25s takes care of an annotation prior to the next trigger
                    sampend = periodend;
                else %otherwise, period ends with end of SPESperiod
                    periodend = spesend;
                    sampend = [];
                end
                % annotations within the specified period
                numannots = find([annots_new{:,1}]<periodend & [annots_new{:,1}] > annots_new{numnotrigger(i),1} ==1);
                
                for j=1:size(numannots,2)
                    % find stimulation pair and other stimulus settings
                    annotsplit = strsplit(annots_new{numannots(j),2},'_');
                    stimnumber = regexp(lower(annotsplit{1}),'\d*','match');
                    stimname = regexp(lower(annotsplit{1}),'[a-z]*','match');
                    %                    currannot = regexp(lower(annots_new{numannots(j),2}),'ma', 'once');
                    
                    if size(annotsplit,2)>1
                        currsplit = strsplit(lower(annotsplit{2}),'ma');
                        stimcurrstr = currsplit{1};
                        stimcurr = str2double(stimcurrstr)/1000;
                    else
                        stimcurr = stimcurdefault/1000;
                    end
                    
                    if size(stimnumber,2) == 2 && size(stimname,2)==2 % when both have size=2, then it should be a stimulus pair
                        
                        [~,stimnum] = findstimpair(stimnumber,stimname,ch_label);
                        
                        % find stimulation samples for each annotation
                        sampstart = annots_new{numannots(j),1};
                        if j==1 && isempty(sampend) % if there are no triggers present
                            sampend = annots_new{numannots(j)+1,1};
                        elseif numannots(j) == size(annots_new,2) % if triggers are not present in the last stimulation pair
                            sampend = spesend;
                        else % if there are triggers present, than sampend is before the next trigger
                            sampend = periodend;                    
                        end
                        
                        samplocs = findtrigger(data,stimnum,fs, sampstart, sampend)+annots_new{numannots(j),1};
                        
                        eventssize = size(eventsannots.type,2) ;
                        for cc = eventssize+1:eventssize+size(samplocs,2)
                            eventsannots.type{cc} = 'electrical_stimulation';
                            eventsannots.sub_type{cc} = 'SPES';
                            eventsannots.stim_type{cc} = 'monophasic';
                            eventsannots.samp_start{cc} = samplocs(cc-eventssize);
                            eventsannots.s_start{cc} = round(samplocs(cc-eventssize)/fs,1); % time in seconds (1 decimal)
                            
                            eventsannots.site_name{cc} = [ch_label{stimnum(1)}, '-', ch_label{stimnum(2)}];
                            eventsannots.site_channum{cc} = num2str([stimnum(1), stimnum(2)]);
                            eventsannots.duration{cc} = 1/1000;
                            eventsannots.s_end{cc} = 'n/a';
                            eventsannots.samp_end{cc} = 'n/a';
                            eventsannots.ch_name_on{cc} = 'n/a';
                            eventsannots.ch_name_off{cc} = 'n/a';
                            eventsannots.stim_cur{cc} = stimcurr;
                            eventsannots.notes{cc} = note;
                            
                        end
                        
                    end
                    
                    
                end
                
            elseif annots_new{numnotrigger(i),1} > esmstart && annots_new{numnotrigger(i),1} < esmend
                % FIX LATER
            elseif annots_new{numnotrigger(i),1} > slowesmstart && annots_new{numnotrigger(i),1} < slowesmend
                % FIX LATER
            end
            %         if ~isempty(trigger)
            %
            %         else
            %         end
        end
    end
    
    s_start         = eventsannots.s_start;
    s_end           = eventsannots.s_end;
    duration        = eventsannots.duration;
    type            = eventsannots.type;
    sub_type        = eventsannots.sub_type;
    ch_name_on      = eventsannots.ch_name_on;
    ch_name_off     = eventsannots.ch_name_off;
    samp_start      = eventsannots.samp_start;
    samp_end        = eventsannots.samp_end;
    stim_type       = eventsannots.stim_type;
    site_name       = eventsannots.site_name;
    site_channum    = eventsannots.site_channum;
    stim_cur        = eventsannots.stim_cur;
    notes           = eventsannots.notes;
    
    % sort to put the no-triggers in the right order of stimulation
    SPEStype        = strcmp(sub_type,'SPES');
    noSPEStype      = ~SPEStype;
    [~,I]           = sort([samp_start{SPEStype}]) ;
    I               = I + find(SPEStype==1,1,'first')-1;
    s_start         = {s_start{noSPEStype},s_start{I}};
    s_end           = {s_end{noSPEStype},s_end{I}};
    duration        = {duration{noSPEStype},duration{I}};
    type            = {type{noSPEStype},type{I}};
    sub_type        = {sub_type{noSPEStype},sub_type{I}};
    ch_name_on      = {ch_name_on{noSPEStype},ch_name_on{I}};
    ch_name_off     = {ch_name_off{noSPEStype},ch_name_off{I}};
    samp_start      = {samp_start{noSPEStype},samp_start{I}};
    samp_end        = {samp_end{noSPEStype},samp_end{I}};
    stim_type       = {stim_type{noSPEStype},stim_type{I}};
    site_name       = {site_name{noSPEStype},site_name{I}};
    site_channum    = {site_channum{noSPEStype},site_channum{I}};
    stim_cur        = {stim_cur{noSPEStype},stim_cur{I}};
    notes           = {notes{noSPEStype},notes{I}};
end


if isempty(s_start)
    s_start= 'n/a';
    s_end = 'n/a';
    duration = 'n/a';
    type = 'n/a';
    sub_type = 'n/a';
    ch_name_on ='n/a';
    ch_name_off = 'n/a';
    samp_start ='n/a';
    samp_end ='n/a';
    stim_type ='n/a';
    site_name='n/a';
    site_channum='n/a';
    stim_cur='n/a';
    notes='n/a';
end

annotation_tsv  = table(s_start', s_end', duration', type', sub_type', ch_name_on', ch_name_off', samp_start', samp_end', stim_type', site_name', site_channum',stim_cur', notes',  ...
    'VariableNames',{'onset', 'offset','duration','trial_type', 'sub_type','electrodes_involved_onset','electrodes_involved_offset','sample_start','sample_end','electrical_stimulation_type','electrical_stimulation_site','electrical_stimulation_site_num','electrical_stimulation_current','notes' });

if ~isempty(annotation_tsv)
    [p, f, x] = fileparts(cfg.outputfile);
    g = strsplit(f,'_ieeg');
    
    filename = fullfile(p, [g{1} '_events.tsv']);
    %filename = replace(filename,'_task-acute','')
    if isfile(filename)
        existing = read_tsv(filename);
    else
        existing = [];
    end % try
    if ~isempty(existing)
        ft_error('existing file is not empty');
    end
    write_tsv(filename, annotation_tsv);
end


%% extract all metadata needed for bid structure

% annots - annotations of the trc file
% ch     - channel labels of all channels in the trc file

function [status,msg,metadata]=extract_metadata_from_annotations(header,annots,ch,trigger,patName,cfg) % used on line 40
try
    status=0;
    metadata=[];
    
    %Codes for start and stop of good segments
    %     BEG_GS=222;
    %     END_GS=223;
    %
    %     ART_Start='xxx';
    %     ART_STOP='yyy';
    %
    %     ODD_Start='vvv';
    %     ODD_STOP='www';
    
    if (~isempty(trigger))
        trig_pos=trigger(1,:);
        trig_v=trigger(2,:);
    else
        trig_pos = [];
        trig_v = [];
    end
    
    %% Check the compulsory fields
    % Included; markers; situation name;Bad;(Bad field can appear more than once)
    % Resected;Edges;Format
    
    % Session
    ses_idx=cellfun(@(x) contains(x,{'Session'}),annots(:,2));
    
    if(sum(ses_idx)~=1)
        %       status=1;
        warning('Missing Session annotation (example "Session;1"), so session1 is set')
        metadata.ses_name='1';
    else
        str2parse=annots{ses_idx,2};
        %metadata.ses_name=strsplit(str2parse,'n');
        C=strsplit(str2parse,';');
        metadata.ses_name=C{2};
    end
    
    % Run
    run_idx=cellfun(@(x) contains(x,{'Run'}),annots(:,2));
    
    if(sum(run_idx)~=1)
        status=1;
        error('Missing run annotation (example "day1") or too many run annotations')
    end
    str2parse=annots{run_idx,2};
    C=strsplit(str2parse,';');
    D = strsplit(C{2},'y');
    if str2num(D{2}) <10
        metadata.run_name=['0' D{2}];
    else
        metadata.run_name=D{2};
    end
    
    % task
    task_idx=cellfun(@(x) contains(x,{'Task'}),annots(:,2));
    
    if(sum(task_idx)~=1)
        status=1;
        error('Missing task annotation (example "Task;SPES") or too many task annotations')
    end
    str2parse=annots{task_idx,2};
    C=strsplit(str2parse,';');
    metadata.task_name=C{2};
    
    % Stimcurr unknown (in SPES)
    stimcur_idx=cellfun(@(x) contains(x,{'Stimcurr'}),annots(:,2));
    
    if(sum(stimcur_idx)~=1)
        metadata.stimcurr = [];
    else
        str2parse=annots{stimcur_idx,2};
        %metadata.ses_name=strsplit(str2parse,'n');
        C=strsplit(str2parse,';');
        metadata.stimcurr=C{2};
    end
    
    % useful channels
    included_idx=cellfun(@(x) contains(x,{'Included'}),annots(:,2));
    metadata.ch2use_included= false(size(ch));
    if(sum(included_idx))
        metadata.ch2use_included=single_annotation(annots,'Included',ch);
        fprintf('File had Included-annotation, so no electrodes.tsv is used\n')
        metadata.incl_exist = 1;
    else % if "Included" is not annotated in the ECoG, there should be a previous ECoG with annoted "Included"
        metadata.incl_exist = 0;
        files_cECoG = dir(fullfile(cfg(2).proj_dirinput,patName,['ses-',metadata.ses_name],'ieeg',[patName, '_ses-',metadata.ses_name,'_electrodes.tsv']));
        if ~isempty(files_cECoG)
            files = files_cECoG;
            fprintf('%s/%s/ses-%s/ieeg/%s_ses-%s_electrodes.tsv is used.\n',cfg(2).proj_dirinput,patName,metadata.ses_name,patName,metadata.ses_name)
        end
        if ~isempty(files)
            elecName = fullfile(files(1).folder, '/',files(1).name);
            cc_elecs = readtable(elecName,'FileType','text','Delimiter','\t');
            elec_incl = false(size(ch));
            elec_soz = false(size(ch));
            elec_silicon = false(size(ch));
            elec_resected = false(size(ch));
            elec_edge = false(size(ch));
            % excluding electrodes other (so no grid, strip, depth)
            elec_name = cc_elecs.name;
            for i=1:size(ch,1)
                stimname = regexp(ch{i},'[a-z_A-Z]*','match');
                stimnum = regexp(ch{i},'[0-9]*','match');
                if ~isempty(stimname) && ~isempty(stimnum)
                    test1 = sprintf('%s0%d',stimname{:}, str2double(stimnum{:}));
                    test2 = sprintf('%s%d',stimname{:}, str2double(stimnum{:}));
                else
                    test1 = ch{i};
                    test2 = ch{i};
                end
                
                if sum(cellfun(@(x) strcmp(x,test1),elec_name)) == 1
                    idx = cellfun(@(x) strcmp(x,test1),elec_name);
                    elec_incl(i) = ~strcmp(cc_elecs.group{idx},'other');
                    elec_soz(i) = strcmp(cc_elecs.soz{idx},'yes');
                    elec_silicon(i) = strcmp(cc_elecs.silicon{idx},'yes');
                    elec_resected(i) = strcmp(cc_elecs.resected{idx},'yes');                    
                    elec_edge(i) = strcmp(cc_elecs.edge{idx},'yes');
                elseif sum(cellfun(@(x) strcmp(x,test2),elec_name)) == 1
                    idx = cellfun(@(x) strcmp(x,test2),elec_name);
                    elec_incl(i) = ~strcmp(cc_elecs.group{idx},'other');                    
                    elec_soz(i) = strcmp(cc_elecs.soz{idx},'yes');
                    elec_silicon(i) = strcmp(cc_elecs.silicon{idx},'yes');
                    elec_resected(i) = strcmp(cc_elecs.resected{idx},'yes');                    
                    elec_edge(i) = strcmp(cc_elecs.edge{idx},'yes');
                else
                    elec_incl(i) = false;
                    elec_soz(i) = false;
                    elec_silicon(i) = false;
                    elec_resected(i) = false;
                    elec_edge(i) = false;
                end
            end
            
            metadata.ch2use_included=logical(elec_incl);
            
        else
            error('There is no ECoG with annotated Included')
        end
    end
    
    % markers start and stop good segments
    %     begins=find(trig_v==BEG_GS);
    %     ends=find(trig_v==END_GS);
    %     if(isempty(begins) || isempty(ends) )
    %         status=1;
    %         error('Missing markers for good segments %i %i',BEG_GS,END_GS);
    %     end
    %     if(length(begins)~=length(ends))
    %         status=1;
    %         error('Missing start or stop Trigger');
    %     end
    
    %     for i=1:numel(begins)
    %         if(~issorted([trig_pos(begins(i)) trig_pos(ends(i))],'ascend'))
    %             status = 1;
    %             error('Trigger are not consecutive')
    %         end
    %     end
    %
    
    %% Look for bad channels
    metadata.ch2use_bad=single_annotation(annots,'Bad;',ch); % without the semicolon, the bad_HF channels are also included in Bad
    
    % cavity and silicon are not onmi present
    %     %% Look for cavity
    %     cavity_idx=cellfun(@(x) contains(x,{'Cavity'}),annots(:,2));
    %     metadata.ch2use_cavity= false(size(ch));
    %     if(sum(cavity_idx))
    %         metadata.ch2use_cavity=single_annotation(annots,'Cavity',ch);
    %     end
    
    %% Look for bad channels in high frequency band
    badhf_idx = cellfun(@(x) contains(x,{'Bad_HF'}),annots(:,2));
    metadata.ch2use_badhf= false(size(ch));
    if(sum(badhf_idx))
        metadata.ch2use_badhf=single_annotation(annots,'Bad_HF',ch);
        if any(contains(annots(:,2),'NB BadHF annotated in avg'))
            metadata.ch2use_badhf.note = 'NB BadHF annotated in avg';
        end
    end
    
    %% Look for silicon
    silicon_idx=cellfun(@(x) contains(x,{'Silicon'}),annots(:,2));
    metadata.ch2use_silicon= false(size(ch));
    if(sum(silicon_idx))
        metadata.ch2use_silicon=single_annotation(annots,'Silicon',ch);
    elseif metadata.incl_exist == 0
        metadata.ch2use_silicon=elec_silicon;
    end
    
    %% look for resected channels
    resected_idx = cellfun(@(x) contains(x,{'RA'}),annots(:,2));
    metadata.ch2use_resected= false(size(ch));
    if(sum(resected_idx))
        metadata.ch2use_resected=single_annotation(annots,'RA',ch);
    elseif metadata.incl_exist == 0
        metadata.ch2use_resected=elec_resected;
    end
    
    %% look for edge channels
    edge_idx = cellfun(@(x) contains(x,{'Edge'}),annots(:,2));
    metadata.ch2use_edge= false(size(ch));
    if(sum(edge_idx))
        metadata.ch2use_edge=single_annotation(annots,'Edge',ch);
    elseif metadata.incl_exist == 0
        metadata.ch2use_edge=elec_edge;
    end
    
    %% look for SOZ channels
    soz_idx = cellfun(@(x) contains(x,{'SOZ'}),annots(:,2));
    metadata.ch2use_soz= false(size(ch));
    if(sum(soz_idx))
        metadata.ch2use_soz=single_annotation(annots,'SOZ',ch);
    elseif metadata.incl_exist == 0
        metadata.ch2use_soz=elec_soz;
    end
    
    %% look for screw channels - only in seeg
    screw_idx = cellfun(@(x) contains(x,{'Screw'}),annots(:,2));
    metadata.ch2use_screw= false(size(ch));
    if(sum(screw_idx))
        metadata.ch2use_screw=single_annotation(annots,'Screw',ch);
    elseif metadata.incl_exist == 0
        if strcmp(fieldnames(cc_elecs),'screw')
            elec_screw= strcmp(cc_elecs.screw,'yes');
            metadata.ch2use_screw=logical(elec_screw);
        end            
    end

    %% look for white matter channels - only in seeg
    wm_idx = cellfun(@(x) contains(x,{'WM'}),annots(:,2));
    metadata.ch2use_wm= false(size(ch));
    if(sum(wm_idx))
        metadata.ch2use_wm=single_annotation(annots,'WM',ch);
    elseif metadata.incl_exist == 0
        if strcmp(fieldnames(cc_elecs),'wm')
            elec_wm= strcmp(cc_elecs.wm,'yes');
            metadata.ch2use_wm=logical(elec_wm);
        end
    end
    
    %% look for gray matter channels - only in seeg
    gm_idx = cellfun(@(x) contains(x,{'GM'}),annots(:,2));
    metadata.ch2use_gm= false(size(ch));
    if(sum(gm_idx))
        metadata.ch2use_gm=single_annotation(annots,'GM',ch);
    elseif metadata.incl_exist == 0
        if strcmp(fieldnames(cc_elecs),'gm')
            elec_gm= strcmp(cc_elecs.gm,'yes');
            metadata.ch2use_gm=logical(elec_gm);
        end
    end
    
    %% look for CSF channels - only in seeg
    csf_idx = cellfun(@(x) contains(x,{'CSF'}),annots(:,2));
    metadata.ch2use_csf= false(size(ch));
    if(sum(csf_idx))
        metadata.ch2use_csf=single_annotation(annots,'CSF',ch);
    elseif metadata.incl_exist == 0
        if strcmp(fieldnames(cc_elecs),'csf')
            elec_csf= strcmp(cc_elecs.csf,'yes');
            metadata.ch2use_csf=logical(elec_csf);
        end
    end
    
    %% look for amygdala channels - only in seeg
    amyg_idx = cellfun(@(x) contains(x,{'Amyg'}),annots(:,2));
    metadata.ch2use_amyg= false(size(ch));
    if(sum(amyg_idx))
        metadata.ch2use_amyg=single_annotation(annots,'Amyg',ch);
    elseif metadata.incl_exist == 0
        if strcmp(fieldnames(cc_elecs),'amyg')
            elec_amyg= strcmp(cc_elecs.amyg,'yes');
            metadata.ch2use_amyg=logical(elec_amyg);
        end
    end
    
    %% look for hippocampal channels - only in seeg
    hipp_idx = cellfun(@(x) contains(x,{'Hipp'}),annots(:,2));
    metadata.ch2use_hipp= false(size(ch));
    if(sum(hipp_idx))
        metadata.ch2use_hipp=single_annotation(annots,'Hipp',ch);
    elseif metadata.incl_exist == 0
        if strcmp(fieldnames(cc_elecs),'hipp')
            elec_hipp= strcmp(cc_elecs.hipp,'yes');
            metadata.ch2use_hipp=logical(elec_hipp);
        end
    end
    
    %% look for lesion channels - only in seeg
    lesion_idx = cellfun(@(x) contains(x,{'Lesion'}),annots(:,2));
    metadata.ch2use_lesion= false(size(ch));
    if(sum(lesion_idx))
        metadata.ch2use_lesion=single_annotation(annots,'Lesion',ch);
    elseif metadata.incl_exist == 0
        if strcmp(fieldnames(cc_elecs),'lesion')
            
            elec_lesion= strcmp(cc_elecs.lesion,'yes');
            metadata.ch2use_lesion=logical(elec_lesion);
        end
    end
     
    %% look for gliosis channels - only in seeg
    gliosis_idx = cellfun(@(x) contains(x,{'Glio'}),annots(:,2));
    metadata.ch2use_gliosis= false(size(ch));
    if(sum(gliosis_idx))
        metadata.ch2use_lesion=single_annotation(annots,'Glio',ch);
    elseif metadata.incl_exist == 0
        if strcmp(fieldnames(cc_elecs),'gliosis')
            
            elec_gliosis= strcmp(cc_elecs.gliosis,'yes');
            metadata.ch2use_gliosis=logical(elec_gliosis);
        end
    end
    
    %% Look for artefacts cECoG
    metadata.artefacts=look_for_annotation_start_stop(annots,'Art_on','Art_off',ch,header);
    
    %% Look for sleep data
    metadata.sleep=look_for_annotation_start_stop(annots,'Sl_on','Sl_off',ch,header);
    
    %% Look for data between rest and sleep
    metadata.slaw_trans = look_for_annotation_start_stop(annots,'Slawtrans_on','Slawtrans_off',ch,header);
    
    %% Look for data eyes opened
    metadata.eyes_open = look_for_annotation_start_stop(annots,'Eyes_open','Eyes_close',ch,header);
    
    %% Look for seizures
    metadata.seizure=look_for_annotation_start_stop(annots,'Sz_on','Sz_off',ch,header);
    
    %% Look for period of stimulation
    metadata.stimulation=look_for_annotation_start_stop(annots,'Stim_on','Stim_off',ch,header);
    
    %% Look for period of SPES
    metadata.spes=look_for_annotation_start_stop(annots,'SPES_on','SPES_off',ch,header);
    
    %% Look for period of Electrical Stimulation Mapping(ESM)
    metadata.esm=look_for_annotation_start_stop(annots,'ESM_on','ESM_off',ch,header);
    
    %% Look for period of motor task
    metadata.motortask=look_for_annotation_start_stop(annots,'Motor_on','Motor_off',ch,header);
    
    %% Look for period of sens task
    metadata.senstask=look_for_annotation_start_stop(annots,'Sens_on','Sens_off',ch,header);
    
    %% Look for period of language task
    metadata.langtask=look_for_annotation_start_stop(annots,'Language_on','Language_off',ch,header);
    
     %% Look for SWS selection - slow wave sleep
    metadata.SWSselection=look_for_annotation_start_stop(annots,'SWS10_on','SWS10_off',ch,header);
   
     %% Look for REM selection - rapid eye movement
    metadata.REMselection=look_for_annotation_start_stop(annots,'REM10_on','REM10_off',ch,header);
    
     %% Look for IIAW selection - inter ictal awake
    metadata.IIAWselection=look_for_annotation_start_stop(annots,'IIAW10_on','IIAW10_off',ch,header);
    
    %% Look for EI selection
    metadata.EIselection=look_for_annotation_start_stop(annots,'EI_on','EI_off',ch,header);

    %% Look for artefacts
    
    %metadata.artefacts_aECoG=look_for_annotation_start_stop(annots,'xxx','yyy',ch);
    
    %% look for odd behaviour in the recordings additional notes
    
    metadata.add_notes=look_for_annotation_start_stop(annots,'vvv','www',ch,header);
    
    %% look for burst suppression
    
    metadata.bsuppression=look_for_burst_suppression(annots);
    
    %% look for Format
    %TODO double check for the syntax
    
    format_idx=cellfun(@(x) contains(x,{'Format'}),annots(:,2));
    if(sum(format_idx)<1)
        %         status=1;
        %         error('Missing Format annotation (example "Format;Gr[5x4];")')
        % load format from another ECoG --> json file
        file = dir(fullfile(cfg(2).proj_dirinput,patName,['ses-',metadata.ses_name],'ieeg','*ieeg.json'));
        if ~isempty(file)
            ieeg_json = jsondecode(fileread([file(1).folder '/' file(1).name]) );
            metadata.format_info = ieeg_json.iEEGElectrodeGroups;
        else
            error('Missing Format annotation, and no other json from other ECoG with format annotation found')
        end
    else
        loc = find(format_idx==1);
        annots_format = cell(1,size(loc,1));
        for i=1:size(loc,1)
            annots_format_all = strsplit(annots{loc(i),2},'Format;');
            annots_format{i} = [annots_format_all{2} ';'];
        end
        
        % putting format in order ECoG,strip,depth
        annotsformatsplit = strsplit([annots_format{:}],';');
        annotsformatsplit = annotsformatsplit(~cellfun(@isempty,annotsformatsplit));
        ecogloc = find(cellfun('length',regexp(lower(annotsformatsplit),'ecog')) == 1);
        striploc = find(cellfun('length',regexp(lower(annotsformatsplit),'strip')) == 1);
        depthloc = find(cellfun('length',regexp(lower(annotsformatsplit),'depth')) == 1);
        seegloc = find(cellfun('length',regexp(lower(annotsformatsplit),'seeg')) == 1);
        
        locs_all = sort([ecogloc, striploc, depthloc,seegloc,size(annotsformatsplit,2)+1]);
        
        % ECoG
        if ~isempty(ecogloc)
            ecogformat = cell(1,size(ecogloc,2));
            for i=1:size(ecogloc,2)
                ecogformat{i} = annotsformatsplit(ecogloc(i)+1: locs_all(find(locs_all==ecogloc(i))+1)-1);
            end
            ecogformatall =  strcat([ecogformat{:}],';');
            ecogformatfin = ['ECoG;' ecogformatall{:}];
        else
            ecogformatfin = [];
        end
        
        % strip
        if ~isempty(striploc)
            stripformat = cell(1,size(striploc,2));
            for i=1:size(striploc,2)
                stripformat{i} = annotsformatsplit(striploc(i)+1: locs_all(find(locs_all==striploc(i))+1)-1);
            end
            stripformatall = strcat([stripformat{:}],';');
            stripformatfin = ['strip;' stripformatall{:}];
        else
            stripformatfin = [];
        end
        
        % depth
        if ~isempty(depthloc)
            depthformat = cell(1,size(depthloc,2));
            for i=1:size(depthloc,2)
                depthformat{i} = annotsformatsplit(depthloc(i)+1: locs_all(find(locs_all==depthloc(i))+1)-1);
            end
            depthformatall = strcat([depthformat{:}],';');
            depthformatfin = ['depth;' depthformatall{:}];
        else
            depthformatfin = [];
        end
        
        % seeg
        if ~isempty(seegloc)
            seegformat = cell(1,size(seegloc,2));
            for i=1:size(seegloc,2)
                seegformat{i} = annotsformatsplit(seegloc(i)+1: locs_all(find(locs_all==seegloc(i))+1)-1);
            end
            seegformatall = strcat([seegformat{:}],';');
            seegformatfin = ['seeg;' seegformatall{:}];
        else
            seegformatfin = [];
        end
        
        metadata.format_info=[ecogformatfin, stripformatfin, depthformatfin, seegformatfin];
    end
    
    %SEEG/ECoG?
    if contains(lower(metadata.format_info),'seeg')
        metadata.elec_info = 'SEEG';
    elseif contains(lower(metadata.format_info),'ecog')
        metadata.elec_info = 'ECoG';
    end
    
    %% add triggers
    if ~isempty(trigger)
        metadata.trigger.pos  = trigger(1,:)  ;
        metadata.trigger.val  = trigger(end,:);
    else
        metadata.trigger.pos  = [];
        metadata.trigger.val  = [];
    end
    
    %% add channel labels
    
    metadata.ch_label = ch;
    
    status = 0 ;
    msg    = '';
    %
catch ME
    status = 1;
    msg = sprintf('%s err:%s --func:%s',deblank(patName'),ME.message,ME.stack(1).name);
    
end

function [artefacts]=look_for_annotation_start_stop(annots,str_start,str_stop,ch,header)

fs = header.Rate_Min;

start_art=find(contains(annots(:,2),str_start));
end_art=find(contains(annots(:,2),str_stop));

% in a seizure, it is possible that the involved electrodes cannot be
% described in one Sz_on or one Sz_off, this part deals with this
if strcmp(str_start,'Sz_on')
    
    sz_cont = find(contains(annots(:,2),'Sz_cont'), 1);
    if ~isempty(sz_cont)
        if annots{sz_cont,1} < 10*fs % Sz_cont is annotated at beginning of file, sz_cont should be added in start_art
                start_art = sort([start_art; sz_cont]);
        else
            end_art = sort([end_art; sz_cont]);
        end
    end
        
    
    if any(diff(end_art) == 1) % if two Sz_offs are after each other
        diffSzoffs = (diff(end_art));
        sampSzOffs = diff([annots{end_art,1}]);
        
        for i=1:size(sampSzOffs,2)
            if diffSzoffs(i) == 1 && sampSzOffs(i) < 1*fs % SzOffs belong to the same seizure
                if strcmp(annots{end_art(i),2}(end),';')
                    annot1 = annots{end_art(i),2};
                else
                    annot1 = [annots{end_art(i),2},';'];
                end
                annotsplit2 = strsplit(annots{end_art(i+1),2},'Sz_off;');
                annots{end_art(i),2} = '';
                annots{end_art(i+1),2} = [annot1, annotsplit2{2}];
            elseif diffSzoffs(i) == 1 && sampSzOffs(i) > 1*fs % two Sz_offs do not belong to the same seizure
                error('starts and ends did not match')
            end
        end
        
        end_art=find(contains(annots(:,2),str_stop));
    end
    
    if any(diff(start_art) == 1) % if two Sz_ons are after each other
        diffSzons = (diff(start_art));
        sampSzOns = diff([annots{start_art,1}]);
        
        for i=1:size(sampSzOns,2)
            if diffSzons(i) == 1 && sampSzOns(i) < 1*fs % SzOffs belong to the same seizure
                if strcmp(annots{start_art(i),2}(end),';')
                    annot1 = annots{start_art(i),2};
                else
                    annot1 = [annots{start_art(i),2},';'];
                end
                annotsplit2 = strsplit(annots{start_art(i+1),2},'Sz_on;');
                annots{start_art(i),2} = '';
                annots{start_art(i+1),2} = [annot1, annotsplit2{2}];
            elseif diffSzons(i) == 1 && sampSzOns(i) > 1*fs % two Sz_offs do not belong to the same seizure
                error('starts and ends did not match')
            end
        end
        start_art=find(contains(annots(:,2),str_start));
    end
end

if (length(start_art)~=length(end_art)) && ~strcmp(str_start,'Eyes_open')
    error('starts and ends did not match')
end

if strcmp(str_start,'Art_on')
    %this code works perfectly when an event is closed in which the same
    %electrodes are involved in 'on' and 'off', like artefact. This is not the situation
    %in case of a seizure. Different electrodes can be involved in the end
    %compared to the beginning. In case of a motortask/sleep, a subtype is
    %mentioned in '.._on' but not in '.._off'
    artefacts=cell(size(start_art));
    
    for i=1:numel(start_art)
        art=struct;
        matched_end=find(contains(annots(:,2),replace(annots{start_art(i),2},str_start,str_stop)));
        if(isempty(matched_end))
            error('start and stop %s does not match',annots{start_art(i),2});
        end
        if(length(matched_end)>1)
            matched_end=matched_end((matched_end-start_art(i))>0);
            [val,idx_closest]=min(matched_end);
            matched_end=matched_end(idx_closest);%take the closest in time
        end
        ch_art_idx=parse_annotation(annots{start_art(i),2},ch);
        
        if isempty({ch{logical(ch_art_idx)}})
            art.ch_names_on={'all'};
            art.ch_names_off={'all'};
        else
            art.ch_names_on={ch{logical(ch_art_idx)}};
            art.ch_names_off={ch{logical(ch_art_idx)}};
        end
        art.pos=[(annots{start_art(i),1}) annots{matched_end,1}];
        art.type = 'n/a';
        if any(contains(annots(:,2),'NB artefacts annotated in avg'))
            art.notes = 'Please note, artefacts annotated in avg!';
        else
            art.notes = 'n/a';
        end
        artefacts{i}=art;
    end
elseif  strcmp(str_start,'Sz_on') % in case of a seizure
    artefacts=cell(size(start_art));
    for i=1:numel(start_art)
        art=struct;
        
        % seizure onset
        annotsplit = strsplit(annots{start_art(i),2},';');
        annotsplit = annotsplit(~cellfun(@isempty,annotsplit));
        if size(annotsplit,2) >2 % both subtype of seizure and channelnames are mentioned 
            if strcmp(annotsplit{2},'clin') || strcmp(annotsplit{2},'cluster') || strcmp(annotsplit{2},'subclin') || strcmp(annotsplit{2},'aura')
                type = annotsplit{2};
            else
                type = 'unknown';
            end
            ch_art_idx = parse_annotation(annots{start_art(i),2},ch);
            ch_names_on = {ch{logical(ch_art_idx)}};
            
        elseif size(annotsplit,2) == 2 % either subtype of seizure or channelnames are mentiond
            if strcmp(annotsplit{2},'clin') || strcmp(annotsplit{2},'cluster') || strcmp(annotsplit{2},'subclin') || strcmp(annotsplit{2},'aura')
                type = annotsplit{2};
                
                if strcmp(annotsplit{1},'Sz_cont') % Sz_cont is annotated, so no channelnames
                    ch_names_on = {'continuation of seizure'};
                else
                    ch_names_on = {'diffuse'};
                end
            else
                type = 'unknown';
                ch_art_idx = parse_annotation(annots{start_art(i),2},ch);
                ch_names_on = {ch{logical(ch_art_idx)}};
            end
        elseif  size(annotsplit,2) ==1 && strcmp(annotsplit{1},'Sz_cont') % for example with Sz_cont;, it is possible that only Sz_cont is mentioned
            type = 'unknown';
            ch_names_on = {'continuation of seizure'};
        end
        
        % seizure offset
        annotsplit = strsplit(annots{end_art(i),2},';');
        annotsplit = annotsplit(~cellfun(@isempty,annotsplit));
        if size(annotsplit,2) > 1 % for example "Sz_off;C[13,14]
            ch_art_idx = parse_annotation(annots{end_art(i),2},ch);
            ch_names_off = {ch{logical(ch_art_idx)}};
        else
            if strcmp(annotsplit{1},'Sz_cont')
                ch_names_off = {'Seizure continues in next file'};
            else
                ch_names_off = {'diffuse'};
            end
        end
        
        art.type = type;
        art.ch_names_on = ch_names_on;
        art.ch_names_off = ch_names_off;
        art.pos = [(annots{start_art(i),1}) (annots{end_art(i),1})]; % each seizure is consecutive, so as long as a seizure's end has the correct annotation, this will work
        artefacts{i} = art;
    end
elseif  strcmp(str_start,'Eyes_open') % in case of eyes open
    all_art = sort([start_art;end_art]);
    artefacts=cell(size(all_art));
    
    for i=1:numel(all_art)
        art=struct;
       
        annotsplit = strsplit(annots{all_art(i),2},';');
        annotsplit = annotsplit(~cellfun(@isempty,annotsplit));
        
        if strcmp(annotsplit{1},'Eyes_open')
            type = 'opened';
        else % if annotsplit is Eyes_close
            type = 'closed';
        end        
        art.type = type;
        
        art.ch_names_on={'all'};
        art.ch_names_off={'all'};
        
        if i<numel(all_art)
            art.pos=[(annots{all_art(i),1}) annots{all_art(i+1),1}];
        elseif i==numel(all_art) % if last annotation is selected, this event takes till the end of the file
            art.pos=[(annots{all_art(i),1}) header.Num_Samples];
        end
        artefacts{i}=art;
    end
    
else 
    artefacts=cell(size(start_art));
    
    for i=1:numel(start_art)
        art=struct;
        
        %matched_end=find(contains(annots(:,2),replace(annots{start_art(i),2},str_start,str_stop)));
        %         if(isempty(matched_end))
        %             error('start and stop %s does not match',annots{start_art(i),2});
        %         end
        %         if(length(matched_end)>1)
        %             matched_end=matched_end((matched_end-start_art(i))>0);
        %             [val,idx_closest]=min(matched_end);
        %             matched_end=matched_end(idx_closest);%take the closest in time
        %         end
        
        annotsplit = strsplit(annots{start_art(i),2},';');
        annotsplit = annotsplit(~cellfun(@isempty,annotsplit));
        
        if size(annotsplit,2) == 1
            if strcmp(str_start,'Eyes_open')
                type = 'n/a';
            else % if str_start is motor/sens/language
                type = 'unknown';
            end
        elseif size(annotsplit,2) == 2
            type = annotsplit{2};
        end
        
        ch_art_idx=parse_annotation(annots{start_art(i),2},ch);
        
        art.type = type;
        if isempty({ch{logical(ch_art_idx)}})
            art.ch_names_on={'all'};
            art.ch_names_off={'all'};
        else
            art.ch_names_on={ch{logical(ch_art_idx)}};
            art.ch_names_off={ch{logical(ch_art_idx)}};
        end
        art.pos=[(annots{start_art(i),1}) annots{end_art(i),1}];
        artefacts{i}=art;
    end
end



function [annots_new, eventsannots ] = add_event2annotation(event,evname,eventsannots, annots_new, header)

fs = header.Rate_Min;

if size(eventsannots.type,2) == 0
    cc = 1;
else
    cc = size(eventsannots.type,2)+1;
end

if(~isempty(event))
    for i=1:numel(event)
        
        eventsannots.type{cc}    = evname                           ;
        
        %         if strcmp(evname,'artefact')
        eventsannots.sub_type{cc} = event{i}.type;
        eventsannots.s_start{cc} = round(event{i}.pos(1)/fs,1); % time in seconds (1 decimal)
        eventsannots.samp_start{cc} = num2str(event{i}.pos(1))          ;
        eventsannots.s_end{cc}   = round(event{i}.pos(end)/fs,1); % time in seconds (1 decimal)
        eventsannots.samp_end{cc} = num2str(event{i}.pos(end))          ;
        eventsannots.duration{cc} = round(eventsannots.s_end{cc} - eventsannots.s_start{cc},3);
        eventsannots.stim_type{cc} = 'n/a';
        eventsannots.site_name{cc} = 'n/a';
        eventsannots.site_channum{cc} = 'n/a';
        eventsannots.stim_cur{cc} = 'n/a';
       
        if isfield(event{i},'notes') && ~isempty(event{i}.notes)
            eventsannots.notes{cc} = event{i}.notes;
        else
            eventsannots.notes{cc} = 'n/a';
        end
        
        if size(event{i}.ch_names_on,2) == 1
            name = event{i}.ch_names_on{1}              ;
        else
            for j=1
                name =[event{i}.ch_names_on{1}];
            end
            
            for j=2:size(event{i}.ch_names_on,2)
                name = [name ,',', event{i}.ch_names_on{j}];
            end
        end
        eventsannots.ch_name_on{cc} = name;
        
        if size(event{i}.ch_names_off,2) == 1
            name = event{i}.ch_names_off{1}              ;
        else
            for j=1
                name =[event{i}.ch_names_off{1}];
            end
            
            for j=2:size(event{i}.ch_names_off,2)
                name = [name ,',', event{i}.ch_names_off{j}];
            end
        end
        
        eventsannots.ch_name_off{cc} = name;
        
        annots_new([annots_new{:,1}]==event{i}.pos(1),:)=[];
        annots_new([annots_new{:,1}]==event{i}.pos(end),:)=[];
        
        cc          = cc + 1                               ;
        
        %         end
    end
end


function samplocs = findtrigger(data, stimnum,fs, sampstart,sampend)

if median(data(stimnum(1),sampstart:sampend)) > median(data(stimnum(2),sampstart:sampend))
    anode = stimnum(2);% with positive saturation right after stimulation
    cathode = stimnum(1);% with negative saturation right after stimulation
else
    anode = stimnum(1);% with positive saturation right after stimulation
    cathode = stimnum(2);% with negative saturation right after stimulation
end

dt_data(1,:) = diff(data(anode,sampstart:sampend));
dt_data(2,:) = diff(data(cathode,sampstart:sampend));
SD = std(dt_data,[],2);

[~,locsdtpos] = findpeaks(dt_data(1,:),'MinPeakDistance',fs,'MinPeakHeight',10*SD(1));
[~,locsdtneg] = findpeaks(-1*dt_data(1,:),'MinPeakDistance',fs,'MinPeakHeight',10*SD(1));
 
[~,locsdt2pos] = findpeaks(dt_data(2,:),'MinPeakDistance',fs,'MinPeakHeight',10*SD(2));
[~,locsdt2neg] = findpeaks(-1*dt_data(2,:),'MinPeakDistance',fs,'MinPeakHeight',10*SD(2));

locsall = sort([locsdtpos locsdtneg, locsdt2pos, locsdt2neg]);
distances = pdist2(locsall',locsall');

within1 = distances < round(0.5*fs); % group all locs that are within 500ms
numelements = numel(locsall);
% use only lower triangle part
t = logical(triu(ones(numelements,numelements),1)); 
within1(t) = 0;
[labeledGroups,numGroups ] = bwlabel(within1,4);

for k = 1 : numGroups
  [rows, columns] = find(labeledGroups == k);
  indexes = unique([rows, columns]);
  if k>1
      if locsall(indexes(1))-samplocs(k-1) < round(6*fs) % WARNING: assumes that stimuli are applied maximally every 5s
          samplocs(k) = min(locsall(indexes));
      else
          break
      end
  else
      samplocs(k) = min(locsall(indexes));
  end
end




function bsuppression=look_for_burst_suppression(annots)

BS_Start='200';
BS_Stop='201';

start_bs=find(startsWith(annots(:,2),BS_Start));
end_bs=find(startsWith(annots(:,2),BS_Stop));

if(length(start_bs)~=length(end_bs))
    error('burst suppression: starts and ends did no match')
end

bsuppression=cell(size(start_bs));

for i=1:numel(start_bs)
    bs=struct;
    matched_end=find(contains(annots(:,2),BS_Stop));
    if(isempty(matched_end))
        error('start and stop %s does not match',annots{start_bs(i),2});
    end
    if(length(matched_end)>1)
        matched_end=matched_end((matched_end-start_bs(i))>0);
        [val,idx_closest]=min(matched_end);
        matched_end=matched_end(idx_closest);%take the closest in time
    end
    
    bs.pos=[(annots{start_bs(i),1}) annots{matched_end,1}];
    bsuppression{i}=bs;
end

function [ch_parsed]=single_annotation(annots,keyWord,ch)


ch_idx=cellfun(@(x) contains(x,{keyWord}),annots(:,2));

if(sum(ch_idx)<1)
    error('Missing annotation (example "%s;Gr01;Gr[3:5]")',keyWord)
end
ch_parsed=zeros(size(ch));
if(sum(ch_idx))
    str2parse={annots{ch_idx,2}};
    for i=1:numel(str2parse)
        C=strsplit(str2parse{i},';');
        C=C(~cellfun(@isempty,C));
        if(numel(C)>1)%TODO better check
            ch_parsed= ch_parsed | parse_annotation(str2parse{i},ch);
        end
    end
end

function mydirMaker(dirname)
if exist(dirname, 'dir')
    warning('%s exist already',dirname)
else
    mkdir(dirname)
end

function [ch_status,ch_status_desc]=status_and_description(metadata)

ch_label = metadata.ch_label                         ;
ch2use_included = metadata.ch2use_included; 
ch2use_silicon = metadata.ch2use_silicon;
% in some patients in some files an extra headbox is included in a later
% stadium because an EMG-recording was necessary (for example).
% if size(metadata.ch2use_included,1) < size(metadata.ch_label,1)
%     ch2use_included = zeros(size(ch_label));
%     ch2use_included(1:size(metadata.ch2use_included,1)) =  metadata.ch2use_included;
%     ch2use_included = logical(ch2use_included);
%     ch2use_silicon = zeros(size(ch_label));
%     ch2use_silicon(1:size(metadata.ch2use_silicon,1)) = metadata.ch2use_silicon;
%     ch2use_silicon = logical(ch2use_silicon);
%     
% elseif size(metadata.ch2use_included,1) == size(metadata.ch_label,1)
%     ch2use_included = metadata.ch2use_included;
%     ch2use_silicon = metadata.ch2use_silicon;
% elseif size(metadata.ch2use_included,1) > size(metadata.ch_label,1)
%     ch2use_included =  metadata.ch2use_included(1:size(metadata.ch_label,1));
%     ch2use_included = logical(ch2use_included);
%     ch2use_silicon = metadata.ch2use_silicon(1:size(metadata.ch_label,1));
%     ch2use_silicon = logical(ch2use_silicon);
%     
% end

ch_status                                                       = cell(size(ch2use_included))      ;
ch_status_desc                                                  = cell(size(ch2use_included))      ;

idx_ecg                                                         = ~cellfun(@isempty,regexpi(ch_label,'ECG'));
idx_ecg                                                         = idx_ecg                                  ;
idx_mkr                                                         = ~cellfun(@isempty,regexpi(ch_label,'MKR'));
idx_mkr                                                         = idx_mkr                                  ;
% channels which are open but not recording
ch_open                                                         = ~(ch2use_included | ...
    metadata.ch2use_bad      | ...
    ch2use_silicon  | ...
    idx_ecg                  | ...
    idx_mkr                    ...
    )                                       ;
%     metadata.ch2use_cavity   | ...

[ch_status{:}]                                                  = deal('good')                              ;

if(any(metadata.ch2use_bad              ...
        )) % removed metadata.ch2use_cavity  | ...
    
    [ch_status{(metadata.ch2use_bad     ...
        )}] = deal('bad'); % removed metadata.ch2use_cavity  | ...
end

% bad in high frequency band
if(any(metadata.ch2use_badhf)) % removed metadata.ch2use_cavity  | ...
    
    [ch_status{(metadata.ch2use_badhf     ...
        )}] = deal('bad_hf'); % removed metadata.ch2use_cavity  | ...
end

if (any(ch_open))
    [ch_status{ch_open}] = deal('bad');
end

%% status description
if(any(ch2use_included))
    [ch_status_desc{ch2use_included}] = deal('included');
end

if(any(metadata.ch2use_bad))
    [ch_status_desc{metadata.ch2use_bad}] = deal('noisy (visual assessment)');
end

if(any(metadata.ch2use_badhf))
    if strcmp(metadata.ch2use_badhf.note,'NB BadHF annotated in avg')
        [ch_status_desc{metadata.ch2use_badhf}] = deal('noisy in high frequency bands >80Hz&<500 Hz(visual assessment). Please note, BadHF was annotated in avg!');
    else
        [ch_status_desc{metadata.ch2use_badhf}] = deal('noisy in high frequency bands >80Hz&<500 Hz(visual assessment)');
    end
end

% if(any(metadata.ch2use_cavity))
%     [ch_status_desc{metadata.ch2use_cavity}] = deal('cavity');
% end
%
% silicon information not in channels.tsv but in electrodes.tsv
% if(any(metadata.ch2use_silicon))
%     [ch_status_desc{metadata.ch2use_silicon}] = deal('silicon');
% end

if(any(ch_open))
    [ch_status_desc{ch_open}] = deal('not recording');
end

if(sum(idx_ecg))
    [ch_status_desc{idx_ecg}] = deal('not included');
end
if(sum(idx_mkr))
    [ch_status_desc{idx_mkr}] = deal('not included');
end

% extract group information
% assumption the included are only grid and strip
function ch_group = extract_group_info(metadata)

ch_label                                    = metadata.ch_label                    ;
ch2use_included = metadata.ch2use_included;

if strcmpi(metadata.elec_info,'SEEG')
    idx_depths = ch2use_included;
    idx_strips = false(size(ch_label));
    idx_grid = false(size(ch_label));
elseif strcmpi(metadata.elec_info,'ECoG')
    C = strsplit(metadata.format_info,{';','['});
    
    id_ngrid = regexpi(C,'1');
    ngridnum = find(cellfun(@isempty,id_ngrid)==0);
    
    id_strip = regexpi(C,'strip');
    stripnum = find(cellfun(@isempty,id_strip)==0);
    id_depth = regexpi(C,'depth');
    depthnum = find(cellfun(@isempty,id_depth)==0);
    
    idx_depths = zeros(size(ch_label));
    idx_strips = zeros(size(ch_label));
    for i=1:size(ngridnum,2)
        if ~isempty(stripnum) && isempty(depthnum)
            if any(ngridnum(i) > stripnum)
                idx_strip = regexpi(ch_label,C{ngridnum(i)-1});
                idx_strip = cellfun(@isempty,idx_strip);
                idx_strip = ~idx_strip;
                idx_strips = idx_strips + idx_strip;
            end
        elseif isempty(stripnum) && ~isempty(depthnum)
            if any(ngridnum(i) > depthnum)
                idx_depth = regexpi(ch_label,C{ngridnum(i)-1});
                idx_depth = cellfun(@isempty,idx_depth);
                idx_depth = ~idx_depth;
                idx_depths = idx_depths + idx_depth;
            end
        elseif ~isempty(stripnum) && ~isempty(depthnum)
            if depthnum < stripnum
                if any(ngridnum(i) > depthnum) && any(ngridnum(i) < stripnum)
                    idx_depth = regexpi(ch_label,C{ngridnum(i)-1});
                    idx_depth = cellfun(@isempty,idx_depth);
                    idx_depth = ~idx_depth;
                    idx_depths = idx_depths + idx_depth;
                elseif ngridnum(i) > stripnum
                    idx_strip = regexpi(ch_label,C{ngridnum(i)-1});
                    idx_strip = cellfun(@isempty,idx_strip);
                    idx_strip = ~idx_strip;
                    idx_strips = idx_strips + idx_strip;
                end
            elseif depthnum >stripnum
                if ngridnum(i) > stripnum && ngridnum(i) < depthnum
                    idx_strip = regexpi(ch_label,C{ngridnum(i)-1});
                    idx_strip = cellfun(@isempty,idx_strip);
                    idx_strip = ~idx_strip;
                    idx_strips = idx_strips + idx_strip;
                    
                elseif ngridnum(i) > depthnum
                    idx_depth = regexpi(ch_label,C{ngridnum(i)-1});
                    idx_depth = cellfun(@isempty,idx_depth);
                    idx_depth = ~idx_depth;
                    idx_depths = idx_depths + idx_depth;
                end
            end
        end
    end
    idx_grid = ~idx_depths & ~idx_strips & ch2use_included;
end

idx_depths = logical(idx_depths);
idx_strips = logical(idx_strips);

ch_group                                    = cell(size(ch2use_included)) ;
if(any(idx_grid))
    [ch_group{idx_grid}]                    = deal('grid')                         ;
end
if(any(idx_strips))
    [ch_group{idx_strips}]                   = deal('strip')                        ;
end
if(any(idx_depths))
    [ch_group{idx_depths}]                   = deal('depth')                        ;
end
if(any(~ch2use_included))
    [ch_group{ ~ch2use_included }] = deal('other')                        ;
end

function write_scans_tsv(cfg,metadata,annotation_tsv)
[p,f] = fileparts(cfg(1).outputfile);
g = strsplit(f,'_');
filename = fullfile(p,[g{1},'_',g{2},'_scans.tsv']);

files = dir(p);
if contains([files(:).name],'scans')
    
    % read existing scans-file
    scans_tsv = read_tsv(filename);
    
    if any(contains(scans_tsv.name,f))
        scansnum = find(contains(scans_tsv.name,f) ==1);
    else
        scansnum = size(scans_tsv,1)+1;
    end
    
    name                    = scans_tsv.name;
    artefact                = scans_tsv.artefact;
    sleep_total             = scans_tsv.sleep_total;
    sleep_rem               = scans_tsv.sleep_rem;
    sleep_nrem              = scans_tsv.sleep_nrem;
    seizure                 = scans_tsv.seizure_total;
    seizuresubclin          = scans_tsv.seizure_subclinical;
    seizureclin             = scans_tsv.seizure_clinical;
    motor                   = scans_tsv.motor;
    spes                    = scans_tsv.spes;
    esm                     = scans_tsv.esm;
    language                = scans_tsv.language;
    sleepwaketransition     = scans_tsv.sleepwaketransition;
    format                  = scans_tsv.format;
 
    if any(contains(scans_tsv.Properties.VariableNames,'sens')) % this was added later, so not allf iles have this
        sens                    = scans_tsv.sens;
    else
        sens = zeros(size(name,1),1);
    end
    
    if any(contains(scans_tsv.Properties.VariableNames,'slowesm')) % this was added later, so not allf iles have this
        slowesm = scans_tsv.slowesm;
    else
        slowesm = zeros(size(name,1),1);
    end
    
    if any(contains(scans_tsv.Properties.VariableNames,'sws_sel')) % this was added later, so not allf iles have this
        sws_sel = scans_tsv.sws_sel;
    else
        sws_sel = zeros(size(name,1),1);
    end
    
    if any(contains(scans_tsv.Properties.VariableNames,'rem_sel')) % this was added later, so not allf iles have this
        rem_sel = scans_tsv.rem_sel;
    else
        rem_sel = zeros(size(name,1),1);
    end
    
    if any(contains(scans_tsv.Properties.VariableNames,'iiaw_sel')) % this was added later, so not allf iles have this
        iiaw_sel = scans_tsv.iiaw_sel;
    else
        iiaw_sel = zeros(size(name,1),1);
    end
    
    if any(contains(scans_tsv.Properties.VariableNames,'EI_sel')) % this was added later, so not allf iles have this
        EI_sel = scans_tsv.EI_sel;
    else
        EI_sel = zeros(size(name,1),1);
    end
    
    
else
    scansnum = 1;
end

name{scansnum,1}                  = f;

% sleep period
id_sleep                          = strcmp(annotation_tsv.trial_type,'sleep');
durationsl_total = 0;
durationsl_rem = 0;
durationsl_nrem = 0;
annotsleep = find(id_sleep==1);

for i=1:sum(id_sleep)
    durationsl_total(i) = annotation_tsv.duration{annotsleep(i)};
    if strcmp(annotation_tsv.sub_type{annotsleep(i)},'REM')
        durationsl_rem(i) = annotation_tsv.duration{annotsleep(i)};
    elseif strcmp(annotation_tsv.sub_type{annotsleep(i)},'nREM')
        durationsl_nrem(i) = annotation_tsv.duration{annotsleep(i)};
    end
end

sleep_total(scansnum,1)           = sum(durationsl_total);
sleep_rem(scansnum,1)             = sum(durationsl_rem);
sleep_nrem(scansnum,1)            = sum(durationsl_nrem);

% motor period
id_motor                          = strcmp(annotation_tsv.trial_type,'motortask');
durationmt_total = 0;
annotmt = find(id_motor==1);

for i=1:sum(id_motor)
    durationmt_total(i) = annotation_tsv.duration{annotmt(i)};
end
motor(scansnum,1)           = sum(durationmt_total);

% language period
id_lang                          = strcmp(annotation_tsv.trial_type,'languagetask');
durationlang_total = 0;
annotlang = find(id_lang==1);

for i=1:sum(id_lang)
    durationlang_total(i) = annotation_tsv.duration{annotlang(i)};
end
language(scansnum,1)           = sum(durationlang_total);

% sensing task period
id_sens                          = strcmp(annotation_tsv.trial_type,'sensing task');
durationsens_total = 0;
annotsens = find(id_sens==1);

for i=1:sum(id_lang)
    durationsens_total(i) = annotation_tsv.duration{annotsens(i)};
end
sens(scansnum,1)           = sum(durationsens_total);


artefact(scansnum,1)              = sum(strcmp(annotation_tsv.trial_type,'artefact'));
seizure(scansnum,1)               = sum(strcmp(annotation_tsv.trial_type,'seizure'));
seizuresubclin(scansnum,1)        = sum(strcmp(annotation_tsv.sub_type,'subclin'));
seizureclin(scansnum,1)           = sum(strcmp(annotation_tsv.sub_type,'clin'));
spes(scansnum,1)                  = sum(strcmpi(annotation_tsv.sub_type,'spes'));
esm(scansnum,1)                   = sum(strcmpi(annotation_tsv.sub_type,'esm'));
slowesm(scansnum,1)               = sum(strcmpi(annotation_tsv.sub_type,'slowesm'));
sleepwaketransition(scansnum,1)   = sum(strcmp(annotation_tsv.trial_type,'sleep-wake transition'));
sws_sel(scansnum,1)               = sum(strcmp(annotation_tsv.trial_type,'sws selection'));
rem_sel(scansnum,1)               = sum(strcmp(annotation_tsv.trial_type,'rem selection'));
iiaw_sel(scansnum,1)               = sum(strcmp(annotation_tsv.trial_type,'iiaw selection'));
EI_sel(scansnum,1)               = sum(strcmp(annotation_tsv.trial_type,'EI selection'));

if metadata.incl_exist == 1
    format{scansnum,1}            = 'included';
else
    format{scansnum,1}            = 'not included';
end

scans_tsv  = table(name, format, artefact, sleep_total, sleep_rem, sleep_nrem,...
    sleepwaketransition, seizure, seizureclin, seizuresubclin, motor, spes, ...
    esm, slowesm, language, sens, sws_sel, rem_sel, iiaw_sel, EI_sel,...
    'VariableNames',{'name', 'format','artefact','sleep_total', 'sleep_rem',...
    'sleep_nrem','sleepwaketransition','seizure_total','seizure_clinical', ...
    'seizure_subclinical','motor','spes','esm','slowesm', 'language','sens',...
    'sws_se','rem_sel','iiaw_sel','EI_sel'});

if ~isempty(scans_tsv)
    write_tsv(filename, scans_tsv);
end

function write_participants_tsv(cfg,header)
[p,f] = fileparts(cfg(1).outputfile);
q = strsplit(p,'/');

% find session number
containsses = q{contains(q,'ses')};
sesnum = str2double(containsses(regexp(containsses,'\d')));

filename = ['/', q{2},'/' q{3},'/','participants.tsv'];

files = dir(['/',q{2},'/',q{3},'/']);
pat_exist = [];
if contains([files(:).name],'participants')
    % read existing scans-file
    participants_tsv = read_tsv(filename);
    
    if any(contains(participants_tsv.name,deblank(header.name))) % look whether the name is already in the participants-table
        partnum = find(contains(participants_tsv.name,deblank(header.name)) ==1 & participants_tsv.session==sesnum); %find patient number and session number
        pat_exist = 1;
    else % if participant is not yet in the table, the number is the last one plus one
        partnum = size(participants_tsv,1)+1;
    end
    
    name = participants_tsv.name;
    age = participants_tsv.age;
    session = participants_tsv.session;
else
    partnum = 1;
end

% set RESPect name and session number
name{partnum,1}   = deblank(header.name);
session(partnum,1) = sesnum;
% set age of RESPect patient (comparing with current participants-table)
if pat_exist == 1
    if age(partnum,1) == header.age && age(partnum,1) ~= 0 % if age in participants.tsv is not equal to 0  and equal to header.age
        age(partnum,1)    = header.age;
    elseif age(partnum,1) ~= 0 && header.age == 0 % if age is not equal to 0 (assumed to be correct)
        
    elseif age(partnum,1) == 0 && header.age ~= 0 % if age is equal to 0 and header.age is not (latter is assumed to be correct)
        age(partnum,1) = header.age;
    elseif age(partnum,1) ~= 0 && header.age ~= 0 && age(partnum,1) ~= header.age % if both ages are not 0 and conflicting, keep current age
        warning('ages between this file and other file are in conflict!')
    elseif age(partnum,1) == 0 && header.age == 0
        warning('age is 0 years... assumed to be incorrect!')
    end
else
    if header.age == 0
        warning('age is 0 years... assumed to be incorrect!')
    end
    age(partnum,1) = header.age;
end

% extract RESPect numbers from RESPect names
numname = zeros(size(name));
for i=1:size(name,1)
    numname(i) = str2double(name{i}(5:end));
end

% sorts table based on RESPect number and session number
[~,I] = sortrows([numname,session]);
 
name_sort = name(I);
age_sort = age(I);
session_sort = session(I);

% makes a table from name, session and age
participants_tsv  = table(name_sort, session_sort, age_sort,  ...
    'VariableNames',{'name','session', 'age'});

% save participants.tsv
if ~isempty(participants_tsv)
    write_tsv(filename, participants_tsv);
end


function eventsannots = add_spestrigger2annotation(stim,evname, eventsannots,annots_new, header,metadata,trigger)
fs = header.Rate_Min;
% ch2use_included = metadata.ch2use_included;
ch_label = metadata.ch_label;
maxchan = max(cellfun('length',ch_label(metadata.ch2use_included)));
if strcmpi(metadata.elec_info,'SEEG')
    stimcurdefault = 2;
elseif strcmpi(metadata.elec_info,'ECoG')
    stimcurdefault = 8;
end

% notification that stimcurr is unknown
if contains(lower(metadata.stimcurr),'unknown')
    note_desc = sprintf('Stimulation intensity is suggested to be %i mA but may differ when applied in eloquent tissue',stimcurdefault);
else
    note_desc = 'n/a';
end

if size(eventsannots.type,2) == 0
    cc = 1;
else
    cc = size(eventsannots.type,2)+1;
end

eventsannots.type{cc} = 'electrical_stimulation';
eventsannots.sub_type{cc} = evname;
eventsannots.stim_type{cc} = 'monophasic';
eventsannots.samp_start{cc} = trigger.pos(stim);
eventsannots.s_start{cc} = round(trigger.pos(stim)/fs,1); % time in seconds (1 decimal)

% stimulation site
[~,numannots]=max(1./(repmat(trigger.pos(stim),size(annots_new,1),1)-[annots_new{:,1}]')); %distance between triggerposition and nearbiest annotation (must be the stim channels then)

% stimpair
annotsplit = strsplit(annots_new{numannots,2},'_');
stimnames = regexp(lower(annotsplit{1}),'[a-z]*','match');
stimnums = regexp(lower(annotsplit{1}),'\d*','match');

% does this have 'ma'/'neg'/'bi'/'current is lower than expected' in the string? (respectively
% negative current, lower pulse current, biphasic instead of monophasic)
negannot = regexp(lower(annots_new{numannots,2}),'neg');
currannot = regexp(lower(annots_new{numannots,2}),'ma');
if ~isempty(currannot)
    annotsplit = strsplit(lower(annots_new{numannots,2}),'_');
    currsplit = strsplit(lower(annotsplit{2}),'ma');
    stimcurrstr = currsplit{1};
    stimcurr = str2double(stimcurrstr)/1000;
else
    stimcurr = stimcurdefault/1000;
end

biannot = regexp(lower(annots_new{numannots,2}),'bi');
low_expect = regexp(lower(annots_new{numannots,2}),'requested');
if ~isempty(low_expect)
    note = annots_new{numannots,2};
elseif numannots < size(annots_new,1)
    if contains(annots_new{numannots+1,2},'requested') % in newer patients, the sentence "Current is lower than requested" is added at the end of the stimulation
        note = annots_new{numannots+1,2};
    else
        note = note_desc;
    end
else
    note = note_desc;
end

if size(stimnums,2) ==2 && size(stimnames,2) == 2 % it is a stimulus pair
    [stimchan,stimnum] = findstimpair(stimnums,stimnames,ch_label);
    if find(stimnum == 0)
        error('one stimulation channel has not been found')
    end
elseif ~isempty(negannot) || ~isempty(currannot) || ~isempty(biannot) || ~isempty(low_expect) %% it is part of stim annotations
    n=1;
    while size(stimnums,2) <2 || size(stimnames,2) <2
        annotsplit = strsplit(annots_new{numannots-n,2},'_');
        stimnames = regexp(lower(annotsplit{1}),'[a-z]*','match');
        stimnums = regexp(lower(annotsplit{1}),'\d*','match');
        n = n+1;
    end
    numannots = numannots-n+1;
    
    annotsplit = strsplit(annots_new{numannots,2},'_');
    stimnames = regexp(lower(annotsplit{1}),'[a-z]*','match');
    stimnums = regexp(lower(annotsplit{1}),'\d*','match');
    
    [stimchan,stimnum] = findstimpair(stimnums,stimnames,ch_label);
    if find(stimnum == 0)
        error('one stimulation channel has not been found')
    end
    
else % if it is not part of stim annotations
    n=1;
    while size(stimnums,2) ~=2 || size(stimnames,2) ~= 2 || ~isempty(negannot) || ~isempty(currannot) || ~isempty(biannot) || ~isempty(low_expect) %% it is part of stim annotations
        annotsplit = strsplit(annots_new{numannots-n,2},'_');
        stimnames = regexp(lower(annotsplit{1}),'[a-z]*','match');
        stimnums = regexp(lower(annotsplit{1}),'\d*','match');
        n = n+1;
    end
    numannots = numannots-n+1;
    
    annotsplit = strsplit(annots_new{numannots,2},'_');
    stimnames = regexp(lower(annotsplit{1}),'[a-z]*','match');
    stimnums = regexp(lower(annotsplit{1}),'\d*','match');
    
    % does this have 'ma'/'neg'/'bi'/'current is lower than expected' in the string? (respectively
    % negative current, lower pulse current, biphasic instead of monophasic)
    negannot = regexp(lower(annots_new{numannots,2}),'neg');
    currannot = regexp(lower(annots_new{numannots,2}),'ma', 'once');
    if ~isempty(currannot)
        annotsplit = strsplit(lower(annots_new{numannots,2}),'_');
        currsplit = strsplit(lower(annotsplit{2}),'ma');
        stimcurrstr = currsplit{1};
        stimcurr = str2double(stimcurrstr)/1000;
    else
        stimcurr = stimcurdefault/1000;
    end
    
    biannot = regexp(lower(annots_new{numannots,2}),'bi');
    low_expect = regexp(lower(annots_new{numannots,2}),'requested', 'once');
    if ~isempty(low_expect)
        note = annots_new{numannots,2};
    elseif numannots < size(annots_new,1)
        if contains(annots_new{numannots+1,2},'requested') % in newer patients, the sentence "Current is lower than requested" is added at the end of the stimulation
            note = annots_new{numannots+1,2};
        else
            note = note_desc;
        end
    else
        note = note_desc;
    end
    if size(stimnums,2) <2 || size(stimnames,2) <2 % if it is a stim annotation, but no stimpair
        n=1;
        while size(stimnums,2) <2 || size(stimnames,2) <2
            annotsplit = strsplit(annots_new{numannots-n,2},'_');
            stimnames = regexp(lower(annotsplit{1}),'[a-z]*','match');
            stimnums = regexp(lower(annotsplit{1}),'\d*','match');
            n = n+1;
        end
        numannots = numannots-n+1;
        
        annotsplit = strsplit(annots_new{numannots,2},'_');
        stimnames = regexp(lower(annotsplit{1}),'[a-z]*','match');
        stimnums = regexp(lower(annotsplit{1}),'\d*','match');
    end
    
    [stimchan,stimnum] = findstimpair(stimnums,stimnames,ch_label);
    if find(stimnum == 0)
        error('one stimulation channel has not been found')
    end
    
end

if ~isempty(negannot)
    stimchantemp = stimchan{1};
    stimchan{1} = stimchan{2};
    stimchan{2} = stimchantemp;
    stimnumtemp = stimnum(1);
    stimnum(1) = stimnum(2);
    stimnum(2) = stimnumtemp;
end

if ~isempty(biannot)
    eventsannots.stim_type{cc} = 'biphasic';
end

eventsannots.site_name{cc} = [stimchan{1} '-' stimchan{2}];
eventsannots.site_channum{cc} = num2str([stimnum(1), stimnum(2)]);
eventsannots.duration{cc} = 1/1000;
eventsannots.s_end{cc} = 'n/a';
eventsannots.samp_end{cc} = 'n/a';
eventsannots.ch_name_on{cc} = 'n/a';
eventsannots.ch_name_off{cc} = 'n/a';
eventsannots.stim_cur{cc} = stimcurr;
eventsannots.notes{cc} = note;

function eventsannots = add_esmtrigger2annotation(stim,eventsannots,annots_new, header,metadata,trigger)  
fs = header.Rate_Min;
ch_label = metadata.ch_label;

if size(eventsannots.type,2) == 0
    cc = 1;
else
    cc = size(eventsannots.type,2)+1;
end

eventsannots.type{cc} = 'electrical_stimulation';
eventsannots.sub_type{cc} = 'ESM';
eventsannots.stim_type{cc} = 'biphasic';
eventsannots.samp_start{cc} = trigger.pos(stim);
eventsannots.s_start{cc} = round(trigger.pos(stim)/fs,1); % time in seconds (1 decimal)

% stimulation site
[~,numannots]=max(1./(repmat(trigger.pos(stim),size(annots_new,1),1)-[annots_new{:,1}]')); %distance between triggerposition and nearbiest annotation (must be the stim channels then)

annotsplit = strsplit(annots_new{numannots,2},char(32));
stimchans = strsplit(annotsplit{1},'-');
stimcur = cellfun(@(x) contains(x,{'ma'}),lower(annotsplit))==1;
stimfreq = cellfun(@(x) contains(x,{'hz'}),lower(annotsplit))==1;
stimwidth = cellfun(@(x) contains(x,{'sec'}),lower(annotsplit))==1;

if size(stimchans,2) == 1
    stimnumber = regexp(lower(stimchans{1}),'\d*','match');
    stimname = regexp(lower(stimchans{1}),'[a-z]*','match');
    for j=1:size(stimnumber,2)
        stimnum(j) = find(strcmpi([stimname{j} stimnumber{j}],ch_label)==1);
    end
    eventsannots.site_name{cc} = [ch_label{stimnum(1)}, '-', ch_label{stimnum(2)}];
    eventsannots.site_channum{cc} = num2str([stimnum(1), stimnum(2)]);
elseif size(stimchans,2) == 2
    for j=1:size(stimchans,2)
        stimnum(j) = find(strcmpi(stimchans{j},ch_label)==1);
    end
    eventsannots.site_name{cc} = annotsplit{1};
    eventsannots.site_channum{cc} = num2str([stimnum(1), stimnum(2)]);
end


% current
cur = str2double(strsplit(lower(annotsplit{stimcur}),'ma'));

% duration
if ismember(65533,double(annotsplit{stimwidth})) % there is a weird letter, translated from microseconds
   dur = str2double(strsplit(annotsplit{stimwidth},char(65533)))/1000; % to change it to miliseconds
end

note = sprintf('Frequency is %s, width is %1.3f sec, but duration of total stimulus is unclear',annotsplit{stimfreq},dur(1)/1000);

eventsannots.duration{cc} = dur(1)/1000;
eventsannots.s_end{cc} = 'unknown';
eventsannots.samp_end{cc} = 'unknown';
eventsannots.ch_name_on{cc} = 'n/a';
eventsannots.ch_name_off{cc} = 'n/a';
eventsannots.stim_cur{cc} = cur(1)/1000;
eventsannots.notes{cc} = note;


function [stimchan,stimnum] = findstimpair(stimnums,stimchans,ch_label)

stimchan = cell(size(stimnums));
stimnum = zeros(size(stimnums));
for j=1:size(stimnums,2)
    test1 = sprintf('%s%d',stimchans{j},str2double(stimnums{j}));
    test2 = sprintf('%s0%d',stimchans{j},str2double(stimnums{j}));
    if sum(strcmpi(test1,ch_label))==1
        stimchan{j} = ch_label{strcmpi(test1,ch_label)};
        stimnum(j) = find(strcmpi(test1,ch_label)==1);
    elseif sum(strcmpi(test2,ch_label))==1
        stimchan{j} = ch_label{strcmpi(test2,ch_label)};
        stimnum(j) = find(strcmpi(test2,ch_label)==1);
    end
end

%% miscellaneous functions from data2bids.m of fieldtrip

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function tsv = read_tsv(filename)
ft_info('reading %s\n', filename);
tsv = readtable(filename, 'Delimiter', 'tab', 'FileType', 'text', 'ReadVariableNames', true);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function write_tsv(filename, tsv)
ft_info('writing %s\n', filename);
writetable(tsv, filename, 'Delimiter', 'tab', 'FileType', 'text');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function str = truefalse(bool)
if bool
    str = 'true';
else
    str = 'false';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function s = remove_empty(s)
fn = fieldnames(s);
fn = fn(structfun(@isempty, s));
s = removefields(s, fn);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function x = mergevector(x, y)
assert(isequal(size(x), size(y)));
for i=1:numel(x)
    if isnumeric(x) && isnumeric(y) && isnan(x(i)) && ~isnan(y(i))
        x(i) = y(i);
    end
    if iscell(x) && iscell(y) && isempty(x{i}) && ~isempty(y{i})
        x{i} = y{i};
    end
    if iscell(x) && isnumeric(y) && isempty(x{i}) && ~isnan(y{i})
        x{i} = y(i);
    end
end

%% check if the configuration struct contains all the required fields
function check_input(cfg,key) % used on line 21

if (isa(cfg, 'struct'))
    
    fn = fieldnames(cfg);
    if ~any(strcmp(key, fn))
        
        error('Provide the configuration struct with all the fields example: cfg.proj_dir  cfg.filename  error: %s missing ', key);
    end
    
else
    error('Provide the configuration struct with all the fields example: cfg.proj_dir  cfg.filename');
end
