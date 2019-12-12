function cfg = create_jsonsidecar(cfg,metadata,header,fieeg_json_name)

ieeg_json    = [];

task_label    = strcat('task-',replace(deblank(metadata.task_name),' ',''));
if contains(lower(task_label),'spes') || contains(lower(task_label),'rec2stim')
    task_desc = 'No task, electrical stimulation is performed. Patient is resting with eyes open/closed. The latter is not always specified.';
elseif contains(lower(task_label),'rest')
    task_desc = 'Patient is resting with eyes open/closed. The latter is not always specified.';
elseif contains(lower(task_label),'sleep')
    task_desc = 'Patient is sleeping.';
elseif contains(lower(task_label),'slawtrans')
    task_desc = 'Patient is trying to fall asleep or is waking up.';
elseif contains(lower(task_label),'motor')
    task_desc = 'Patient is doing a motor task.';
elseif contains(lower(task_label),'esm')
    task_desc = 'Electrical stimulation mapping is performed to delineate functional areas.';
elseif contains(lower(task_label),'sens')
    task_desc = 'Patient is doing a sensing task.';
elseif contains(lower(task_label),'Language')
    task_desc = 'Patient is doing a language task.';
else
    task_desc = 'Not specified';
    warning('Task description is not specified!')
end

ieeg_json.TaskName                    = task_label;
ieeg_json.TaskDescription             = task_desc;
ieeg_json.Manufacturer                = 'Micromed';
ieeg_json.ManufacturersModelName      = header.acquisition_eq;%sprintf('Acqui.eq:%i  File_type:%i',header.acquisition_eq,header.file_type);
ieeg_json.DeviceSerialNumber          = '';
ieeg_json.SoftwareVersions            = num2str(header.Header_Type);
ieeg_json.SoftwareFilters             = 'n/a';
ieeg_json.InstitutionName             = 'University Medical Center Utrecht';
ieeg_json.InstitutionalDepartmentName = 'Clinical Neurophysiology Department';
ieeg_json.InstitutionAddress          = 'Heidelberglaan 100, 3584 CX Utrecht';
if strfind(ieeg_json.ManufacturersModelName,'LTM') ~=0
    ieeg_json.HardwareFilters.HighpassFilter.CutoffFrequency =             0.15;
    if header.Rate_Min/2.21 < 468
        ieeg_json.HardwareFilters.LowpassFilter.CutoffFrequency = round(header.Rate_Min/2.21);
    else
        ieeg_json.HardwareFilters.LowpassFilter.CutoffFrequency  =             468;
    end
elseif strcmp(ieeg_json.ManufacturersModelName,'SD128')
    ieeg_json.HardwareFilters.HighpassFilter.CutoffFrequency =             0.15;
    ieeg_json.HardwareFilters.LowpassFilter.CutoffFrequency  =             round(header.Rate_Min/3.81);
elseif strcmp(ieeg_json.ManufacturersModelName,'SD64')
    ieeg_json.HardwareFilters.HighpassFilter.CutoffFrequency =             0.15;
    ieeg_json.HardwareFilters.LowpassFilter.CutoffFrequency  =             round(header.Rate_Min/3.81);
end

cfg(1).HardwareFilters.HighpassFilter.CutoffFrequency = ieeg_json.HardwareFilters.HighpassFilter.CutoffFrequency;
cfg(1).HardwareFilters.LowpassFilter.CutoffFrequency = ieeg_json.HardwareFilters.LowpassFilter.CutoffFrequency;


%% IEEG inherited fields used
if contains(lower(metadata.format_info),'ecog')
    ieeg_json.ECOGChannelCount             = sum(metadata.ch2use_included);
    ieeg_json.SEEGChannelCount             = 0;
elseif contains(lower(metadata.format_info),'seeg')
    ieeg_json.ECOGChannelCount             = 0;
    ieeg_json.SEEGChannelCount             = sum(metadata.ch2use_included);
end

%ieeg_json.EEGChannelCount              =
%ieeg_json.EOGChannelCount              =
ieeg_json.ECGChannelCount              = sum(~cellfun(@isempty,regexpi(metadata.ch_label,'ECG')));
%ieeg_json.EMGChannelCount              =
ieeg_json.RecordingDuration            = header.Num_Samples/header.Rate_Min;
ieeg_json.RecordingType                = 'continuous';
ieeg_json.EpochLength                  = 0;


%% IEEG specific fields
ieeg_json.SamplingFrequency            = header.Rate_Min;
ieeg_json.PowerLineFrequency           = 50;
ieeg_json.iEEGReference                = 'probably mastoid';
ieeg_json.ElectrodeManufacturer        = 'AdTech';
ieeg_json.iEEGPlacementScheme          = metadata.hemisphere;
ieeg_json.iEEGElectrodeGroups          = metadata.format_info;
if strfind(ieeg_json.TaskName,'SPES') ~=0
    ieeg_json.ElectricalStimulation        = 'true';
end

for i=1:size(cfg(1).ieeg_dir,2)
filename = fullfile(cfg(1).ieeg_dir{i},fieeg_json_name);

if ~isempty(filename)
    if isfile(filename)
        existing = read_json(filename);
    else
        existing = [];
    end
    write_json(filename, mergeconfig(existing, ieeg_json))
    %     json_options.indent = ' ';
    %     jsonwrite(filename, mergeconfig(existing, ieeg_json), json_options)
end

end
