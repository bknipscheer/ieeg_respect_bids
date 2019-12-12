% Example  anonymize the TRC file and fill in the RESPect number
% author Dorien van Blooijs
% date: 2019

config_anonymization

%% choose the eeg-file and anonymize this file/these files
files = dir(cfg.proj_dirinput);
for i=1:size(files,1)
    if contains(files(i).name,'EEG_')
        filename = files(i).name;
        pathname = cfg.proj_dirinput;
        
        fileName = [pathname, filename];
        % anonymize the TRC file
        [status,msg] = anonymized_asRecorded(fileName,respName)
        
    end
end