% Example  anonymize the TRC file and fill in the RESPect number
% author Dorien van Blooijs
% date: 2019

%% set paths 
% clone https://github.com/dvanblooijs/ieeg_respect_bids
addpath('git_rep/ieeg_respect_bids/anonymization/')
addpath('git_rep/ieeg_respect_bids/example/')

clear

% respect-folder on bulkstorage
cfg.proj_dirinput = '/home/dorien/Desktop/temp_ecog/';

% type patient decoding-number
tempName = input('Respect(name (e.g. [RESP0733]): ','s');

% check whether RESP-number is entered correctly
if strcmp(tempName,'') && ~isempty(respName)
    
elseif contains(tempName,'RESP')
    respName = tempName;
else
    error('Respect name is not correct')
end

% choose the eeg-file and anonymize this file/these files
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