%% config_anonymization

%% set paths 
% clone https://github.com/dvanblooijs/ieeg_respect_bids
addpath(genpath('git_rep/ieeg_respect_bids/anonymization/'))

clear

% respect-folder on bulkstorage
cfg.proj_dirinput = '/home/dorien/Desktop/temp_ecog/';

%% type patient decoding-number
tempName = input('RESPectname/REC2Stimname (e.g. [RESP0733] / [REC2Stim01] / [PRIOS01]): ','s');

% check whether RESP-number is entered correctly
if strcmp(tempName,'') && ~isempty(respName)
    
elseif contains(tempName,'RESP') || contains(tempName,'REC2Stim') || contains(tempName,'PRIOS')
    respName = tempName;
else
    error('RESPect/REC2Stim name is not correct')
end