function convertTRC2brainvision(cfg,ieeg_dir, fieeg_name)

filename = cfg(1).filename;

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

temp = [];
temp.dataset                     = filename;
temp.continuous = 'yes';
data2write = ft_preprocessing(temp);

temp = [];
temp.outputfile                  = fileVHDR{1};

temp.mri.writesidecar       = 'no';
temp.meg.writesidecar        = 'no';
temp.eeg.writesidecar        = 'no';
temp.ieeg.writesidecar       = 'no';
temp.channels.writesidecar   = 'no';
temp.events.writesidecar     = 'no';

% write .vhdr, .eeg, .vmrk
data2bids(temp, data2write)

% to each output-file in fileVHDRcopy, the data should be copied
if ~isempty(fileVHDRcopy{1})
    for i=1:size(fileVHDRcopy,2)
        copyfile(temp.outputfile,fileVHDRcopy{i});
        fprintf('Copy to %s\n', fileVHDRcopy{i})
        file_eeginp = replace(temp.outputfile,'.vhdr','.eeg');
        file_eegoutp = replace(fileVHDRcopy{i},'.vhdr','.eeg');
        copyfile(file_eeginp,file_eegoutp);
        fprintf('Copy to %s\n',file_eegoutp)
        file_vmrkinp = replace(temp.outputfile,'.vhdr','.vmrk');
        file_vmrkoutp = replace(fileVHDRcopy{i},'.vhdr','.vmrk');
        copyfile(file_vmrkinp,file_vmrkoutp);
        fprintf('Copy to %s\n',file_vmrkoutp)
    end
end


end