%% description of this code
% This scripts defines the patient study name and converts ieeg data to
% BIDS with accompanying meta-files

% Dorien van Blooijs, Willemiek Zweiphenning 2019

clear 

%% patient characteristics

cfg.sub_labels = {['sub-' input('Patient number (RESPXXXX)/(REC2StimXX)/(PRIOSXX): ','s')]};

%% set paths
cfg = setLocalDataPath(cfg);

%% 1a) TRC to bids - run all files in patient-folder
files = dir(cfg(1).pathname);
runall = struct;

if size(files,1)<1
    error('Pathname is wrong, no files found')
end

% run all files within your input directory
for i=1:size(files,1) 
    runall(i).file = files(i).name;
    if contains(files(i).name,'EEG_')
        
        cfg(1).filename = [cfg(1).pathname,files(i).name];
        
        pathsplit = strsplit(cfg(1).pathname,{'/'});
        patient = pathsplit{end-1};
        filesplit = strsplit(files(i).name,{'_','.TRC'});
        file = filesplit{end-1};
        
        fprintf('Running %s, writing EEG: %s to BIDS \n', patient,file)
        [runall(i).status,runall(i).msg,runall(i).metadata,runall(i).annots] = annotatedTRC2bids(cfg);
    end
end

if any([runall(:).status])
    disp('All runs are done, but some still have errors. Fix them manually!')
else
    disp('All runs are completed')
end

%% 1b) run files which gave errors again 

for i=1:size(runall,2)
    
    if runall(i).status ==1 
        cfg(1).filename = [pathname,runall(i).file];
        
        pathsplit = strsplit(pathname,{'/'});
        patient = pathsplit{end-1};
        filesplit = strsplit(runall(i).file,{'_','.TRC'});
        file = filesplit{end-1};
        
        fprintf('Running %s, writing EEG: %s to BIDS \n', patient,file)
        [runall(i).status,runall(i).msg,runall(i).metadata,runall(i).annots] = annotatedTRC2bids(cfg);

    end
end

if any([runall(:).status])
    disp('All runs are done, but some still have errors')
else
    disp('All runs are completed')
end

%% 2) run one single file instead of all files within the input directory
      
files = dir(pathname);
eegfiles = {files(contains({files(:).name},'EEG')==1).name};
string = [repmat('%s, ',1,size(eegfiles,2)-1),'%s'];

if size(files,1) <1
    error('Pathname does not contain any files')
else
    fileinput = input(sprintf(['Select one of these files [',string,']: \n'],eegfiles{:}),'s');
end

cfg(1).filename = [pathname,fileinput];

pathsplit = strsplit(pathname,{'/'});
patient = pathsplit{end-1};
filesplit = strsplit(fileinput,{'_','.TRC'});
file = filesplit{end-1};

fprintf('Running %s, writing EEG: %s to BIDS \n', patient,file)
[status,msg,metadata,annots] = annotatedTRC2bids(cfg);

if status
    disp('Run is done, but still had an error')
else
    disp('Run is completed')
end

