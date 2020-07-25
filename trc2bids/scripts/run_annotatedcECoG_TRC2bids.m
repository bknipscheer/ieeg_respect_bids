%% description of this code
% This scripts defines the patient study name and converts ieeg data to
% BIDS with accompanying meta-files

% Dorien van Blooijs, Willemiek Zweiphenning 2019


%% patient characteristics

clear cfg
cfg.sub_labels = {['sub-' input('Patient number (RESPXXXX)/(REC2StimXX)/(PRIOSXX): ','s')]};
cfg.mode = 'bidsconversion';

% set paths
cfg = setLocalDataPath(cfg);

%% run all patients in database
pats = dir(cfg(1).proj_dirinput);
runpat = struct;
for pat = 1:size(pats,1)
    
    if contains(pats(pat).name,'PAT')
        cfg(1).pathname = [fullfile(cfg(1).proj_dirinput,pats(pat).name),'/'];
        
        files = dir(cfg(1).pathname);
        
        if size(files,1)<1
            error('Pathname is wrong, no files found')
        end
        
        % run all files within your input directory
        for i=1:size(files,1)
            runpat(pat).runall(i).file = files(i).name;
            if contains(files(i).name,'EEG_')
                
                cfg(1).filename = [cfg(1).pathname,files(i).name];
                
                pathsplit = strsplit(cfg(1).pathname,{'/'});
                patient = pathsplit{end-1};
                filesplit = strsplit(files(i).name,{'_','.TRC'});
                file = filesplit{end-1};
                
                fprintf('Running %s, writing EEG: %s to BIDS \n', patient,file)
                [runpat(pat).runall(i).status,runpat(pat).runall(i).msg,runpat(pat).runall(i).metadata,runpat(pat).runall(i).annots] = annotatedTRC2bids(cfg);
            
            end
        end
        
        if any([runpat(pat).runall(:).status])
            disp('All runs are done, but some still have errors. Fix them manually!')
        else
            disp('All runs are completed')
        end
    end
end

% check which patients do not run without errors
if contains(fieldnames(runpat),'status')
    runpat = rmfield(runpat, 'status');
end

for i=1:size(runpat,2)
    
    if ~isempty(runpat(i).runall)
        
        if any(vertcat(runpat(i).runall(:).status) == 1)
            
            runpat(i).status = 1;
            
            
        end
    end
end
    
sum([runpat(:).status])

%% 2a) TRC to bids - run all files in patient-folder

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

%% 2b) run files which gave errors again

for i=1:size(runall,2)
    
    if runall(i).status ==1
        cfg(1).filename = [cfg(1).pathname,runall(i).file];
        
        pathsplit = strsplit(cfg(1).pathname,{'/'});
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

%% 3) run one single file instead of all files within the input directory

files = dir(cfg(1).pathname);
eegfiles = {files(contains({files(:).name},'EEG')==1).name};
string = [repmat('%s, ',1,size(eegfiles,2)-1),'%s'];

if size(files,1) <1
    error('Pathname does not contain any files')
else
    fileinput = input(sprintf(['Select one of these files [',string,']: \n'],eegfiles{:}),'s');
end

cfg(1).filename = [cfg(1).pathname,fileinput];

pathsplit = strsplit(cfg(1).pathname,{'/'});
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

