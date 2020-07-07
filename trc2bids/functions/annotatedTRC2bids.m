%  Convert annotated (see annotation scheme in docs) micromed file (.TRC) to Brain Imaging Data Structure (BIDS)
%  it generate all the required directory structure and files
%
%  cfg.proj_dir - directory name where to store the files
%  cfg.filename - name of the micromed file to convert
%
%  output structure with some information about the subject
%  output.subjName - name of the subject
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

function [status,msg,metadata,annots] = annotatedTRC2bids(cfg)

try
    % check whether all required fields are in cfg
    check_input(cfg,'proj_dirinput');
    check_input(cfg,'proj_diroutput');
    check_input(cfg,'filename');
    
    msg = '';
    filename  = cfg(1).filename;
    
    proj_diroutput = {cfg.proj_diroutput};
    proj_diroutput = proj_diroutput(~cellfun(@isempty, proj_diroutput));
    
    % obtain information from the header of the trc-file
    [header,data,data_time,trigger,annots] = read_TRC_HDR_DATA_TRIGS_ANNOTS(filename);
    
    if(isempty(header) || isempty(data) || isempty(data_time) || isempty(annots))
        error('TRC reading failed')  ;
    end
    
    ch_label = deblank({header.elec.Name}');
    sub_label = strcat('sub-',deblank(header.name));
    
    [status,msg,metadata] = extract_metadata_from_annotations(header,annots,ch_label,trigger,sub_label,cfg);
    
    if(status==0)
        
        ses_label     = strcat('ses-',deblank(metadata.ses_name),' ','');
        run_label     = strcat('run-',deblank(metadata.run_name),header.hour,header.min,' ','');
        task_label    = strcat('task-',deblank(metadata.task_name),' ','');
        
        % make directories
        sub_dir       = fullfile(proj_diroutput,sub_label);
        ses_dir       = fullfile(proj_diroutput,sub_label,ses_label);
        ieeg_dir      = fullfile(proj_diroutput,sub_label,ses_label,'ieeg');
        ieeg_file     = strcat(sub_label,'_',ses_label,'_',task_label,'_',run_label);
        cfg(1).ieeg_dir = ieeg_dir;
        
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
        
        % make names of files to be constructed
        fieeg_name = strcat(sub_label,'_',ses_label,'_',task_label,'_',run_label,'_','ieeg.TRC');
        fieeg_json_name = strcat(sub_label,'_',ses_label,'_',task_label,'_',run_label,'_','ieeg','.json');
        fchannels_name = strcat(sub_label,'_',ses_label,'_',task_label,'_',run_label,'_','channels','.tsv');
        felectrodes_name = strcat(sub_label,'_',ses_label,'_','electrodes','.tsv');
        fevents_name = strcat(sub_label,'_',ses_label,'_',task_label,'_',run_label,'_','events','.tsv');
        fscans_name = strcat(sub_label,'_',ses_label,'_scans','.tsv');
        
        %% create Brainvision format from TRC
        
        convertTRC2brainvision(cfg,ieeg_dir, fieeg_name)
        
        %% create json sidecar for ieeg file
        
        cfg = create_jsonsidecar(cfg,metadata,header,fieeg_json_name);
        
        %% create _channels.tsv
        
        create_channels_tsv(cfg,metadata,header,fchannels_name)
        
        %% create _electrodes.tsv
        
        if metadata.incl_exist == 1
            create_electrodes_tsv(cfg,metadata,header,felectrodes_name)
        end
        
        %% write annotations of the TRC
        
        annotations_tsv = write_annotations_tsv(cfg,metadata,header,annots,fevents_name);
        
        %% write scans-file
        
        write_scans_tsv(cfg,metadata,annotations_tsv,fscans_name,fieeg_json_name)
        
        %% write participants-file
        
        write_participants_tsv(cfg,header,metadata)
               
        %% write dataset descriptor
        
        for i=1:size(proj_diroutput,2)
            create_datasetDesc(proj_diroutput{i},sub_label)
        end

        %% write event descriptor
        
        for i=1:size(proj_diroutput,2)
            create_eventDesc(proj_diroutput{i})
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





%% FUNCTIONS %%

% make directory
function mydirMaker(dirname)

if exist(dirname, 'dir')
    warning('%s exist already',dirname)
else
    mkdir(dirname)
end


% check if the configuration struct contains all the required fields
function check_input(cfg,key) % used on line 21

if (isa(cfg, 'struct'))
    
    fn = fieldnames(cfg);
    if ~any(strcmp(key, fn))
        
        error('Provide the configuration struct with all the fields example: cfg.proj_dir  cfg.filename  error: %s missing ', key);
    end
    
else
    error('Provide the configuration struct with all the fields example: cfg.proj_dir  cfg.filename');
end
