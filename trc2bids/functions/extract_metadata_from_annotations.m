%% extract all metadata needed for bid structure

% annots - annotations of the trc file
% ch     - channel labels of all channels in the trc file


function [status,msg,metadata]=extract_metadata_from_annotations(header,annots,ch,trigger,patName,cfg) 
try
    status=0;
    metadata=[];
    
    %% Check the compulsory fields
    
    metadata.sub_label = patName;
    metadata.ch = ch;
    
    %% ---------- SESSION ----------
    ses_idx=cellfun(@(x) contains(x,{'Session'}),annots(:,2));
    if(sum(ses_idx)~=1)
        warning('Missing Session annotation (example "Session;1"), so ses-1 is set')
        metadata.ses_name='1';
    else
        str2parse=annots{ses_idx,2};
        %metadata.ses_name=strsplit(str2parse,'n');
        C=strsplit(str2parse,';');
        metadata.ses_name=C{2};
    end
    
    %% ---------- RUN ----------
    run_idx=cellfun(@(x) contains(x,{'Run'}),annots(:,2));
    if(sum(run_idx)~=1)
        status=1;
        error('Missing run annotation (example "day1") or more than one run annotation')
    end
    str2parse=annots{run_idx,2};
    C=strsplit(str2parse,';');
    D = strsplit(C{2},'y');
    if str2double(D{2}) <10
        metadata.run_name=['0' deblank(D{2})];
    else
        metadata.run_name=deblank(D{2});
    end
    
    %% ---------- TASK ----------
    task_idx=cellfun(@(x) contains(x,{'Task'}),annots(:,2));
    if(sum(task_idx)~=1)
        status=1;
        error('Missing task annotation (example "Task;SPES") or more than one task annotation')
    end
    str2parse=annots{task_idx,2};
    C=strsplit(str2parse,';');
    metadata.task_name=deblank(C{2});
      
    if contains(metadata.task_name,'SPES')
        % Stimcurr unknown (in SPES)
        stimcur_idx=cellfun(@(x) contains(x,{'Stimcurr'}),annots(:,2));
        
        if(sum(stimcur_idx)~=1)
            metadata.stimcurr = [];
            error('Missing (un)known stimulation current')
        else
            str2parse=annots{stimcur_idx,2};
            C=strsplit(str2parse,';');
            metadata.stimcurr=C{2};
        end
    end
    
    %% ---------- FORMAT ----------
    format_idx=cellfun(@(x) contains(x,{'Format'}),annots(:,2));
    if(sum(format_idx)<1)
        file = dir(fullfile(cfg(2).proj_dirinput,patName,['ses-',metadata.ses_name],'ieeg','*ieeg.json'));
        if ~isempty(file)
            ieeg_json = jsondecode(fileread([file(1).folder '/' file(1).name]) );
            metadata.format_info = ieeg_json.iEEGElectrodeGroups;
        else
            error('Missing Format annotation, and no other json from other ECoG with format annotation found')
        end
    else
        metadata = look_for_format(metadata,format_idx,annots);
    end
    
    %% ---------- INCLUDED ---------- 
    included_idx=cellfun(@(x) contains(x,{'Included'}),annots(:,2));
    metadata.ch2use_included= false(size(ch));
    if (sum(included_idx))
        metadata.ch2use_included=single_annotation(annots,'Included',ch);
        fprintf('File had Included-annotation\n')
        metadata.incl_exist = 1;

    else % if "Included" is not annotated in the ECoG, there should be a previous ECoG with annoted "Included"
        metadata.incl_exist = 0;
        
        files = dir(fullfile(cfg(2).proj_dirinput,patName,['ses-',metadata.ses_name],'ieeg',[patName, '_ses-',metadata.ses_name,'_electrodes.tsv']));
        if ~isempty(files)

            metadata = load_chanInfo(cfg,metadata,files);
            
        else
            error('There is no ECoG with annotated Included')
        end
        
    end
    
    %% ---------- ELECTRODE MODEL ---------- 
    elecmodel_idx=cellfun(@(x) contains(x,{'Elec_model'}),annots(:,2));
    if(sum(elecmodel_idx)<1)
        file = dir(fullfile(cfg(2).proj_dirinput,patName,['ses-',metadata.ses_name],'ieeg','*ieeg.json'));
        if ~isempty(file)
            ieeg_json = jsondecode(fileread([file(1).folder '/' file(1).name]) );
            metadata.electrode_manufacturer = ieeg_json.ElectrodeManufacturer;
        else
            warning('Missing Electrode_manufacturer annotation, and no other json from other ECoG with electrode_manufacturer annotation found')
            metadata.electrode_manufacturer = 'supposed to be AdTech, but needs to be checked'; 
        end
    else
        metadata = look_for_electrode_manufacturer(metadata,elecmodel_idx,annots);
    end
    
     %% ---------- GENDER ----------
    gender_idx=cellfun(@(x) contains(x,{'Gender'}),annots(:,2));
    if(sum(gender_idx)~=1)
        warning('Gender is missing')
        metadata.gender = 'unknown';
    else
        str2parse=annots{gender_idx,2};
        C=strsplit(str2parse,';');
        metadata.gender=C{2};
    end
    
    %% Look for bad channels
    metadata.ch2use_bad=single_annotation(annots,'Bad;',ch); % without the semicolon, the bad_HF channels are also included in Bad
    
    %% Look for bad channels in high frequency band
    badhf_idx = cellfun(@(x) contains(x,{'Bad_HF'}),annots(:,2));
    metadata.ch2use_badhf= false(size(ch));
    if(sum(badhf_idx))
        metadata.ch2use_badhf=single_annotation(annots,'Bad_HF',ch);
        if any(contains(annots(:,2),'NB BadHF annotated in avg'))
            metadata.ch2use_badhf_note = 'NB BadHF annotated in avg';
        else
             metadata.ch2use_badhf_note = '';
        end
    end
    
    %% Look for silicon
    silicon_idx=cellfun(@(x) contains(x,{'Silicon'}),annots(:,2));
    metadata.ch2use_silicon= false(size(ch));
    if(sum(silicon_idx))
        metadata.ch2use_silicon=single_annotation(annots,'Silicon',ch);
    elseif metadata.incl_exist == 0
        % done in load_chanInfo
    end
    
    %% look for resected channels
    resected_idx = cellfun(@(x) contains(x,{'RA'}),annots(:,2));
    metadata.ch2use_resected= false(size(ch));
    if(sum(resected_idx))
        metadata.ch2use_resected=single_annotation(annots,'RA',ch);
    elseif metadata.incl_exist == 0
        % done in load_chanInfo
    end
    
    %% look for edge channels
    edge_idx = cellfun(@(x) contains(x,{'Edge'}),annots(:,2));
    metadata.ch2use_edge= false(size(ch));
    if(sum(edge_idx))
        metadata.ch2use_edge=single_annotation(annots,'Edge',ch);
    elseif metadata.incl_exist == 0
        % done in load_chanInfo
    end
    
    %% look for SOZ channels
    soz_idx = cellfun(@(x) contains(x,{'SOZ'}),annots(:,2));
    metadata.ch2use_soz= false(size(ch));
    if(sum(soz_idx))
        metadata.ch2use_soz=single_annotation(annots,'SOZ',ch);
    elseif metadata.incl_exist == 0
        % done in load_chanInfo
    end
    
    %% only in seeg channels
    if contains(metadata.format_info,'seeg') 
        
        % look for screw channels - only in seeg
        screw_idx = cellfun(@(x) contains(x,{'Screw'}),annots(:,2));
        metadata.ch2use_screw= false(size(ch));
        if(sum(screw_idx))
            metadata.ch2use_screw=single_annotation(annots,'Screw',ch);
        elseif metadata.incl_exist == 0
            % done in load_chanInfo
        end
        
        % look for white matter channels - only in seeg
        wm_idx = cellfun(@(x) contains(x,{'WM'}),annots(:,2));
        metadata.ch2use_wm= false(size(ch));
        if(sum(wm_idx))
            metadata.ch2use_wm=single_annotation(annots,'WM',ch);
        elseif metadata.incl_exist == 0
            % done in load_chanInfo            
        end
        
        % look for gray matter channels - only in seeg
        gm_idx = cellfun(@(x) contains(x,{'GM'}),annots(:,2));
        metadata.ch2use_gm= false(size(ch));
        if(sum(gm_idx))
            metadata.ch2use_gm=single_annotation(annots,'GM',ch);
        elseif metadata.incl_exist == 0
            % done in load_chanInfo            
        end
        
        % look for CSF channels - only in seeg
        csf_idx = cellfun(@(x) contains(x,{'CSF'}),annots(:,2));
        metadata.ch2use_csf= false(size(ch));
        if(sum(csf_idx))
            metadata.ch2use_csf=single_annotation(annots,'CSF',ch);
        elseif metadata.incl_exist == 0
            % done in load_chanInfo            
        end
        
        % look for amygdala channels - only in seeg
        amyg_idx = cellfun(@(x) contains(x,{'Amyg'}),annots(:,2));
        metadata.ch2use_amyg= false(size(ch));
        if(sum(amyg_idx))
            metadata.ch2use_amyg=single_annotation(annots,'Amyg',ch);
        elseif metadata.incl_exist == 0
            % done in load_chanInfo            
        end
        
        % look for hippocampal channels - only in seeg
        hipp_idx = cellfun(@(x) contains(x,{'Hipp'}),annots(:,2));
        metadata.ch2use_hipp= false(size(ch));
        if(sum(hipp_idx))
            metadata.ch2use_hipp=single_annotation(annots,'Hipp',ch);
        elseif metadata.incl_exist == 0
            % done in load_chanInfo            
        end
        
        % look for lesion channels - only in seeg
        lesion_idx = cellfun(@(x) contains(x,{'Lesion'}),annots(:,2));
        metadata.ch2use_lesion= false(size(ch));
        if(sum(lesion_idx))
            metadata.ch2use_lesion=single_annotation(annots,'Lesion',ch);
        elseif metadata.incl_exist == 0
            % done in load_chanInfo            
        end
        
        % look for gliosis channels - only in seeg
        gliosis_idx = cellfun(@(x) contains(x,{'Glio'}),annots(:,2));
        metadata.ch2use_gliosis= false(size(ch));
        if(sum(gliosis_idx))
            metadata.ch2use_lesion=single_annotation(annots,'Glio',ch);
        elseif metadata.incl_exist == 0
            % done in load_chanInfo            
        end
    end
    
    %% Look for artefacts cECoG
    metadata.artefacts=look_for_annotation_start_stop(annots,'Art_on','Art_off',ch,header);
    
    if any(contains(annots(:,2),'NB artefacts annotated in avg'))
        metadata.artefacts_note = 'NB artefacts annotated in avg';
    end
    
    %% Look for sleep data
    metadata.sleep=look_for_annotation_start_stop(annots,'Sl_on','Sl_off',ch,header);
    
    %% Look for data between rest and sleep
    metadata.slaw_trans = look_for_annotation_start_stop(annots,'Slawtrans_on','Slawtrans_off',ch,header);
    
    %% Look for data eyes opened
    metadata.eyes_open = look_for_annotation_start_stop(annots,'Eyes_open','Eyes_close',ch,header);
    
    %% Look for seizures
    metadata.seizure=look_for_annotation_start_stop(annots,'Sz_on','Sz_off',ch,header);
    
    %% Look for period of stimulation (SPES,ESM,CHOCS,TRENI,REC2Stim)
    metadata.stimulation=look_for_annotation_start_stop(annots,'Stim_on','Stim_off',ch,header);
       
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
   
    %% Hemisphere where electrodes are placed
    hemisphere_idx = cellfun(@(x) contains(x,{'Hemisphere'}),annots(:,2));
    if (sum(hemisphere_idx) ~= 1)
        if exist('ieeg_json','var') % if ieeg_json is loaded in determining Format
            if isfield(ieeg_json,'iEEGPlacementScheme') % if iEEGPlacementScheme is in ieeg_json
                metadata.hemisphere = ieeg_json.iEEGPlacementScheme;
            else
                warning('Hemisphere where electrodes are implanted is not mentioned')
                metadata.hemisphere='unknown';
            end
        else
            warning('Hemisphere where electrodes are implanted is not mentioned')
            metadata.hemisphere='unknown';
        end
    else
        str2parse=annots{hemisphere_idx,2};
        C=strsplit(str2parse,{'; ',';'});
        metadata.hemisphere=C{2};
        if size(C,2) >2
            warning('Annotation in "Hemisphere" might be incorrect')
        end
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
