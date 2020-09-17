function write_scans_tsv(cfg,metadata,annotation_tsv,fscans_name,fieeg_json_name)

f = replace(fieeg_json_name,'.json','');

for j=1:size(cfg(1).ses_dir,2)
    filename = fullfile(cfg(1).ses_dir{j},fscans_name);
    
    files = dir(cfg(1).ses_dir{j});
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
        
        if any(contains(scans_tsv.Properties.VariableNames,'rec2stim')) % this was added later, so not allf iles have this
            rec2stim = scans_tsv.rec2stim;
        else
            rec2stim = zeros(size(name,1),1);
        end
        
         
        if any(contains(scans_tsv.Properties.VariableNames,'chocs')) % this was added later, so not allf iles have this
            chocs = scans_tsv.chocs; 
        else
            chocs = zeros(size(name,1),1);
        end
        
        
    else
        scansnum = 1;
    end
    
    name{scansnum,1}                  = f;
    
    % sleep period
    id_sleep                          = strcmp(annotation_tsv.trial_type,'sleep');
    annotsleep = find(id_sleep==1);
    
    % pre-allocation
    durationsl_total = zeros(sum(id_sleep),1); durationsl_rem = zeros(sum(id_sleep),1); durationsl_nrem = zeros(sum(id_sleep),1);
    for i=1:sum(id_sleep)
        durationsl_total(i,1) = annotation_tsv.duration{annotsleep(i)};
        if strcmp(annotation_tsv.sub_type{annotsleep(i)},'REM')
            durationsl_rem(i,1) = annotation_tsv.duration{annotsleep(i)};
        elseif strcmp(annotation_tsv.sub_type{annotsleep(i)},'nREM')
            durationsl_nrem(i,1) = annotation_tsv.duration{annotsleep(i)};
        end
    end
    
    sleep_total(scansnum,1)           = sum(durationsl_total);
    sleep_rem(scansnum,1)             = sum(durationsl_rem);
    sleep_nrem(scansnum,1)            = sum(durationsl_nrem);
    
    % motor period
    id_motor                          = strcmp(annotation_tsv.trial_type,'motortask');
    annotmt = find(id_motor==1);
    
    durationmt_total = zeros(sum(id_motor),1);
    for i=1:sum(id_motor)
        durationmt_total(i,1) = annotation_tsv.duration{annotmt(i)};
    end
    motor(scansnum,1)           = sum(durationmt_total);
    
    % language period
    id_lang                          = strcmp(annotation_tsv.trial_type,'languagetask');
    annotlang = find(id_lang==1);
    durationlang_total = zeros(sum(id_lang),1);
    
    for i=1:sum(id_lang)
        durationlang_total(i,1) = annotation_tsv.duration{annotlang(i)};
    end
    language(scansnum,1)           = sum(durationlang_total);
    
    % sensing task period
    id_sens                          = strcmp(annotation_tsv.trial_type,'sensing task');
    durationsens_total = zeros(sum(id_sens),1);
    annotsens = find(id_sens==1);
    
    for i=1:sum(id_sens)
        durationsens_total(i) = annotation_tsv.duration{annotsens(i)};
    end
    sens(scansnum,1)           = sum(durationsens_total);
    
    artefact(scansnum,1)              = sum(strcmp(annotation_tsv.trial_type,'artefact'));
    seizure(scansnum,1)               = sum(strcmp(annotation_tsv.trial_type,'seizure'));
    seizuresubclin(scansnum,1)        = sum(strcmp(annotation_tsv.sub_type,'subclin'));
    seizureclin(scansnum,1)           = sum(strcmp(annotation_tsv.sub_type,'clin'));
    spes(scansnum,1)                  = sum(contains(lower(annotation_tsv.sub_type),'spes'));
    rec2stim(scansnum,1)              = sum(strcmpi(annotation_tsv.sub_type,'rec2stim'));
    esm(scansnum,1)                   = sum(strcmpi(annotation_tsv.sub_type,'esm'));
    chocs(scansnum,1)                 = sum(contains(lower(annotation_tsv.sub_type),'chocs'));
    sleepwaketransition(scansnum,1)   = sum(strcmp(annotation_tsv.trial_type,'sleep-wake transition'));
    sws_sel(scansnum,1)               = sum(strcmp(annotation_tsv.trial_type,'sws selection'));
    rem_sel(scansnum,1)               = sum(strcmp(annotation_tsv.trial_type,'rem selection'));
    iiaw_sel(scansnum,1)              = sum(strcmp(annotation_tsv.trial_type,'iiaw selection'));
    EI_sel(scansnum,1)                = sum(strcmp(annotation_tsv.trial_type,'EI selection'));
    
    if metadata.incl_exist == 1
        format{scansnum,1}            = 'included';
    else
        format{scansnum,1}            = 'not included';
    end
    
    scans_tsv  = table(name, format, artefact, sleep_total, sleep_rem, sleep_nrem,...
        sleepwaketransition, seizure, seizureclin, seizuresubclin, motor, spes, rec2stim, ...
        esm, chocs, language, sens, sws_sel, rem_sel, iiaw_sel, EI_sel,...
        'VariableNames',{'name', 'format','artefact','sleep_total', 'sleep_rem',...
        'sleep_nrem','sleepwaketransition','seizure_total','seizure_clinical', ...
        'seizure_subclinical','motor','spes','rec2stim','esm','chocs', 'language','sens',...
        'sws_se','rem_sel','iiaw_sel','EI_sel'});
    
    if ~isempty(scans_tsv)
        
        write_tsv(filename, scans_tsv);
    end
end