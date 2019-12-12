function eventsannots = add_stimnotrigger2annotation(cfg,metadata,annots_new,fevents_name)

% UPDATE THIS FILE WHEN DATA WITHOUT TRIGGERS IS USED!

% load ECoG
dataName = [cfg(1).ieeg_dir{1},'/', replace(fevents_name,'_events.tsv','_ieeg.vhdr')];
data_raw = ft_read_data(dataName,'dataformat','brainvision_eeg');
data = data_raw(metadata.ch2use_included,:);
% ch_label = {ch_label{ch2use_included}};

% annotation sample
numnotrigger = find(cellfun(@(x) contains(x,{'No trigger'}),annots_new(:,2))==1);
for i=1:size(numnotrigger,1) % for each 'No trigger'-annotation
    if annots_new{numnotrigger(i),1} > spesstart && annots_new{numnotrigger(i),1} < spesend % if 'No trigger' is within SPESperiod
        if ~isempty(trigger.pos) % if triggers are present, then period ends with next trigger
            periodend = trigger.pos(find(annots_new{numnotrigger(i),1}<trigger.pos,1,'first'))-round(0.25*fs); %-0.25s takes care of an annotation prior to the next trigger
            sampend = periodend;
        else %otherwise, period ends with end of SPESperiod
            periodend = spesend;
            sampend = [];
        end
        % annotations within the specified period
        numannots = find([annots_new{:,1}]<periodend & [annots_new{:,1}] > annots_new{numnotrigger(i),1} ==1);
        
        for j=1:size(numannots,2)
            % find stimulation pair and other stimulus settings
            annotsplit = strsplit(annots_new{numannots(j),2},'_');
            stimnumber = regexp(lower(annotsplit{1}),'\d*','match');
            stimname = regexp(lower(annotsplit{1}),'[a-z]*','match');
            %                    currannot = regexp(lower(annots_new{numannots(j),2}),'ma', 'once');
            
            if size(annotsplit,2)>1
                currsplit = strsplit(lower(annotsplit{2}),'ma');
                stimcurrstr = currsplit{1};
                stimcurr = str2double(stimcurrstr)/1000;
            else
                stimcurr = stimcurdefault/1000;
            end
            
            if size(stimnumber,2) == 2 && size(stimname,2)==2 % when both have size=2, then it should be a stimulus pair
                
                [~,stimnum] = findstimpair(stimnumber,stimname,ch_label);
                
                % find stimulation samples for each annotation
                sampstart = annots_new{numannots(j),1};
                if j==1 && isempty(sampend) % if there are no triggers present
                    sampend = annots_new{numannots(j)+1,1};
                elseif numannots(j) == size(annots_new,2) % if triggers are not present in the last stimulation pair
                    sampend = spesend;
                else % if there are triggers present, than sampend is before the next trigger
                    sampend = periodend;
                end
                
                samplocs = findtrigger(data,stimnum,fs, sampstart, sampend)+annots_new{numannots(j),1};
                
                eventssize = size(eventsannots.type,2) ;
                for cc = eventssize+1:eventssize+size(samplocs,2)
                    eventsannots.type{cc} = 'electrical_stimulation';
                    eventsannots.sub_type{cc} = 'SPES';
                    eventsannots.stim_type{cc} = 'monophasic';
                    eventsannots.samp_start{cc} = samplocs(cc-eventssize);
                    eventsannots.s_start{cc} = round(samplocs(cc-eventssize)/fs,1); % time in seconds (1 decimal)
                    
                    eventsannots.site_name{cc} = [ch_label{stimnum(1)}, '-', ch_label{stimnum(2)}];
                    eventsannots.site_channum{cc} = num2str([stimnum(1), stimnum(2)]);
                    eventsannots.duration{cc} = 1/1000;
                    eventsannots.s_end{cc} = 'n/a';
                    eventsannots.samp_end{cc} = 'n/a';
                    eventsannots.ch_name_on{cc} = 'n/a';
                    eventsannots.ch_name_off{cc} = 'n/a';
                    eventsannots.stim_cur{cc} = stimcurr;
                    eventsannots.notes{cc} = note;
                    eventsannots.freq{cc} = freq;
                    
                end 
            end 
        end
    end
end

if strcmpi(metadata.elec_info,'SEEG')
    stimcurdefault = 2;
elseif strcmpi(metadata.elec_info,'ECoG')
    stimcurdefault = 8;
end

if contains(lower(metadata.stimcurr),'unknown')
    note = sprintf('Stimulation intensity is suggested to be %i mA but may differ when applied in eloquent tissue, triggers were added automatically in Matlab',stimcurdefault);
else
    note = 'Triggers were added automatically in Matlab';
end


% sort to put the no-triggers in the right order of stimulation
SPEStype        = strcmp(sub_type,'SPES');
noSPEStype      = ~SPEStype;
[~,I]           = sort([samp_start{SPEStype}]) ;
I               = I + find(SPEStype==1,1,'first')-1;
s_start         = {s_start{noSPEStype},s_start{I}};
s_end           = {s_end{noSPEStype},s_end{I}};
duration        = {duration{noSPEStype},duration{I}};
type            = {type{noSPEStype},type{I}};
sub_type        = {sub_type{noSPEStype},sub_type{I}};
ch_name_on      = {ch_name_on{noSPEStype},ch_name_on{I}};
ch_name_off     = {ch_name_off{noSPEStype},ch_name_off{I}};
samp_start      = {samp_start{noSPEStype},samp_start{I}};
samp_end        = {samp_end{noSPEStype},samp_end{I}};
stim_type       = {stim_type{noSPEStype},stim_type{I}};
site_name       = {site_name{noSPEStype},site_name{I}};
site_channum    = {site_channum{noSPEStype},site_channum{I}};
stim_cur        = {stim_cur{noSPEStype},stim_cur{I}};
notes           = {notes{noSPEStype},notes{I}};
freq            = {freq{noSPEStype},freq{I}};

end