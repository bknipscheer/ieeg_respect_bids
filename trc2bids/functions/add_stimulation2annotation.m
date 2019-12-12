function eventsannots = add_stimulation2annotation(metadata,header,annots_new,eventsannots)

trigger = metadata.trigger;

% FIX THIS: this code now assumes that only 1 type of stimulation is applied in
% one trc-file. 

% determine period of stimulation
stimperiod = metadata.stimulation{1}.pos;
stimtype = metadata.stimulation{1}.type;

if ~isempty(trigger.pos)
    idx_start = find(trigger.val >1000); % with cortical stimulation, triggers are added automatically with a number >1000
    for i=1:numel(idx_start)
        
        eventsannots = add_stimtrigger2annotation(idx_start(i),stimtype,eventsannots,annots_new,header,metadata,trigger);
                
    end
end

% in older ECoGs, there are no triggers, but stimulation NEED TO BE FIXED
if sum(cellfun(@(x) contains(x,{'No trigger'}),annots_new(:,2)))>0
    
    eventsannots = add_stimnotrigger2annotation(cfg,metadata,annots_new,fevents_name);   
    
end