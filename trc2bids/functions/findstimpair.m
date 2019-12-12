function [stimchan,stimnum] = findstimpair(stimnums,stimchans,ch_label)

stimchan = cell(size(stimnums));
stimnum = zeros(size(stimnums));
for j=1:size(stimnums,2)
    test1 = sprintf('%s%d',stimchans{j},str2double(stimnums{j}));
    test2 = sprintf('%s0%d',stimchans{j},str2double(stimnums{j}));
    if sum(strcmpi(test1,ch_label))==1
        stimchan{j} = ch_label{strcmpi(test1,ch_label)};
        stimnum(j) = find(strcmpi(test1,ch_label)==1);
    elseif sum(strcmpi(test2,ch_label))==1
        stimchan{j} = ch_label{strcmpi(test2,ch_label)};
        stimnum(j) = find(strcmpi(test2,ch_label)==1);
    end
end
