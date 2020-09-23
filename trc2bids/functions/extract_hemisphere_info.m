% extract manufacturer info
function ch_hemi = extract_hemisphere_info(metadata)

ch = metadata.ch;
ch2use_included = metadata.ch2use_included;
ch_hemi = cell(size(ch));
[ch_hemi{:}] = deal('n/a');

if contains(lower(metadata.hemisphere),'right') || contains(lower(metadata.hemisphere),'left')
    C = strsplit(metadata.hemisphere,{';'});
    C = C(~cellfun(@isempty,C));
    
    if size(C,2) == 1
        if contains(lower(C{1}),'right')
            side = 'R';
        elseif contains(lower(C{1}),'left')
            side = 'L';
        end
        [ch_hemi{ch2use_included}] = deal(side);

    else
        
        for i=1:size(C,2)
            if contains(lower(C{i}),'right')
                side = 'R';
            elseif contains(lower(C{i}),'left')
                side = 'L';
            elseif contains(C{i},'[')
                
                ch_subset_str=parse_ch_subset(C{i},ch);
                ch2use_temp = zeros(size(ch));
                for chan = 1:size(ch_subset_str,1)
                    ch2use_temp(:,chan) = strcmpi(ch,ch_subset_str{chan});
                end
                
                ch2use = sum(ch2use_temp,2);
                idx_ch2use = boolean(ch2use);
                [ch_hemi{idx_ch2use}] = deal(side);
            end
        end
    end
    
else
    warning('Hemisphere is not annotated. "right" is now filled in in electrodes_tsv.')
    [ch_hemi{ch2use_included}] = deal('right');
end

           
           
 
    
    
   