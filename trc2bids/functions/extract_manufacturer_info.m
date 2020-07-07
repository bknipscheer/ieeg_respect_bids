% extract manufacturer info
function ch_manufacturer = extract_manufacturer_info(metadata)

ch = metadata.ch;
ch2use_included = metadata.ch2use_included;
ch_manufacturer = cell(size(ch));
[ch_manufacturer{:}] = deal('n/a');

if contains(lower(metadata.electrode_manufacturer),'adtech;') || contains(lower(metadata.electrode_manufacturer),'pmt;') || contains(lower(metadata.electrode_manufacturer),'dixi;')
    C = strsplit(metadata.electrode_manufacturer,{';'});
    
    for i=1:size(C,2)
       if contains(lower(C{i}),'adtech')
           pat = 'AdTech';
       elseif contains(lower(C{i}),'dixi')
           pat = 'Dixi';
       elseif contains(lower(C{i}),'pmt')
           pat = 'PMT';
       elseif contains(C{i},'[')
           
           ch_subset_str=parse_ch_subset(C{i},ch);
           ch2use_temp = zeros(size(ch));
           for chan = 1:size(ch_subset_str,1)
               ch2use_temp(:,chan) = strcmpi(ch,ch_subset_str{chan});
           end
           
           ch2use = sum(ch2use_temp,2);
           idx_ch2use = boolean(ch2use);
           [ch_manufacturer{idx_ch2use}] = deal(pat);
       end
    end
    
else
    warning('Electrode model is not annotated. AdTech is now filled in in electrodes_tsv.')
    [ch_manufacturer{ch2use_included}] = deal('AdTech');
end

           
           
 
    
    
   