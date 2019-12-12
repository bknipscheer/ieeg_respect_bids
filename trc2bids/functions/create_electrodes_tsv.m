function create_electrodes_tsv(cfg,metadata,header,felectrodes_name)

temp.electrodes = [];

temp.electrodes.name             = ft_getopt(temp.electrodes, 'name'               , nan);
temp.electrodes.x                = ft_getopt(temp.electrodes, 'x'                  , nan);
temp.electrodes.y                = ft_getopt(temp.electrodes, 'y'                  , nan);
temp.electrodes.z                = ft_getopt(temp.electrodes, 'z'                  , nan);
temp.electrodes.size             = ft_getopt(temp.electrodes, 'size'               , nan);
temp.electrodes.group            = ft_getopt(temp.electrodes, 'group'              , nan);
temp.electrodes.material         = ft_getopt(temp.electrodes, 'material'           , nan);
temp.electrodes.manufacturer     = ft_getopt(temp.electrodes, 'manufacturer'       , nan);
temp.electrodes.silicon          = ft_getopt(temp.electrodes, 'silicon'            , nan);
temp.electrodes.soz              = ft_getopt(temp.electrodes, 'soz'                , nan);
temp.electrodes.ra               = ft_getopt(temp.electrodes, 'ra'                 , nan);
temp.electrodes.edge             = ft_getopt(temp.electrodes, 'edge'               , nan);

fn = {'name' 'x' 'y' 'z' 'size' 'group' 'material' 'manufacturer' 'silicon' 'soz' 'ra' 'edge'};
for i=1:numel(fn)
    if numel(temp.electrodes.(fn{i}))==1
        temp.electrodes.(fn{i}) = repmat(temp.electrodes.(fn{i}), header.Num_Chan, 1);
    end
end

name                                = mergevector({header.elec(:).Name}', temp.electrodes.name)                                   ;
x                                         = repmat({0},header.Num_Chan,1)                                                            ;
y                                         = repmat({0},header.Num_Chan,1)                                                            ;
z                                         = repmat({0},header.Num_Chan,1)                                                            ;
e_size                                    = repmat({'n/a'},header.Num_Chan,1)                                                          ; %TODO ask
group                               = extract_group_info(metadata)                                                             ;
material                                  = repmat({'n/a'},header.Num_Chan,1)                                                          ; %TODO ask
manufacturer                              = repmat({'n/a'},header.Num_Chan,1)                                                          ; %TODO ask
silicon                                   = repmat({'no'},header.Num_Chan,1)                                                          ; %TODO ask
soz                                       = repmat({'no'},header.Num_Chan,1)                                                          ; %TODO ask
resected                                  = repmat({'no'},header.Num_Chan,1)                                                          ; %TODO ask
edge                                      = repmat({'no'},header.Num_Chan,1)                                                          ; %TODO ask

if contains(lower(metadata.format_info),'seeg')
    screw                                     = repmat({'no'},header.Num_Chan,1)                                                          ; %TODO ask
    csf                                       = repmat({'no'},header.Num_Chan,1)                                                          ; %TODO ask
    lesion                                    = repmat({'no'},header.Num_Chan,1)                                                          ; %TODO ask
    gm                                        = repmat({'no'},header.Num_Chan,1)                                                          ; %TODO ask
    wm                                        = repmat({'no'},header.Num_Chan,1)                                                          ; %TODO ask
    hipp                                      = repmat({'no'},header.Num_Chan,1)                                                          ; %TODO ask
    amyg                                      = repmat({'no'},header.Num_Chan,1)                                                          ; %TODO ask
    gliosis                                   = repmat({'no'},header.Num_Chan,1)                                                          ; %TODO ask
end

if(any(metadata.ch2use_included))
    [e_size{metadata.ch2use_included}]        = deal('4.2')                                                                                ;
end

if(any(metadata.ch2use_included))
    [material{metadata.ch2use_included}]      = deal('Platinum')                                                                                ;
end

if(any(metadata.ch2use_included))
    [manufacturer{metadata.ch2use_included}]  = deal('AdTech')                                                                                ;
end

if(any(metadata.ch2use_silicon))
    [silicon{metadata.ch2use_silicon}]  = deal('yes')                                                                                ;
end

if(any(metadata.ch2use_soz))
    [soz{metadata.ch2use_soz}]  = deal('yes')                                                                                ;
end

if(any(metadata.ch2use_resected))
    [resected{metadata.ch2use_resected}]  = deal('yes')                                                                                ;
end

if(any(metadata.ch2use_edge))
    [edge{metadata.ch2use_edge}]  = deal('yes')                                                                                ;
end

if contains(lower(metadata.format_info),'seeg')
    if(any(metadata.ch2use_lesion))
        [lesion{metadata.ch2use_lesion}]  = deal('yes')                                                                                ;
    end
    if(any(metadata.ch2use_wm))
        [wm{metadata.ch2use_wm}]  = deal('yes')                                                                                ;
    end
    if(any(metadata.ch2use_gm))
        [gm{metadata.ch2use_gm}]  = deal('yes')                                                                                ;
    end
    if(any(metadata.ch2use_csf))
        [csf{metadata.ch2use_csf}]  = deal('yes')                                                                                ;
    end
    if(any(metadata.ch2use_screw))
        [screw{metadata.ch2use_screw}]  = deal('yes')                                                                                ;
    end
    if(any(metadata.ch2use_hipp))
        [hipp{metadata.ch2use_hipp}]  = deal('yes')                                                                                ;
    end
    
    if(any(metadata.ch2use_amyg))
        [amyg{metadata.ch2use_amyg}]  = deal('yes')                                                                                ;
    end
    if(any(metadata.ch2use_gliosis))
        [gliosis{metadata.ch2use_gliosis}]  = deal('yes')                                                                                ;
    end
    
    electrodes_tsv                            = table(name, x , y, z, e_size, ...
        group, material, manufacturer, silicon, soz, resected, edge ,...
        screw,csf,wm,gm,hipp,amyg,lesion,gliosis,...
        'VariableNames',{'name', 'x', 'y', 'z', 'size', ...
        'group', 'material', 'manufacturer', 'silicon' 'soz','resected','edge',...
        'screw','csf','whitematter','graymatter','hippocampus','amygdala','lesion','gliosis'})     ;
    
else
    electrodes_tsv                            = table(name, x , y, z, e_size, group, material, manufacturer, silicon, soz, resected, edge ,...
        'VariableNames',{'name', 'x', 'y', 'z', 'size', 'group', 'material', 'manufacturer', 'silicon' 'soz','resected','edge'})     ;
end

if ~isempty(electrodes_tsv)
    for i=1:size(cfg(1).ieeg_dir,2)
        filename = fullfile(cfg(1).ieeg_dir{i},felectrodes_name);
        
        if isfile(filename)
            
            cc_elec_old = readtable(filename,'FileType','text','Delimiter','\t');
            
            struct1 = table2struct(cc_elec_old);
            struct2 = table2struct(electrodes_tsv);
            if ~isequal(struct1,struct2)
                fprintf('%s exists!\n',filename)
                n=1;
                while isfile(filename)
                    nameminelec = strsplit(filename,'electrodes');
                    filename = [nameminelec{1} 'electrodes_' num2str(n) '.tsv'];
                    n=n+1;
                end
            end
        end
        write_tsv(filename, electrodes_tsv);
    end
end

end
