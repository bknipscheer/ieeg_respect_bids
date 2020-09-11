%% visualisation of specific electrodes
% This script renders a brian surface with Destrieux maps
% Can potentially add electrodes on top
% dhermes & jvanderaar & dvanblooijs 2019, UMC Utrecht


function check_atlas_elec_MRI(cfg,tb_elecs)


if contains(fieldnames(cfg),'transparency')
    transparency = cfg.transparency;
else
    transparency = 1;
end

% pick a viewing angle:
v_dirs = [90 0; 270 0]; %[90 0;90 -60;270 -60;0 0;270 0; 270 60]; %zij, onder, .., voor, zij, zijboven

for i=1:size(cfg.hemisphere,2)
    
    % gifti file name:
    dataGiiName = fullfile(cfg.surface_directory,...
        [cfg.sub_labels{:} '_' cfg.ses_label '_T1w_pial.' cfg.hemisphere{i} '.surf.gii']);
    % load gifti:
    g.(cfg.hemisphere{i}) = gifti(dataGiiName);
    
    % surface labels
    if strcmp(cfg.atlas,'DKT')
        surface_labels_name = fullfile(cfg.freesurfer_directory,'label',...
            [cfg.hemisphere{i} 'h.aparc.DKTatlas.annot']);
    elseif strcmp(cfg.atlas,'Destrieux')
        surface_labels_name = fullfile(cfg.freesurfer_directory,'label',...
            [cfg.hemisphere{i} 'h.aparc.a2009s.annot']);
    end
    % surface_labels = MRIread(surface_labels_name);
    [~, label, colortable] = read_annotation(surface_labels_name);
    vert_label.(cfg.hemisphere{i}) = label; % these labels are strange and do not go from 1:76, but need to be mapped to the colortable
    % mapping labels to colortable
    for kk = 1:size(colortable.table,1) % 76 are labels
        vert_label.(cfg.hemisphere{i})(label==colortable.table(kk,5)) = kk;
    end
    
    % make a colormap for the labels
    cmap = colortable.table(:,1:3)./256;
end

% electrode locations name:
if isempty(tb_elecs)
    dataLocName = dir(fullfile(cfg.ieeg_directory,...
        [cfg.sub_labels{:},'_',cfg.ses_label '_electrodes.tsv']));
    dataLocName = fullfile(dataLocName(1).folder,dataLocName(1).name);
    % load electrode locations
    tb_elecs = readtable(dataLocName,'FileType','text','Delimiter','\t','TreatAsEmpty','n/a');
end
log_elec_incl = ~strcmp(tb_elecs.group,'other');
tb_elecs = tb_elecs(log_elec_incl,:);
if iscell(tb_elecs.x)
    elecmatrix = [cell2mat(tb_elecs.x) cell2mat(tb_elecs.y) cell2mat(tb_elecs.z)];
else
    elecmatrix = [tb_elecs.x tb_elecs.y tb_elecs.z];
end

%% figure with rendering for different viewing angles
for k = 1:size(v_dirs,1) % loop across viewing angles
    v_d = v_dirs(k,:);
    
    figure('units','normalized','position',[0.01 0.01 0.9 0.9],'color',[1 1 1]);
    
    subplot('position', [0.05 0.25 0.9 0.7])
    
    for i=1:size(cfg.hemisphere,2)
        
        if i == size(cfg.hemisphere,2)
            setLight = 1;
        else
            setLight = 0;
        end
        
        if strcmp(cfg.view_atlas,'yes')
            ecog_RenderGiftiLabels(g.(cfg.hemisphere{i}),vert_label.(cfg.hemisphere{i}),cmap,colortable.struct_names,setLight)
        else
            ecog_RenderGifti(g.(cfg.hemisphere{i}),transparency,setLight) % render
        end
    end
    ecog_ViewLight(v_d(1),v_d(2)) % change viewing angle
    
    if strcmp(cfg.view_elec,'yes')
        
        if strcmp(cfg.elec_offset,'yes')
            % make sure electrodes pop out
            a_offset = 0.1*max(abs(elecmatrix(:,1)))*[cosd(v_d(1)-90)*cosd(v_d(2)) sind(v_d(1)-90)*cosd(v_d(2)) sind(v_d(2))];
            els = elecmatrix+repmat(a_offset,size(elecmatrix,1),1);
        else
            els = elecmatrix;
        end
        
        % add electrode numbers
        if strcmp(cfg.show_labels,'yes')
            ecog_Label(els,tb_elecs.name,30,12) % [electrodes, electrode labels, MarkerSize, FontSize]
            
            % add all electrodes with black dots
            ccep_el_add(els,[0.1 0.1 0.1],20) % [electrodes, MarkerColor, MarkerSize]
        end
    end
    
    set(gcf,'PaperPositionMode','auto')
    
    
end
end