%% visualisation of specific electrodes
% This script renders a brian surface with Destrieux maps
% Can potentially add electrodes on top
% dhermes & jvanderaar & dvanblooijs 2019, UMC Utrecht


function check_atlas_elec_MRI_pre_post(cfg,tb_elecs)

if any(contains(fieldnames(cfg),'transparency'))
   transparency = cfg.transparency;
else
    transparency = 1;
end

if any(strcmp(fieldnames(cfg),'atlas'))
    
else
       cfg.atlas = 'DKT';
end

% pick a viewing angle:
v_dirs = [90 0; 270 0]; %[90 0;90 -60;270 -60;0 0;270 0; 270 60]; %zij, onder, .., voor, zij, zijboven

% gifti file name:
dataGiiName = replace(fullfile(cfg.surface_directory,...
    [cfg.sub_labels{:} '_' cfg.ses_label '_T1w_pial.' cfg.hemisphere '.surf.gii']),'_post','');
dataGiiName_post = fullfile(cfg.surface_directory,...
    [cfg.sub_labels{:} '_' cfg.ses_label '_T1w_pial.' cfg.hemisphere '.surf.gii']);
% load gifti:
g = gifti(dataGiiName);
g_post = gifti(dataGiiName_post);

% surface labels
if strcmp(cfg.atlas,'DKT')
    surface_labels_name = fullfile(cfg.freesurfer_directory,'label',...
        [cfg.hemisphere 'h.aparc.DKTatlas.annot']);
elseif strcmp(cfg.atlas,'Destrieux')
    surface_labels_name = fullfile(cfg.freesurfer_directory,'label',...
        [cfg.hemisphere 'h.aparc.a2009s.annot']);
end
% surface_labels = MRIread(surface_labels_name);
[~, label, colortable] = read_annotation(surface_labels_name);
vert_label = label; % these labels are strange and do not go from 1:76, but need to be mapped to the colortable
% mapping labels to colortable
for kk = 1:size(colortable.table,1) % 76 are labels
    vert_label(label==colortable.table(kk,5)) = kk;
end

% make a colormap for the labels
cmap = colortable.table(:,1:3)./256;

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
    
    subplot('position', [0.05 0.25 0.45 0.7])
    
    if strcmp(cfg.view_atlas,'yes')
        ecog_RenderGiftiLabels(g,vert_label,cmap,colortable.struct_names)
    else
        ecog_RenderGifti(g,transparency) % render
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
            
            % add all electrodes with  dots
            ccep_el_add(els,[0 1 0],20) % [electrodes, MarkerColor, MarkerSize]
        end
    end
    title(sprintf('%s MRI & electrodes pre-resection',cfg(1).sub_labels{:}))
    
    subplot('position', [0.55 0.25 0.45 0.7])
    
    if strcmp(cfg.view_atlas,'yes')
        ecog_RenderGiftiLabels(g_post,vert_label,cmap,colortable.struct_names)
    else
        ecog_RenderGifti(g_post,transparency) % render
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
            
            % add all electrodes with dots
            ccep_el_add(els,[0 1 0],20) % [electrodes, MarkerColor, MarkerSize]
        end
    end
    title(sprintf('%s MRI and electrodes post-resection',cfg(1).sub_labels{:}))
    
    set(gcf,'PaperPositionMode','auto')
    
    
end
end