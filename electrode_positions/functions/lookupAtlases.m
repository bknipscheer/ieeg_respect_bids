function [electrodes_tableWithlabels, destrieux_labels, DKT_labels] = lookupAtlases(cfg,tb_elecs)
%
% This function looks up DKT and Destrieux atlas labels for electrodes and
% writes them to the _electrodes.tsv file
% for every electrode, it looks up the gifti vertices within mm_distance
% and extracts the labels from the corresponding atlas. It takes the most
% represented brain area within 3 mm.
%
% Inputs:
% dataRootPath
% subj:             subject(s) number
% ses:              sessios(s) number
% freesurfer_dir:   directory with freesurfer output + Benson & Kastner maps
% hemi_small:       hemisphere to look up labels (smaller case: l or h)
%
% Output:
% electrode_tsv table with labels added for every electrode
%
% D Hermes, G Castegnaro, J van der Aar and D van Blooijs, UMC Utrecht, 2019

% load gifti file and electrodes.tsv file
% g = gifti(fullfile(dataRootPath,'derivatives','surfaces',['sub-' subj],...
%     ['sub-' subj '_T1w_pial.' hemi_cap '.surf.gii']));
g = gifti(fullfile(cfg.surface_directory,...
    [cfg.sub_labels{:}, '_', cfg.ses_label, '_T1w_pial.' cfg.hemisphere{1} '.surf.gii']));

% electrodes_tsv = [cfg.ieeg_directory,  '_electrodes.tsv'];

% define name output file
% output_file = 'electrodes.tsv';

% load the electrodes.tsv file:
% loc_info = readtable(electrodes_tsv,'FileType','text','Delimiter','\t','TreatAsEmpty',{'N/A','n/a'});
% elecmatrix = [loc_info.x loc_info.y loc_info.z];
elecmatrix = [tb_elecs.x tb_elecs.y tb_elecs.z];

%%%% LOAD SURFACE LABELS FOR ALL THE MAPS

% DESTRIEUX ATLAS
% load Destrieux map
surface_labels_Destrieux = fullfile(cfg.freesurfer_directory,'label',...
    sprintf('%sh.aparc.a2009s.annot',cfg.hemisphere{1}));
[~, vert_label_Destrieux, colortable_Destrieux] = read_annotation(surface_labels_Destrieux);

% REORDER FREESURFER DESTRIEUX LABELS TO DESTRIEUX LABELS FROM Destrieux et al. NeuroImage2011
% destrieux_labels and attach the corresponding number
% Medial wall is in the middle of the labels but should be at the end
destrieux_labels = colortable_Destrieux.struct_names;
idx = contains(destrieux_labels,'Medial_wall');
destrieux_labels = destrieux_labels(:).';
destrieux_labels = [destrieux_labels(~idx), destrieux_labels(idx)];
destrieux_labels = destrieux_labels(:);

% unknown = 0, S_temporal_transvers = 74
numlabels = 0:size(destrieux_labels,1)-1;
destrieux_labels = [destrieux_labels, num2cell(numlabels')];

% mapping labels to colortable
for kk = 1:size(colortable_Destrieux.table,1) %76 are labels
    corr_numlabel = contains(destrieux_labels(:,1),colortable_Destrieux.struct_names{kk})==1;
    
    vert_label_Destrieux(vert_label_Destrieux==colortable_Destrieux.table(kk,5)) = destrieux_labels{corr_numlabel,2}; % unknown
end

% DKT ATLAS
% load DKT map
surface_labels_DKT = fullfile(cfg.freesurfer_directory,'label',...
    sprintf('%sh.aparc.DKTatlas.annot',cfg.hemisphere{1})); % corresponds to DKT40 atlas
[~, vert_label_DKT, colortable_DKT] = read_annotation(surface_labels_DKT);

% REORDER FREESURFER DKT LABELS TO DKT LABELS FROM Klein & Tourville Front. Neuroscience 2012
DKT_labels = colortable_DKT.struct_names;
idx1 = contains(DKT_labels,'corpuscallosum');
idx2 = contains(DKT_labels,'temporalpole');
idx3 = contains(DKT_labels,'bankssts');
idx4 = contains(DKT_labels,'frontalpole');
idx = idx1+idx2+idx3+idx4;
DKT_labels = DKT_labels(:).';
DKT_labels = [DKT_labels(~idx), DKT_labels(idx1), DKT_labels(idx2), DKT_labels(idx3), DKT_labels(idx4)];
DKT_labels = DKT_labels(:);

% unknown = 0
numlabels = 0:size(DKT_labels,1)-1;
DKT_labels = [DKT_labels, num2cell(numlabels')];

% mapping labels to colortable
for kk = 1:size(colortable_DKT.table,1) %35 are labels
    corr_numlabel = contains(DKT_labels(:,1),colortable_DKT.struct_names{kk})==1;
    
    vert_label_DKT(vert_label_DKT==colortable_DKT.table(kk,5)) = DKT_labels{corr_numlabel,2}; 
end

%load Wang map
Wang_ROI_Names = {...
    'V1v' 'V1d' 'V2v' 'V2d' 'V3v' 'V3d' 'hV4' 'VO1' 'VO2' 'PHC1' 'PHC2' ...
    'TO2' 'TO1' 'LO2' 'LO1' 'V3B' 'V3A' 'IPS0' 'IPS1' 'IPS2' 'IPS3' 'IPS4' ...
    'IPS5' 'SPL1' 'FEF'};
surface_labels_name = fullfile(cfg.freesurfer_directory,'surf',...
    sprintf('%sh.wang15_mplbl.mgz',cfg.hemisphere{1}));
if exist(surface_labels_name,'file')
    surface_labels = MRIread(surface_labels_name);
    vert_label_Wang = surface_labels.vol(:);
end
clear surface_labels_name

%load Benson map
Benson_Area_Names = {'V1','V2','V3','hV4','V01','V02','L01','L02','T01','T02','V3b','V3a'};
surface_labels_name = fullfile(cfg.freesurfer_directory,'surf',...
    sprintf('%sh.benson14_varea.mgz',cfg.hemisphere{1}));
if exist(surface_labels_name,'file')
    surface_labels_B = MRIread(surface_labels_name);
    vert_label_Benson = surface_labels_B.vol(:);
end
clear surface_labels_name

% load Benson Eccen
surface_labels_name = fullfile(cfg.freesurfer_directory,'surf',...
    sprintf('%sh.benson14_eccen.mgz',cfg.hemisphere{1}));
if exist(surface_labels_name,'file')
    surface_labels_B = MRIread(surface_labels_name);
    vert_eccen_label = surface_labels_B.vol(:);
    clear surface_labels surface_labels_name
end

%load Benson Angle
surface_labels_name = fullfile(cfg.freesurfer_directory,'surf',...
    sprintf('%sh.benson14_angle.mgz',cfg.hemisphere{1}));
if exist(surface_labels_name,'file')
    surface_labels_B = MRIread(surface_labels_name);
    vert_angle_label_B = surface_labels_B.vol(:);
    clear surface_labels surface_labels_name
end

% load Benson Sigma
surface_labels_name = fullfile(cfg.freesurfer_directory,'surf',...
    sprintf('%sh.benson14_sigma.mgz',cfg.hemisphere{1}));
if exist(surface_labels_name,'file')
    surface_labels_B = MRIread(surface_labels_name);
    vert_sigma_label = surface_labels_B.vol(:);
    clear surface_labels surface_labels_name
end

%%% DEFINE THE OUTPUT
DKTatlas_label = NaN(size(elecmatrix,1),1);
DKTatlas_label_text = cell(size(elecmatrix,1),1);
Destrieux_label = NaN(size(elecmatrix,1),1);
Destrieux_label_text = cell(size(elecmatrix,1),1);
Wang_label = NaN(size(elecmatrix,1),1);
Wang_label_text = cell(size(elecmatrix,1),1);
Benson_label = NaN(size(elecmatrix,1),1);
Benson_label_text = cell(size(elecmatrix,1),1);
Benson_eccen = NaN(size(elecmatrix,1),1);
Benson_polarangle = NaN(size(elecmatrix,1),1);
Benson_sigma = NaN(size(elecmatrix,1),1);

%%% LOOP THROUGH ELECTRODES AND ASSIGN LABELS

% define the range in which the code searches for the most represented
% area. So in this case it will look within 3 mm how which brain regions
% are found most, it takes the most represented one
electrode_to_vertex_dist = 3; % in mm

for elec = 1:size(elecmatrix,1) % loop across electrodes
    
    % dist from electrode to vertices
    b = sqrt(sum((g.vertices-repmat(elecmatrix(elec,:),size(g.vertices,1),1)).^2,2));
    [B,I] = sort(b,'ascend'); % sort distances from electrode to vertices with minimum on top
    
    %%%% DESTRIEUX:    
    if ~isnan(B(1))
        %take the label of the area with the shortest distance
        localized_electrodes = vert_label_Destrieux(I(1));
    else
        localized_electrodes = NaN;
    end

    % put the labels (vert_label) back in the matrix
    if ~isnan(localized_electrodes)
        Destrieux_label(elec,1) =  localized_electrodes;
        loc_destrieuxname = [destrieux_labels{:,2}]==localized_electrodes;
        Destrieux_label_text{elec,1} = destrieux_labels{loc_destrieuxname,1};
    else
        Destrieux_label(elec,1) =  NaN;
        Destrieux_label_text{elec,1} = 'n/a';
    end
    
    %%%% DKT:
    % take the mode of the labels within X mm
%     localized_electrodes = mode(vert_label_DKT(b<electrode_to_vertex_dist));
    if ~isnan(B(1))
        %take the label of the area with the shortest distance
        localized_electrodes = vert_label_DKT(I(1));
    else
        localized_electrodes = NaN;
    end
    
    % put the labels (vert_label) back in the matrix
    if localized_electrodes~=0 && ~isnan(localized_electrodes)
        
        DKTatlas_label(elec,1) =  localized_electrodes;
        loc_DKTname = [DKT_labels{:,2}]==localized_electrodes;        
        DKTatlas_label_text{elec,1} = DKT_labels{loc_DKTname,1};
        
    else
        DKTatlas_label(elec,1) =  NaN;
        DKTatlas_label_text{elec,1} = 'n/a';
    end
    
    %%%% WANG:
    % take the mode of the labels within X mm
    if exist('vert_label_Wang','var')
        area_of_electrode = mode(vert_label_Wang(b<electrode_to_vertex_dist));
        % put the labels (vert_label) back in the matrix
        Wang_label(elec,1) = area_of_electrode;
        if area_of_electrode > 0
            Wang_label_text{elec,1} = Wang_ROI_Names{area_of_electrode};
        else
            Wang_label_text{elec,1} = 'n/a';
        end
    else 
        Wang_label(elec,1) = NaN;
        Wang_label_text{elec,1} = 'n/a';
    end
    
    %%%% BENSON AREA:
    % take the mode of the labels within X mm
    if exist('vert_label_Benson','var')
        area_of_electrode = mode(vert_label_Benson(b<electrode_to_vertex_dist));
        % put the labels (vert_label) back in the matrix
        Benson_label(elec,1) =  area_of_electrode;
        if area_of_electrode>0
            Benson_label_text{elec,1} = Benson_Area_Names{area_of_electrode};
        else
            Benson_label_text{elec,1} = 'n/a';
        end
    else
        Benson_label(elec,1) = NaN;
        Benson_label_text{elec,1} = 'n/a';
    end
    
    % Only add eccentricity, angle and size if there is a Benson area
    if Benson_label(elec,1)>0
        %%%% BENSON ECCENTRICITY:
        % take the mean of the eccentricity within 3 mm
        eccen_of_electrode = mean(vert_eccen_label(b<electrode_to_vertex_dist));
        %put the labels (eccentricity) back in the matrix
        Benson_eccen(elec,1) = eccen_of_electrode;
        
        %%%% BENSON ANGLE:
        % take the mean of the polar angle within 3 mm
        angle_of_electrode = mean(vert_angle_label_B(b<electrode_to_vertex_dist));
        %put the labels (polar angle) back in the matrix
        Benson_polarangle(elec,1) = angle_of_electrode;
        
        %%%% BENSON SIGMA:
        % take the mean of the sigma within 3 mm
        sigma_of_electrode = mean(vert_sigma_label(b<electrode_to_vertex_dist));
        %put the labels (sigma) back in the matrix
        Benson_sigma(elec,1) = sigma_of_electrode;
    end
end


%%%% Integrating electrode positions with TSV-file of coordinates

% Make a new table with added variables
t = table(...
    Destrieux_label, Destrieux_label_text,...
    DKTatlas_label,DKTatlas_label_text,...
    Wang_label, Wang_label_text,...
    Benson_label,Benson_label_text, Benson_eccen,...
    Benson_polarangle, Benson_sigma);

% concatenate the table to what is already in loc_info
electrodes_tableWithlabels = horzcat(tb_elecs,t);

electrodes_tableWithlabels = bids_tsv_nan2na(electrodes_tableWithlabels);

% write table
% if ~exist(fullfile(cfg.ieeg_directory, ...
%         [cfg.sub_labels{:} '_' cfg.ses_label '_' output_file]),'file')
%     disp(['writing output ' output_file])
%     writetable(electrodes_tableWithlabels, ...
%         fullfile(cfg.ieeg_directory, ...
%         [cfg.sub_labels{:} '_' cfg.ses_label '_' output_file]),...
%         'Filetype','text','Delimiter','\t');
% else
%     disp(['ERROR: can not overwrite, output file already exists ' output_file])
% end
