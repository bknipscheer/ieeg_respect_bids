
function create_elecDesc(proj_dir,cfg)

destrieux_labels = cfg.destrieux_labels;
DKT_labels = cfg.DKT_labels;

elecdesc_json.name                  = 'Name of the electrode';
elecdesc_json.x_y_z                 = 'X, y and z position of the electrode on the brain of the subject';
elecdesc_json.size                  = 'Surface size in mm2 of the electrode';
elecdesc_json.group                 = 'Group to which electrode belongs, this can be grid, strip, depth or other';
elecdesc_json.material              = 'Material of the electrode. This is platinum in most situations';
elecdesc_json.manufacturer          = 'Manufacturer of the electrode';
elecdesc_json.silicon               = 'When an electrode is overlapping with another electrode, the electrode on top is recording silicon. These electrodes should be excluded from analysis!';
elecdesc_json.soz                   = 'An electrode located on the seizure onset zone';
elecdesc_json.ra                    = 'An electrode located on tissue that is resected after recordings during surgery (Resected Area)';
elecdesc_json.edge                  = 'An electrode located on the edge of the resected area. It is unclear whether it is resected or not.';

strings_destrieux = cell(1,size(destrieux_labels,1));
for i=1:size(destrieux_labels,1)
    if i< size(destrieux_labels,1)
        strings_destrieux{i} = sprintf('%d = %s, ',destrieux_labels{i,2},destrieux_labels{i,1})    ;
    elseif i == size(destrieux_labels,1)
        strings_destrieux{i} = sprintf('%d = %s',destrieux_labels{i,2},destrieux_labels{i,1})    ;
    end
end

strings_DKT = cell(1,size(DKT_labels,1));
for i=1:size(DKT_labels,1)
    if i< size(DKT_labels,1)
        strings_DKT{i} = sprintf('%d = %s, ',DKT_labels{i,2},DKT_labels{i,1})    ;
    elseif i == size(destrieux_labels,1)
        strings_DKT{i} = sprintf('%d = %s',DKT_labels{i,2},DKT_labels{i,1})    ;
    end
end

elecdesc_json.Destrieux_label       = 'Electrode location in a region according to Destrieux et al. NeuroImage2011';
elecdesc_json.Destrieux_label_text  = ['Electrode location in a region according to Destrieux et al. NeuroImage2011: ',[strings_destrieux{:}]];
elecdesc_json.DKT_label             = 'Electrode location in a region according to Klein & Tourville Front. Neuroscience 2012';
elecdesc_json.DKT_label_text        = ['Electrode location in a region according to Klein & Tourville Front. Neuroscience 2012: ',[strings_DKT{:}]];
elecdesc_json.Wang_label            = 'Electrode location in a visual cortex region according to Wang';
elecdesc_json.Wang_label_text       = 'Electrode location in a visual cortex region according to Wang';
elecdesc_json.Benson_label          = 'Electrode location in a visual cortex region according to Benson';
elecdesc_json.Benson_label_text     = 'Electrode location in a visual cortex region according to Benson';
elecdesc_json.Benson_eccen          = 'Electrode location in a visual cortex region according to Benson';
elecdesc_json.Benson_polarangle     = 'Electrode location in a visual cortex region according to Benson';
elecdesc_json.Benson_sigma          = 'Electrode location in a visual cortex region according to Benson';

if ~isempty(elecdesc_json)
    
    filename = fullfile(proj_dir,'electrode_description.json');
    if isfile(filename)
        existing = read_json(filename);
    else
        existing = [];
    end
    write_json(filename, mergeconfig(existing, elecdesc_json))
end
end

function json = read_json(filename)
ft_info('reading %s\n', filename);
if ft_hastoolbox('jsonlab', 3)
    json = loadjson(filename);
else
    fid = fopen(filename, 'r');
    str = fread(fid, [1 inf], 'char=>char');
    fclose(fid);
    json = jsondecode(str);
end
end

function write_json(filename, json)
json = remove_empty(json);
ft_info('writing %s\n', filename);
if ft_hastoolbox('jsonlab', 3)
    savejson('', json, filename);
else
    str = jsonencode(json);
    fid = fopen(filename, 'w');
    fwrite(fid, str);
    fclose(fid);
end
end

function s = remove_empty(s)
fn = fieldnames(s);
fn = fn(structfun(@isempty, s));
s = removefields(s, fn);
end