function create_json_mri(fmri_json_name)

mri_json    = [];

mri_json.Instruction                    = 'Lie still';
mri_json.EyesOpenClosed                 = 'Not defined';

%% write ieeg.json


filename = fmri_json_name;

if ~isempty(filename)
    if isfile(filename)
        existing = read_json(filename);
    else
        existing = [];
    end
    write_json(filename, mergeconfig(existing, mri_json))
end


end