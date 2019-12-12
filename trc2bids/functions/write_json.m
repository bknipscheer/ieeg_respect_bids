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
