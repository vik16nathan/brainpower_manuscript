%Go from Schaeffer 200 17 atlas --> 7 atlas
schaeffer_7_index=11;
schaeffer_17_index=10;
s200_17_v_dict = containers.Map();

s200_7_v_dict = containers.Map();

for region=surface_file_aud.Atlas(schaeffer_17_index).Scouts %change atlas in brainstorm
    region_name=region.Label;
    vertices_in_rgn = region.Vertices;
    s200_17_v_dict(region_name) = vertices_in_rgn;
end


for region=surface_file_aud.Atlas(schaeffer_7_index).Scouts %change atlas in brainstorm
    region_name=region.Label;
    vertices_in_rgn = region.Vertices;
    s200_7_v_dict(region_name) = vertices_in_rgn;
end

