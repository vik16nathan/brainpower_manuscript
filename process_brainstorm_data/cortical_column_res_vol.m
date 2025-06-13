%Map of L/A (resistance) and L*A (volume) across the brains of the resting
%individuals

%load('../thickness_areas_rest.mat');
load('../../thickness_areas_omega_fem.mat');
avg_gm_res = 3.5; %ohm*meters
%SubjectNames = {'sub-0002', 'sub-0003', ...
%    'sub-0004', 'sub-0006', 'sub-0007'};
SubjectNames = {'sub-0002', ...
    'sub-0004', 'sub-0006', 'sub-0007'};
% freq_band_names = {'delta', 'theta', 'alpha', ...
%     'beta', 'gamma1', 'gamma2'};
%atlas='s600';
atlas='s600';
load(sprintf('%s_17_to_7.mat', atlas));
load(sprintf('../%s_region_names.mat', atlas));
band=''; %default - when working with non-band-pass filtered data
%also remove the for loop below
dk_index = 2;
%schaeffer_index = 4; %schaeffer 200, 17 networks
schaeffer_index=6; %schaefer 600, 17 networks

s17_names_to_ignore = {'Background+FreeSurfer_Defined_Medial_Wall L', ...
 'Background+FreeSurfer_Defined_Medial_Wall R'   };

%These are AVERAGES for each vertex in the parcel.
sub_rgn_vol_dicts = containers.Map(); %thickness x support area
sub_rgn_res_dicts = containers.Map(); %L/A, assuming uniform resistivity across all regions
sub_rgn_thick_dicts = containers.Map(); %mean cortical thickness
sub_rgn_sa_dicts = containers.Map();    %mean dipole support area (note: this is not the same as the parcel area)
for s=1:numel(SubjectNames)
    subject=SubjectNames{s}; 
    sub_anat_filename = sprintf(['/export02/data/vikramn/TutorialOmega/anat/' ...
        '%s/tess_cortex_central_low.mat'], subject);
    support_areas = subject_support_areas(subject);
    thickness = subject_thicknesses(subject);
    surface_file = load(sub_anat_filename);
    region_vol_dict = containers.Map();
    region_res_dict = containers.Map();
    region_thick_dict = containers.Map();
    region_sa_dict = containers.Map();
    for region=surface_file.Atlas(schaeffer_index).Scouts %change atlas in brainstorm
        %IGNORE MEDIAL WALL

        region_name_s17=region.Label;
        if ismember(region_name_s17, s17_names_to_ignore)
            continue;
        end

        %region_name = s600_17_to_7(region_name_s17); %7 network AUTOMATE
        %LATER
        region_name = s600_17_to_7(region_name_s17);
        vertices_in_rgn = region.Vertices;
        num_v = length(vertices_in_rgn);
        region_res_dict(region_name) = avg_gm_res*sum(thickness(vertices_in_rgn) ./ ...
            support_areas(vertices_in_rgn))/num_v;
        %normalize by parcel area before correlating with power/current?
        region_vol_dict(region_name) = sum(support_areas(vertices_in_rgn) .* ...
            thickness(vertices_in_rgn))/num_v;
        region_thick_dict(region_name) = sum(thickness(vertices_in_rgn))/num_v;
        region_sa_dict(region_name) = sum(support_areas(vertices_in_rgn))/num_v;
        
    end
    sub_rgn_vol_dicts(subject) = region_vol_dict;
    sub_rgn_res_dicts(subject) = region_res_dict;
    sub_rgn_thick_dicts(subject) = region_thick_dict;
    sub_rgn_sa_dicts(subject) = region_sa_dict;
end

%Average across subjects

avg_rgn_vol_dict = containers.Map();
avg_rgn_res_dict = containers.Map();
avg_rgn_thick_dict = containers.Map();
avg_rgn_sa_dict = containers.Map();

for region_list=region_names
    region_name=region_list{1};
    rgn_total_vol = 0;
    rgn_total_res = 0;
    rgn_total_thick=0;
    rgn_total_sa=0;

    for s=1:numel(SubjectNames)
        sub_vol_dict = sub_rgn_vol_dicts(subject);
        rgn_vol = sub_vol_dict(region_name);
        rgn_total_vol = rgn_total_vol + rgn_vol;

        sub_res_dict = sub_rgn_res_dicts(subject);
        rgn_res = sub_res_dict(region_name);
        rgn_total_res = rgn_total_res + rgn_res;

        sub_thick_dict = sub_rgn_thick_dicts(subject);
        rgn_thick = sub_thick_dict(region_name);
        rgn_total_thick = rgn_total_thick + rgn_thick;

        sub_sa_dict = sub_rgn_sa_dicts(subject);
        rgn_sa = sub_sa_dict(region_name);
        rgn_total_sa = rgn_total_sa + rgn_sa;
    end

    sub_avg_vol = rgn_total_vol/numel(SubjectNames);
    sub_avg_res = rgn_total_res/numel(SubjectNames);
    sub_avg_thick = rgn_total_thick/numel(SubjectNames);
    sub_avg_sa = rgn_total_sa/numel(SubjectNames);

    avg_rgn_vol_dict(region_name) = sub_avg_vol;
    avg_rgn_res_dict(region_name) = sub_avg_res;
    avg_rgn_thick_dict(region_name) = sub_avg_thick;
    avg_rgn_sa_dict(region_name) = sub_avg_sa;
end

save(sprintf('%s_avg_rgn_vol_res_dicts.mat',atlas),'avg_rgn_vol_dict','avg_rgn_res_dict', ...
    'avg_rgn_thick_dict', 'avg_rgn_sa_dict');
produceOutputCSV(avg_rgn_vol_dict, sprintf('%s_subavg_column_vol.csv', atlas));
produceOutputCSV(avg_rgn_res_dict, sprintf('%s_subavg_column_res.csv', atlas));
produceOutputCSV(avg_rgn_thick_dict, sprintf('%s_subavg_column_thick.csv', atlas));
produceOutputCSV(avg_rgn_sa_dict, sprintf('%s_subavg_column_sa.csv', atlas));