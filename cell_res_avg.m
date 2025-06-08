SubjectNames = {'sub-0002', ...
    'sub-0004', 'sub-0006', 'sub-0007'};

output_dir = '/export02/data/vikramn/brainstorm3/resistances/';
thickness_areas = load('thickness_areas_omega_fem.mat');
%Use default anatomy for now
surface_file_paths = {
    '/export02/data/vikramn/TutorialOmega/anat/sub-0002/tess_cortex_central_low.mat',
    '/export02/data/vikramn/TutorialOmega/anat/sub-0004/tess_cortex_central_low.mat' ,
    '/export02/data/vikramn/TutorialOmega/anat/sub-0006/tess_cortex_central_low.mat' ,
    '/export02/data/vikramn/TutorialOmega/anat/sub-0007/tess_cortex_central_low.mat'
    };

num_v=15002;
dend_r=1.816*10^7;
gm_res=3.5; %ohm.m, Abboud et al. 2021
%Load the conversion from DKT regions to areal neuron densities
%Load surface file from brainstorm (oddball?)
rgn_neur_dens = readtable("./cell_props_filtered_redo.csv");
for s=1:numel(SubjectNames)
    subject= string(SubjectNames(s));
    surface_file = in_bst_data(char(surface_file_paths(s)));
    thickness = thickness_areas.subject_thicknesses(subject);
    support_areas = thickness_areas.subject_support_areas(subject);
    vertex_res_map_old=zeros(num_v, 1);
    vertex_resistance_map=zeros(num_v,1);
    for dkt_region=surface_file.Atlas(2).Scouts
        region=dkt_region.Label;
      
        % Find the row index where "brainstorm_name" is equal to "entorhinal L"
        rowIndex = find(strcmp(rgn_neur_dens.brainstorm_name, region));
        if(isempty(rowIndex))
            vertex_res_map_old(v) = gm_res*thickness(v)/support_areas(v);
            vertex_resistance_map(v) = gm_res*thickness(v)/support_areas(v);

            disp(region);
            continue
        end

        % Extract the corresponding value in the "neuron_area_density" column
        neuronDensityValue = rgn_neur_dens.neuron_area_density(rowIndex);
    
        vertices_in_rgn = dkt_region.Vertices;
        for i=1:numel(vertices_in_rgn)
            v=vertices_in_rgn(i);
            N=neuronDensityValue*(1000)^2*support_areas(v); %convert from neurons/mm^2 --> neurons/m^2
            R=dend_r/N;
            vertex_res_map_old(v) = gm_res*thickness(v)/support_areas(v);
            vertex_resistance_map(v)=R;
             
            
        end
    end
    save(sprintf('%s%s_region_res_cell_column.mat', output_dir, subject), 'vertex_res_map_old', ...
                'vertex_resistance_map');
end