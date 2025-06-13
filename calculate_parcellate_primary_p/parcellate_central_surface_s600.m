load('s600_17_to_7.mat');
atlas='s600';
output_dir='./cent_surf_fem_fwd/'
load(sprintf('../%s_region_names.mat', atlas));
schaeffer_index=6; %which atlas in brainstorm
SubjectNames={'sub-0002','sub-0004', 'sub-0006', 'sub-0007'};
load('../../thickness_areas_omega_fem.mat');
s17_names_to_ignore = {'Background+FreeSurfer_Defined_Medial_Wall L', ...
     'Background+FreeSurfer_Defined_Medial_Wall R'   };
    sub_rgn_p_dicts = containers.Map();
    sub_rgn_i_dicts = containers.Map();
    for s=1:numel(SubjectNames)
        subject=SubjectNames{s}; 
        sub_anat_filename = sprintf(['/export02/data/vikramn/TutorialOmega/anat/' ...
            '%s/tess_cortex_central_low.mat'], subject); %NOT pial
        support_areas = subject_support_areas(subject);
        surface_file = load(sub_anat_filename);
        region_v_dict = containers.Map();
        region_area_dict = containers.Map();
        for region=surface_file.Atlas(schaeffer_index).Scouts %change atlas in brainstorm
            %IGNORE MEDIAL WALL
    
            region_name_s17=region.Label;
            if ismember(region_name_s17, s17_names_to_ignore)
                continue;
            end
    
            region_name = s600_17_to_7(region_name_s17); %7 network
            vertices_in_rgn = region.Vertices;
            region_v_dict(region_name) = vertices_in_rgn;
            region_area_dict(region_name) = sum(support_areas(vertices_in_rgn));
        end
      
        produceOutputCSV(region_area_dict,sprintf('%s_rgn_area_dict.csv',subject));

        sub_fem_pdip_sec=load(sprintf('../../../duneuro_malte/%s_fem_p_0s.csv',subject));
        %sub_fem_pdip_pri=load(sprintf('%s_dip_p_0s.mat',subject));

        %load L2 norm of 0s secondary current in all three dimensions
        sub_fem_idip_pri=load(sprintf('%s_dip_i_0s.mat',subject));

        %dip_avg_sec_p = findAndAverageFiles(input_dir, subject, band, 'p'); %avg. power of each second for 100 s
        %dip_avg_sec_i = findAndAverageFiles(input_dir, subject, band, 'i');
    
        region_dip_pvt_dict = containers.Map();
        region_dip_ivt_dict = containers.Map();
        for region_list=region_names
            region_name = region_list{1};
            vertices_in_rgn = region_v_dict(region_name);
    
            dipole_p_v_t = sub_fem_pdip_sec; %REDO WITH PRIMARY
            %dipole_p_v_t = sub_fem_pdip_pri.dip_p;
            dipole_i_v_t = sub_fem_idip_pri.dip_i;
    
            region_dip_pvt_dict(region_name) = dipole_p_v_t(vertices_in_rgn,:); %NOTE - this is average power per TRIAL
            region_dip_ivt_dict(region_name) = dipole_i_v_t(vertices_in_rgn,:);
                
        end
    
        region_pvt_dict = containers.Map();
        region_ivt_dict = containers.Map();
    
        fprintf("Subject: %s\n", subject);
        for region_list=region_names
            region_name=region_list{1};
            total_p = sum(region_dip_pvt_dict(region_name),1);
            region_pvt_dict(region_name) = total_p;
            total_i = sum(region_dip_ivt_dict(region_name),1);
            region_ivt_dict(region_name) = total_i;
       
        end
    
        save(sprintf('%s%s_%s_region_0s_sec_dip_p_i.mat', output_dir, atlas, subject), 'region_dip_pvt_dict', ...
            'region_dip_ivt_dict');
    
        save(sprintf('%s%s_%s_region_0s_sec_p_i.mat', output_dir, atlas, subject), 'region_pvt_dict', ...
            'region_ivt_dict');
    
        produceOutputCSV(region_ivt_dict, sprintf('%s%s_%s_i_0s_2min_rest.csv', output_dir, atlas, subject))
        produceOutputCSV(region_pvt_dict, sprintf('%s%s_%s_p_sec_0s_2min_rest.csv', output_dir, atlas, subject))
    
        sub_rgn_p_dicts(subject) = region_pvt_dict;
        sub_rgn_i_dicts(subject) = region_ivt_dict;
    
    end
    
    %Subject average - to compre with PET averaged across individuals
    %note that sample sizes are different
    avg_rgn_p_dict = containers.Map();
    avg_rgn_i_dict = containers.Map();
    for region_list=region_names
        region_name=region_list{1};
        rgn_total_power_vec = zeros(1, 1);
        rgn_total_current_vec = zeros(1, 1);
    
        for s=1:numel(SubjectNames)
            sub_p_dict = sub_rgn_p_dicts(subject);
            rgn_power = sub_p_dict(region_name);
            rgn_total_power_vec = rgn_total_power_vec + rgn_power;
    
            sub_i_dict = sub_rgn_i_dicts(subject);
            rgn_curr = sub_i_dict(region_name);
            rgn_total_current_vec = rgn_total_current_vec + rgn_curr;
        end
    
        sub_avg_p_ts = rgn_total_power_vec/numel(SubjectNames);
        sub_avg_i_ts = rgn_total_current_vec/numel(SubjectNames);
    
        avg_rgn_p_dict(region_name) = sub_avg_p_ts;
        avg_rgn_i_dict(region_name) = sub_avg_i_ts;
    end
   
    %save(sprintf('%s%s_avg_rgn_p_i_dicts.mat',output_dir, atlas),'avg_rgn_p_dict','avg_rgn_i_dict');
    produceOutputCSV(avg_rgn_p_dict, sprintf('%s%s_p_sec_0s_2min_rest_subavg.csv', output_dir, atlas));
    produceOutputCSV(avg_rgn_i_dict, sprintf('%s%s_i_0s_2min_rest_subavg.csv', output_dir, atlas));
