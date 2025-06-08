%Prerequisite: 
%1. Run power_calc_pre_post_stim.m
%2. Export surface_file from Matlab
%3. Run preprocess_avg_p_good_trials.m
%load('sta_dip_avg_trial_p.mat');
%load('dev_dip_avg_trial_p.mat');

%NOTE - instead of loading a poststimulus average, you can just as easily do this
%with prestimulus/poststimulus powers for run 1/run 2 (with slight
%modifications) - replace load files above

%input_dir = '../resting_power/bands';
input_dir = '../resting_power/';
output_dir = './no_div_a/';
load('../thickness_areas_rest.mat');
load('s600_17_to_7.mat');

SubjectNames = {'sub-0002', 'sub-0003', ...
    'sub-0004', 'sub-0006', 'sub-0007'};
% freq_band_names = {'delta', 'theta', 'alpha', ...
%     'beta', 'gamma1', 'gamma2'};
atlas='s600';
load(sprintf('../%s_region_names.mat', atlas));
num_seconds = 100;

region_dip_pvt_dict = containers.Map();
region_dip_ivt_dict = containers.Map();

band=''; %default - when working with non-band-pass filtered data
%also remove the for loop below
dk_index = 2;
%schaeffer_index = 4; %schaeffer 200, 17 networks
schaeffer_index=6;
% for b=1:numel(freq_band_names) %can always remove later
%     band=freq_band_names{b};
    %Load surface file for each subject
    s17_names_to_ignore = {'Background+FreeSurfer_Defined_Medial_Wall L', ...
     'Background+FreeSurfer_Defined_Medial_Wall R'   };
    sub_rgn_p_dicts = containers.Map();
    sub_rgn_i_dicts = containers.Map();
    for s=1:numel(SubjectNames)
        subject=SubjectNames{s}; 
        sub_anat_filename = sprintf(['/export02/data/vikramn/TutorialOmega/anat/' ...
            '%s/tess_cortex_pial_low.mat'], subject);
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
       
        dip_avg_sec_p = findAndAverageFiles(input_dir, subject, band, 'p'); %avg. power of each second for 100 s
        dip_avg_sec_i = findAndAverageFiles(input_dir, subject, band, 'i');
    
        region_dip_pvt_dict = containers.Map();
        region_dip_ivt_dict = containers.Map();
        for region_list=region_names
            region_name = region_list{1};
            vertices_in_rgn = region_v_dict(region_name);
    
            dipole_p_v_t = cell2mat(dip_avg_sec_p);
            dipole_i_v_t = cell2mat(dip_avg_sec_i);
    
            region_dip_pvt_dict(region_name) = dipole_p_v_t(vertices_in_rgn,:); %NOTE - this is average power per TRIAL
            region_dip_ivt_dict(region_name) = dipole_i_v_t(vertices_in_rgn,:);
                
        end
    
        region_pvt_dict = containers.Map();
        region_ivt_dict = containers.Map();
    
        fprintf("Subject: %s\n", subject);
        max_power = 0;
        max_power_region = '';
        max_curr = 0;
        max_curr_region = '';
        for region_list=region_names
            region_name=region_list{1};
            %a = region_area_dict(region_name);
            %p_per_m2 = sum(region_dip_pvt_dict(region_name),1)/a;
            total_p = sum(region_dip_pvt_dict(region_name),1);
            %region_pvt_dict(region_name) = p_per_m2;
            region_pvt_dict(region_name) = total_p;
            %i_per_m2 =  sum(region_dip_ivt_dict(region_name),1)/a;
            total_i = sum(region_dip_ivt_dict(region_name),1);
            %region_ivt_dict(region_name) = i_per_m2;
            region_ivt_dict(region_name) = total_i;
        
%             if max(p_per_m2) > max_power
%                 max_power = max(p_per_m2);
%                 max_power_region = region_name;
%             end
%             if max(i_per_m2) > max_curr
%                 max_curr = max(i_per_m2);
%                 max_curr_region = region_name;
%             end
        
        end
%     
%         fprintf("Subject: %s\n", subject);
%         fprintf('Maximum Power/m^2: %s, Region: %s\n', num2str(max_power), max_power_region);
%         fprintf('Maximum Curr/m^2: %s, Region: %s\n', num2str(max_curr), max_curr_region);
%         fprintf("*************************************************************\n");
    
        save(sprintf('%s%s_%s_region_100s_dip_p_i.mat', output_dir, atlas, subject), 'region_dip_pvt_dict', ...
            'region_dip_ivt_dict');
    
        save(sprintf('%s%s_%s_region_100s_p_i.mat', output_dir, atlas, subject), 'region_pvt_dict', ...
            'region_ivt_dict');
%     
%         save(sprintf('%s_%s_%s_region_100s_dip_p_i.mat', atlas, subject, band), 'region_dip_pvt_dict', ...
%             'region_dip_ivt_dict');
%     
%         save(sprintf('%s_%s_%s_region_100s_p_i.mat', atlas, subject, band), 'region_pvt_dict', ...
%             'region_ivt_dict');
    
        produceOutputCSV(region_ivt_dict, sprintf('%s%s_%s_i_100s_rest.csv', output_dir, atlas, subject))
        produceOutputCSV(region_pvt_dict, sprintf('%s%s_%s_p_100s_rest.csv', output_dir, atlas, subject))
% 
%           produceOutputCSV(region_ivt_dict, sprintf('%s_%s_%s_i_100s_rest.csv', atlas, subject, band))
%           produceOutputCSV(region_pvt_dict, sprintf('%s_%s_%s_p_100s_rest.csv', atlas, subject, band))
    
        sub_rgn_p_dicts(subject) = region_pvt_dict;
        sub_rgn_i_dicts(subject) = region_ivt_dict;
    
    end
    
    avg_rgn_p_dict = containers.Map();
    avg_rgn_i_dict = containers.Map();
    for region_list=region_names
        region_name=region_list{1};
        rgn_total_power_vec = zeros(1, num_seconds);
        rgn_total_current_vec = zeros(1, num_seconds);
    
        for s=1:numel(SubjectNames)
            sub_p_dict = sub_rgn_p_dicts(subject);
            rgn_power_100s = sub_p_dict(region_name);
            rgn_total_power_vec = rgn_total_power_vec + rgn_power_100s;
    
            sub_i_dict = sub_rgn_i_dicts(subject);
            rgn_curr_100s = sub_i_dict(region_name);
            rgn_total_current_vec = rgn_total_current_vec + rgn_curr_100s;
        end
    
        sub_avg_p_ts = rgn_total_power_vec/numel(SubjectNames);
        sub_avg_i_ts = rgn_total_current_vec/numel(SubjectNames);
    
        avg_rgn_p_dict(region_name) = sub_avg_p_ts;
        avg_rgn_i_dict(region_name) = sub_avg_i_ts;
    end
    
    %save(sprintf('%s_avg_rgn_p_i_dicts.mat',band),'avg_rgn_p_dict','avg_rgn_i_dict');
    save(sprintf('%s%s_avg_rgn_p_i_dicts.mat',output_dir, atlas),'avg_rgn_p_dict','avg_rgn_i_dict');
   % produceOutputCSV(avg_rgn_p_dict, sprintf('%s_%s_p_100s_rest_subavg.csv', atlas, band));
   % produceOutputCSV(avg_rgn_i_dict, sprintf('%s_%s_i_100s_rest_subavg.csv', atlas, band));
     produceOutputCSV(avg_rgn_p_dict, sprintf('%s%s_p_100s_rest_subavg.csv', output_dir, atlas));
     produceOutputCSV(avg_rgn_i_dict, sprintf('%s%s_i_100s_rest_subavg.csv', output_dir, atlas));
%end
%Repeat the above code and further decompose by frequency (see above)

%HELPER FUNCTIONS
