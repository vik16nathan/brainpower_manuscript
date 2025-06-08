%Goal: modify power_calc_vec to store the per-dipole power consumption at
%each time point, and track how it changes in each region of the brain
%during the oddball task (tutorial data)

import containers.Map;

%Initialize global variables
nrow = 45006;
ncol = 361;
total_t = 361;
num_v = nrow/3;
gm_res = 3.5; %grey matter resistivity in ohms*meters (from Bouda et al., 2021)
%estimate from Lennie et al., 2003 (whole-brain)

% 1. Get the x, y, and z currents at each source. → current dipole
%15k vertices
scaling_factor=1;

%Repeat for both runs
%run_1_filenames = {'link|Subject01/S01_AEF_20131218_01_600Hz_notch/results_MN_MEG_KERNEL_230531_1033.mat|Subject01/S01_AEF_20131218_01_600Hz_notch/data_deviant_average_230528_1851.mat',...
    %'link|Subject01/S01_AEF_20131218_01_600Hz_notch/results_MN_MEG_KERNEL_230531_1033.mat|Subject01/S01_AEF_20131218_01_600Hz_notch/data_standard_average_230528_1851.mat'};
%run_2_filenames = {'link|Subject01/S01_AEF_20131218_02_600Hz_notch/results_MN_MEG_KERNEL_230626_1246.mat|Subject01/S01_AEF_20131218_02_600Hz_notch/data_deviant_average_230528_1851.mat',...
    %'link|Subject01/S01_AEF_20131218_02_600Hz_notch/results_MN_MEG_KERNEL_230626_1247.mat|Subject01/S01_AEF_20131218_02_600Hz_notch/data_standard_average_230528_1851.mat'};

%Redo with the same number of standard/deviant trials before averaging - 35 trials
run_1_filenames = {'link|Subject01/S01_AEF_20131218_01_600Hz_notch/results_MN_MEG_KERNEL_230704_1525.mat|Subject01/S01_AEF_20131218_01_600Hz_notch/data_deviant_average_230704_1519.mat',...
                   'link|Subject01/S01_AEF_20131218_01_600Hz_notch/results_MN_MEG_KERNEL_230704_1526.mat|Subject01/S01_AEF_20131218_01_600Hz_notch/data_standard_average_230704_1519.mat'};

run_2_filenames = {'link|Subject01/S01_AEF_20131218_02_600Hz_notch/results_MN_MEG_KERNEL_230704_1540.mat|Subject01/S01_AEF_20131218_02_600Hz_notch/data_deviant_average_230704_1519.mat',...
    'link|Subject01/S01_AEF_20131218_02_600Hz_notch/results_MN_MEG_KERNEL_230704_1541.mat|Subject01/S01_AEF_20131218_02_600Hz_notch/data_standard_average_230704_1519.mat'};

run_filenames = {run_1_filenames, run_2_filenames};

sta_p_v_t_full = {};
sta_dip_v_t_full = {};
dev_p_v_t_full = {};
dev_dip_v_t_full = {};

%for run=[1,2]
run=1; %Focus only on run 1 for now (can always repeat with run = 2)
curr_dev = in_bst_results(run_filenames{run}{1}, 1); 
curr_sta = in_bst_results(run_filenames{run}{2}, 1);

curr_sta_amp = scaling_factor*curr_sta.ImageGridAmp;
curr_dev_amp = scaling_factor*curr_dev.ImageGridAmp;

curr_sta_xyz = convertTo3DArray(curr_sta_amp, nrow, ncol);
curr_dev_xyz = convertTo3DArray(curr_dev_amp, nrow, ncol);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%2. Get the cortical thickness at each source. → wire length
%ISSUE - ARBITRARY UNITS FOR THICKNESS

%projected to 15k sources
thickness_file = in_bst_results('C:\Users\vik16\OneDrive\Documents\brainstorm_db\TutorialIntroduction\data\Group_analysis\CAT12\results_surface_thickness_230601_1116_Subject01.mat',1);

%full 200k sources
%thickness_file = in_bst_results('C:\Users\vik16\OneDrive\Documents\brainstorm_db\TutorialIntroduction\data\Subject01\CAT12\results_surface_thickness_230601_1116.mat',1);
thickness = thickness_file.ImageGridAmp(:,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Need to export pial_15002V as a variable called surface_file beforehand. 
% How do I automate this 

%3. Get all triangles surrounding each source in the tessellated edge map.
triangles = surface_file.Faces;
vertices_xyz = surface_file.Vertices;

%Create a dictionary of triangle areas attached to each vertex
v_area_dict = create_tri_area_dict(triangles, vertices_xyz);

%For each vertex (source in our current source model), we assume a single "wire" in the cortex with 
%length equal to the cortical thickness and support area equal to 1/3*the
%area of the attached triangles, where current is uniformly distributed. 

%We use this 'wire' to calculate power at each dipole.
%Look at unconstrained currents only (x, y, z rather than normal to cortex)
%Calculate power at each vertex using Ohm's law and efficient matrix
%operations
thickness_in_m = thickness / 1000;
support_areas = oneThirdSumAreas(num_v, v_area_dict);

sta_power_calc = calculatePower(curr_sta_xyz, thickness_in_m, support_areas, gm_res, num_v, total_t);
sta_power_v_t = sta_power_calc{1};
sta_dipole_power_v_t = sta_power_calc{2};

dev_power_calc = calculatePower(curr_dev_xyz, thickness_in_m, support_areas, gm_res, num_v, total_t);
dev_power_v_t = dev_power_calc{1};
dev_dipole_power_v_t = dev_power_calc{2};

%sta_p_v_t_full = [sta_p_v_t_full, sta_power_v_t];
%sta_dip_v_t_full = [sta_dip_v_t_full, sta_dipole_power_v_t];
%dev_p_v_t_full = [dev_p_v_t_full, dev_power_v_t];
%dev_dip_v_t_full = [dev_dip_v_t_full, dev_dipole_power_v_t];   


dk_region_names={};
region_v_dict = containers.Map();
region_area_dict = containers.Map();
region_sta_dip_pvt_dict = containers.Map();
region_dev_dip_pvt_dict = containers.Map();

for dk_region=surface_file.Atlas(2).Scouts
    region=dk_region.Label;
    dk_region_names = [dk_region_names, region];
    vertices_in_rgn = dk_region.Vertices;
    region_v_dict(region) = vertices_in_rgn;
    region_area_dict(region) = sum(support_areas(vertices_in_rgn));
    region_sta_dip_pvt_dict(region) = sta_dipole_power_v_t(vertices_in_rgn,:);
    region_dev_dip_pvt_dict(region) = dev_dipole_power_v_t(vertices_in_rgn,:);
end

%Calculate the power (normalized by area) for each region
region_sta_pvt_dict = containers.Map();
region_dev_pvt_dict = containers.Map();
figure;
max_sta_power = 0;
max_sta_power_region = '';
max_dev_power = 0;
max_dev_power_region = '';

for region_list=dk_region_names
    region=region_list{1};
    a = region_area_dict(region);
    sta_pwr = sum(region_sta_dip_pvt_dict(region),1)/a;
    region_sta_pvt_dict(region) = sum(region_sta_dip_pvt_dict(region),1)/a;
    dev_pwr =  sum(region_dev_dip_pvt_dict(region),1)/a;
    region_dev_pvt_dict(region) = dev_pwr;

    if max(sta_pwr) > max_sta_power
        max_sta_power = max(sta_pwr);
        max_sta_power_region = region;
    end
    if max(dev_pwr) > max_dev_power
        max_dev_power = max(dev_pwr);
        max_dev_power_region = region;
    end
end

fprintf('Maximum Standard Power: %s, Region: %s\n', num2str(max_sta_power), max_sta_power_region);
fprintf('Maximum Deviant Power: %s, Region: %s\n', num2str(max_dev_power), max_dev_power_region);

%Find DK atlas regions that consume largest power (can easily rerun for
%smallest power)
%label the regions with the five highest average overall power consumptions
rgn_avg_sta_pwr = containers.Map();
rgn_avg_dev_pwr = containers.Map();
for region_list=dk_region_names
    region=region_list{1};
    sta_pvt =  region_sta_pvt_dict(region);
    rgn_avg_sta_pwr(region) = mean(sta_pvt);
    dev_pvt = region_dev_pvt_dict(region);
    rgn_avg_dev_pwr(region) = mean(dev_pvt);
end

%Find the regions with the largest average standard/deviant power values
rgn_avg_sta_pwr_desc = sortMapDescend(rgn_avg_sta_pwr);
rgn_avg_dev_pwr_desc = sortMapDescend(rgn_avg_dev_pwr);

%Visualization 
subplot(2,1,1);
%highlight top five regions for power consumption 
top_5_sta_regions = {};
for i=1:size(rgn_avg_sta_pwr_desc,1)
    region = rgn_avg_sta_pwr_desc{i,1};
    if i <= 5
        plot(1/600*[1:total_t], region_sta_pvt_dict(region), 'DisplayName', region, 'LineWidth', 2);
        top_5_sta_regions = [top_5_sta_regions, region];
    
    else
        plot(1/600*[1:total_t], region_sta_pvt_dict(region));
        ylim([2*10^-11, 2.5*10^-9]);
    end
    hold on;
    i=i+1;
end
legend([gca().Children(1) gca().Children(2) gca().Children(3) gca().Children(4) gca().Children(5)], top_5_sta_regions); 
title('Average Regional Power Consumption Across Standard Tones')
xlabel('Time (s)')
ylabel('Normalized power (W/m^2)')
 
subplot(2,1,2);
%highlight top five regions for power consumption 
top_5_dev_regions = {};
for i=1:size(rgn_avg_dev_pwr_desc,1)
    region = rgn_avg_dev_pwr_desc{i,1};
    if i <= 5
        plot(1/600*[1:total_t], region_dev_pvt_dict(region), 'DisplayName', region, 'LineWidth', 2);
        top_5_dev_regions = [top_5_dev_regions, region];
    
    else
        plot(1/600*[1:total_t], region_dev_pvt_dict(region));
        ylim([2*10^-11, 2.5*10^-9]);
    end
    hold on;
    i=i+1;
end
legend([gca().Children(1) gca().Children(2) gca().Children(3) gca().Children(4) gca().Children(5)], top_5_dev_regions); 
title('Average Regional Power Consumption Across Deviant Tones')
xlabel('Time (s)')
ylabel('Normalized power (W/m^2)')

region_sta_pvt_split = splitMapLR(region_sta_pvt_dict);
region_dev_pvt_split = splitMapLR(region_dev_pvt_dict);
produceOutputCSV(region_sta_pvt_split{1}, 'dk_sta_pvt_35trial_l.csv')
produceOutputCSV(region_sta_pvt_split{2}, 'dk_sta_pvt_35trial_r.csv')
produceOutputCSV(region_dev_pvt_split{1}, 'dk_dev_pvt_35trial_l.csv')
produceOutputCSV(region_dev_pvt_split{2}, 'dk_dev_pvt_35trial_r.csv')

%HELPER FUNCTIONS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function sortedArray=sortMapDescend(myMap)
    % Convert the containers.Map object to a cell array of key-value pairs
    keys = myMap.keys();
    values = myMap.values();
    cellArray = cell(68, 2); %there are 68 DK atlas parcellations
    cellArray(:,1) = keys;
    cellArray(:,2) = values;
    % Sort the cell array based on the values in descending order
    sortedArray = sortrows(cellArray, -2);
    
    % Convert the sorted cell array back to a containers.Map object
    disp(sortedArray);
end

function out=produceOutputCSV(myMap, filename)
    keysCell = keys(myMap);
    valuesCell = values(myMap);
    dataCell = [keysCell', valuesCell'];
    
    % Convert cell array to table
    dataTable = cell2table(dataCell, 'VariableNames', {'Keys', 'Values'});
    
    % Write table to CSV file
    writetable(dataTable, filename, 'Delimiter', ',');
    
    disp('CSV file exported successfully.');
end

function output=splitMapLR(myMap)
    leftMap = containers.Map();
    rightMap = containers.Map();
    keysCell = keys(myMap);
    valuesCell = values(myMap);
    for i = 1:numel(keysCell)
        key = keysCell{i};
        value = valuesCell{i};
        
        % Check if the key ends with 'L'
        if endsWith(key, ' L')
            newKey = strrep(key, ' L', '');
            leftMap(newKey) = value;
        end
        
        % Check if the key ends with 'R'
        if endsWith(key, ' R')
            newKey = strrep(key, ' R', '');
            rightMap(newKey) = value;
        end
    end
    
    % Display the new containers.Map objects
    disp("Left Map:");
    disp(leftMap);
    
    disp("Right Map:");
    disp(rightMap);
    output={leftMap,rightMap};
end
