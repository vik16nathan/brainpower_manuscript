import containers.Map;
%Initialize global variables
nrow = 45006;
%nrow = 595209;
ncol = 361;
%total_t = 2;
num_v = nrow/3;
gm_res = 3.5; %grey matter resistivity in ohms*meters (from Bouda et al., 2021)
neur_vol_density = 40000*(1000^3); %40,000 neurons/mm^3 * (1000^3 mm^3/m^3) --> ISSUE - this is a big assumption
dendrite_area = 10^-12; %roughly 1 square micron in cross-sectional area
%estimate from Lennie et al., 2003 (whole-brain)

% 1. Get the x, y, and z currents at each source. → current dipole
%15k vertices
scaling_factor=1;
%curr_dev = in_bst_results('link|Subject01/S01_AEF_20131218_01_600Hz_notch/results_MN_MEG_KERNEL_230531_1033.mat|Subject01/S01_AEF_20131218_01_600Hz_notch/data_deviant_average_230528_1851.mat', 1); 
curr_sta = in_bst_results('link|Subject01/S01_AEF_20131218_01_600Hz_notch/results_MN_MEG_KERNEL_230531_1033.mat|Subject01/S01_AEF_20131218_01_600Hz_notch/data_standard_average_230528_1851.mat', 1);

%SNR 10 - observed larger currents
%curr_sta = in_bst_results('link|Subject01/S01_AEF_20131218_01_600Hz_notch/results_MN_MEG_KERNEL_230621_1509.mat|Subject01/S01_AEF_20131218_01_600Hz_notch/data_standard_average_230528_1851.mat', 1);
%200k vertices
%curr_sta = in_bst_results('link|Subject01/S01_AEF_20131218_01_600Hz_notch/results_MN_MEG_KERNEL_230620_1136.mat|Subject01/S01_AEF_20131218_01_600Hz_notch/data_standard_average_230528_1851.mat',1);
% Get x, y, z coordinates; define the number of rows and columns

%Beamformer
%scaling_factor = 8.0886*10^(-11);
%curr_sta = in_bst_results('link|Subject01/S01_AEF_20131218_01_600Hz_notch/results_PNAI_MEG_KERNEL_230621_1437.mat|Subject01/S01_AEF_20131218_01_600Hz_notch/data_standard_average_230528_1851.mat', 1);

%sLORETA 
%scaling_factor = 0.036561;
%curr_sta = in_bst_results('link|Subject01/S01_AEF_20131218_01_600Hz_notch/results_sLORETA_MEG_KERNEL_230621_1643.mat|Subject01/S01_AEF_20131218_01_600Hz_notch/data_standard_average_230528_1851.mat', 1);

%%%
%curr_dev_amp = curr_dev.ImageGridAmp;
curr_sta_amp = scaling_factor*curr_sta.ImageGridAmp;
%curr_dev_xyz = convertTo3DArray(curr_dev_amp, nrow, ncol);cke
curr_sta_xyz = convertTo3DArray(curr_sta_amp, nrow, ncol);

%Repeat with constrained currents - normal to cortex
curr_sta_norm_file = in_bst_results('link|Subject01/S01_AEF_20131218_01_600Hz_notch/results_MN_MEG_KERNEL_230528_2055.mat|Subject01/S01_AEF_20131218_01_600Hz_notch/data_standard_average_230528_1851.mat',1);
%Repeat for all 200k sources
%curr_sta_norm_file = in_bst_results('link|Subject01/S01_AEF_20131218_01_600Hz_notch/results_MN_MEG_KERNEL_230620_1055.mat|Subject01/S01_AEF_20131218_01_600Hz_notch/data_standard_average_230528_1851.mat',1);
curr_sta_norm = curr_sta_norm_file.ImageGridAmp;
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

%Create a dictionary of triangles attached to each vertex
v_t_dict = containers.Map('KeyType', 'int32', 'ValueType', 'any');
%cellWrapper = struct('data', {}); %allow insert into v_t_dict
for i = 1:size(triangles,1)
    %Iterate through each vertex of each triangle to add triangle to dict.
    for j = 1:3
        ver = triangles(i,j);
        if ~isKey(v_t_dict, ver)
           %cellWrapper.data = {triangles(i,:)};
           v_t_dict(ver) = {triangles(i,:)};
        else
           v_t_dict(ver) = [v_t_dict(ver), triangles(i,:)];
        end
    end
end

%For each vertex (source in our current source model), we want an estimate
%for cross-sectional area of a "wire" in the cortex with length equal to
%the cortical thickness. This wire represents tens of thousands of parallel
%neurons. See averageTriangleAreas function at bottom 

%Power estimate for each vertex at t=1... eventually repeat for more times
overall_P_unc = 0;
overall_P_con = 0;

support_areas = [];
num_neurons = [];
t=1;
resistances = [];
currents = [];
dip_powers = [];
for v=1:num_v

    thickness_in_m = thickness(v)/1000; %ASSUMPTION: displayed cortical thicknesses are in mm
    I_x = curr_sta_xyz(v,t,1)/thickness_in_m;
    I_y = curr_sta_xyz(v,t,2)/thickness_in_m;
    I_z = curr_sta_xyz(v,t,3)/thickness_in_m;

    %repeat for constrained source model
    I_norm = curr_sta_norm(v,t)/thickness_in_m;

    %Convert from resistivity to resistance
    %option 1 - use R = rho*L/A, where A is the area of the "support"
    %around each current dipole
    support_area = oneThirdSumAreas(v, vertices_xyz, v_t_dict);
    support_areas = [support_areas, support_area];
    R = gm_res*thickness_in_m/support_area;
    resistances = [resistances, R];
    %option 2 - treat resistance like current dipole and simply divide by
    %length (ignore cortical area) --> IGNORE THIS FOR NOW
    %R = gm_res/thickness_in_m;

    %option 3 - assume the cross-sectional area of a "cortical wire" is
    %simply that of an apical dendrite
    %cross sectional area of ~1 square micrometer (10^-12 square meters)
    %R = gm_res*thickness_in_m/dendrite_area;

    %When using this method and assuming summed neuronal currents,
    % a 1/N correction factor from dipole current --> summed neuronal
    % currents is needed (see calculations, ask Mark), where N is the
    % inferred number of neurons contributing to the current
    volume = support_area*thickness_in_m;
    num_neur = neur_vol_density*volume;
    num_neurons = [num_neurons, num_neur];
    correction_factor = 1; %--> use this for options 1 and 2
    %correction_factor = 1/num_neur; %use this for option 3
    %P_x = I_x*I_x*R;
    %P_y = I_y*I_y*R;
    %P_z = I_z*I_z*R;
    %overall_P_unc = overall_P_unc + sqrt(P_x^2 + P_y^2 + P_z^2);
    overall_P_con = overall_P_con + correction_factor*I_norm^2*R;
    I_equiv = sqrt(I_x^2+I_y^2+I_z^2);
    currents = [currents, I_equiv];
    overall_P_unc = overall_P_unc + correction_factor*I_equiv^2*R;
    dip_powers = [dip_powers, correction_factor*I_equiv^2*R];
end
%unc_power_vs_t = [unc_power_vs_t, overall_P_unc];
%con_power_vs_t = [con_power_vs_t, overall_P_con];

fprintf("Overall cortical power consumption (unconstrained source model) in W is: %s\n", num2str(mean(unc_power_vs_t)));
fprintf("Overall cortical power consumption (constrained source model) in W is: %s\n", num2str(mean(con_power_vs_t)));

%figure;
%subplot(2, 1, 1);
%plot(1:total_t, unc_power_vs_t);
%xlabel('Time');
%ylabel('Power (W)');
%title(sprintf('Unconstrained Power Consumption vs. Time'));

%subplot(2, 1, 2);
%plot(1:total_t, con_power_vs_t);
%xlabel('Time');
%ylabel('Power (W)');
%title(sprintf('Constrained Power Consumption vs. Time'));
    

%Other debugging: make sure currents are in reasonable units
%Convert current dipoles to absolute value and nA.m
curr_abs = abs(unique(curr_sta_xyz));
curr_abs_norm = abs(unique(curr_sta_norm));
curr_abs_1_10_nAm = find(curr_abs >= 1*10^-9 & curr_abs <= 10*10^-9);
curr_abs_norm_1_10_nAm = find(curr_abs_norm >=1*10^-9 & curr_abs_norm <= 10*10^-9);

curr_abs_1_10_pAm = find(curr_abs >= 1*10^-12 & curr_abs <= 10*10^-12);
curr_abs_norm_1_10_pAm = find(curr_abs_norm >=1*10^-12 & curr_abs_norm <= 10*10^-12);
fprintf("Proportion of xyz currents between 1-10 nA.m: %s\n",num2str(length(curr_abs_1_10_nAm)/length(curr_abs)));
fprintf("Proportion of constrained currents between 1-10 nA.m: %s\n", num2str(length(curr_abs_norm_1_10_nAm)/length(curr_abs_norm)));

fprintf("Proportion of xyz currents between 1-10 pA.m: %s\n",num2str(length(curr_abs_1_10_pAm)/length(curr_abs)));
fprintf("Proportion of constrained currents between 1-10 pA.m: %s\n", num2str(length(curr_abs_norm_1_10_pAm)/length(curr_abs_norm)));
% Histogram for support areas
support_areas_no_outliers = support_areas(support_areas < quantile(support_areas, 0.95));

figure;
histogram(support_areas_no_outliers);
xlabel('Area (m^2)');
ylabel('Number of dipoles');
title('Cross-sectional cortical column areas across dipoles');
saveas(gcf, 'support_areas_no_out.png');

% Histogram for resistances
resistances_no_outliers = resistances(resistances < quantile(resistances, 0.95));

figure;
histogram(resistances_no_outliers);
xlabel('Resistance (ohms)');
ylabel('Number of dipoles');
title('R=rho*L/A across cortical columns');
saveas(gcf, 'resistances_no_out.png');

% Histogram for currents
currents_no_outliers = currents(currents < quantile(currents, 0.95));

figure;
histogram(currents_no_outliers);
xlabel('Current (A)');
ylabel('Number of dipoles');
title('Current (current dipole/thickness) across cortical columns');
saveas(gcf, 'currents_no_out.png');

% Histogram for dipole powers
dip_powers_no_outliers = dip_powers(dip_powers < quantile(dip_powers, 0.95));

figure;
histogram(dip_powers_no_outliers);
xlabel('Power (W)');
ylabel('Number of dipoles');
title('Power across cortical columns');
saveas(gcf, 'powers_no_out.png');

% Histogram for thickness
thickness_no_outliers = thickness(thickness < quantile(thickness, 0.95));

figure;
histogram(thickness_no_outliers);
xlabel('Cortical thickness (mm)');
ylabel('Number of dipoles');
title('Thickness across cortical columns');
saveas(gcf, 'thickness_no_out.png');

%hist(num_neurons);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Wrangle current dipoles extracted from brainstorm into a single 3D matrix
%with three-nested 2D matrices for current * time in the x, y, and z
%dimensions


function output = convertTo3DArray(matrix, nrow, ncol)  
    % Initialize the 3D array
    output = zeros(nrow/3, ncol, 3);
    
    % Convert the matrix to as 3D array
    for n = 0:nrow/3-1
        output(n+1,:,1) = matrix(3*n+1,:,:);
        output(n+1,:,2) = matrix(3*n+2,:,:);
        output(n+1,:,3) = matrix(3*n+3,:,:);
    end
end

%Area of a triangle given coordinates of 3D points - sanity check for
%heronsFormula
function area = calculateTriangleArea3D(x1, y1, z1, x2, y2, z2, x3, y3, z3)
    v1 = [x2-x1, y2-y1, z2-z1];
    v2 = [x3-x1, y3-y1, z3-z1];
    area = 0.5 * norm(cross(v1, v2));
end

function area = heronsFormula(v1_xyz, v2_xyz, v3_xyz)
    d_12 = pdist2(v1_xyz, v2_xyz);
    d_23 = pdist2(v2_xyz, v3_xyz);
    d_13 = pdist2(v1_xyz, v3_xyz);
    s = (d_12 + d_23 + d_13)/2;
    area = sqrt(s*(s-d_12)*(s-d_23)*(s-d_13));
end

%Estimate the cross-sectional area of a cortical wire by looking at the
%triangles attached to a dipole
function avg_area = oneThirdSumAreas(v, vertices_xyz, v_t_dict)
    attached_tris = v_t_dict(v);
    areas = [];
    for i=1:length(attached_tris) %find vertices of all triangles attached to a vertex
        v1_xyz = vertices_xyz(attached_tris{i}(1),:);
        v2_xyz = vertices_xyz(attached_tris{i}(2),:);
        v3_xyz = vertices_xyz(attached_tris{i}(3),:);
        area = heronsFormula(v1_xyz, v2_xyz, v3_xyz);
        areas = [areas, area];
    end
    avg_area = 1/3*sum(areas);
end

