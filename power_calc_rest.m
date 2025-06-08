import containers.Map;
%Initialize global variables
nrow = 45006;
%nrow = 595209;
ncol = 2400;
total_t = 1200; %sampling frequency
num_v = nrow/3;
gm_res = 3.5; %grey matter resistivity in ohms*meters (from Bouda et al., 2021)
output_dir='/export02/data/vikramn/brainstorm3/brainpower/resting_power/';
load('../thickness_areas_rest.mat'); %store matlab variables for support areas/
%cortical thickness attached to each vertex for EACH individual of this
%dataset (extend containers.Map object somewhat)

%Single run for each individual - loop across individuals instead
SubjectNames = {'sub-0002', 'sub-0003', ...
    'sub-0004', 'sub-0006', 'sub-0007'};

mn_run_times = {'1304','1304','1304','1305','1305'};
num_seconds = 100;


for j=1:num_seconds %run 5 operations at a time because 12 
    %nodes are available, but we don't want to overload memory
    xyz_currs = containers.Map();
    numDigits=3;
    trial_num = sprintf('%0*d', numDigits, j);

    for s=1:numel(SubjectNames)
        run_time = mn_run_times{s};            
        subject=SubjectNames{s};
        thickness_in_m = subject_thicknesses(subject);
        support_areas = subject_support_areas(subject);
        data = in_bst_results(sprintf(['link|%s/%s_ses-01_task-rest_run-01_' ...
            'meg_notch_high/results_MN_MEG_KERNEL_230801_%s.mat|' ...
            '%s/%s_ses-01_task-rest_run-01_meg_notch_high/data_block%s_02.mat'], ...
            subject, subject, run_time, subject, subject, trial_num), 1);
    
        amp = data.ImageGridAmp;
        curr_xyz = convertTo3DArray(amp, nrow, ncol);
        curr_xyz_downsample = extractEveryOtherColumn(curr_xyz);
        xyz_currs(subject) = curr_xyz_downsample;
        
    end

    parfor s=1:numel(SubjectNames)
        subject=SubjectNames{s};
        curr_xyz_downsample = xyz_currs(subject); %look into this
        
        power_vs_t_all = calculatePower(curr_xyz_downsample, thickness_in_m, support_areas, gm_res, num_v, total_t);
        %cell array with whole-brain time series + dipole-resolved time
        %series
        dip_p = power_vs_t_all{2};
        dip_i = vecnorm(curr_xyz, 2, 3);
        parsave(sprintf('%s%s_%s_dip_p.mat', output_dir, subject, trial_num), dip_p);
        parsave(sprintf('%s%s_%s_dip_i.mat', output_dir, subject, trial_num), dip_i);
           
    end
end
%%%%%%%%%%HELPER FUNCTIONS*********************************
function resultMatrix = extractEveryOtherColumn(inputMatrix) %3D matrix 
    resultMatrix = inputMatrix(:, 1:2:end,:);
end
