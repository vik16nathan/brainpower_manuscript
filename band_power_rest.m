import containers.Map;
%Initialize global variables
nrow = 45006;
%nrow = 595209;
ncol = 2400;
total_t = 1200; %sampling frequency
num_v = nrow/3;
gm_res = 3.5; %grey matter resistivity in ohms*meters (from Bouda et al., 2021)
output_dir='/export02/data/vikramn/brainstorm3/brainpower/resting_power/bands/';
load('../thickness_areas_rest.mat'); %store matlab variables for support areas/
%cortical thickness attached to each vertex for EACH individual of this
%dataset (extend containers.Map object somewhat)

%Single run for each individual - loop across individuals instead

%Manually update this - there has to be a better way
subject_band_runtimes = containers.Map();
subject_band_runtimes(['sub-0002','delta']) = '1212';
subject_band_runtimes(['sub-0002','theta']) = '1212';
subject_band_runtimes(['sub-0002','alpha']) = '1213';
subject_band_runtimes(['sub-0002','beta']) = '1213';
subject_band_runtimes(['sub-0002','gamma1']) = '1213';
subject_band_runtimes(['sub-0002','gamma2']) = '1214';

subject_band_runtimes(['sub-0003','delta']) = '1218';
subject_band_runtimes(['sub-0003','theta']) = '1218';
subject_band_runtimes(['sub-0003','alpha']) = '1219';
subject_band_runtimes(['sub-0003','beta']) = '1219';
subject_band_runtimes(['sub-0003','gamma1']) = '1219';
subject_band_runtimes(['sub-0003','gamma2']) = '1219';

subject_band_runtimes(['sub-0004','delta']) = '1219';
subject_band_runtimes(['sub-0004','theta']) = '1219';
subject_band_runtimes(['sub-0004','alpha']) = '1219';
subject_band_runtimes(['sub-0004','beta']) = '1219';
subject_band_runtimes(['sub-0004','gamma1']) = '1219';
subject_band_runtimes(['sub-0004','gamma2']) = '1220';

subject_band_runtimes(['sub-0006','delta']) = '1220';
subject_band_runtimes(['sub-0006','theta']) = '1220';
subject_band_runtimes(['sub-0006','alpha']) = '1220';
subject_band_runtimes(['sub-0006','beta']) = '1220';
subject_band_runtimes(['sub-0006','gamma1']) = '1220';
subject_band_runtimes(['sub-0006','gamma2']) = '1220';


subject_band_runtimes(['sub-0007','delta']) = '1220';
subject_band_runtimes(['sub-0007','theta']) = '1220';
subject_band_runtimes(['sub-0007','alpha']) = '1220';
subject_band_runtimes(['sub-0007','beta']) = '1221';
subject_band_runtimes(['sub-0007','gamma1']) = '1221';
subject_band_runtimes(['sub-0007','gamma2']) = '1221';

SubjectNames = {'sub-0002', 'sub-0003', ...
    'sub-0004', 'sub-0006', 'sub-0007'};
freq_band_filenames = {'band', 'band_02', 'band_03', ...
    'band_04', 'band_05', 'band_06'};
freq_band_filenames = flip(freq_band_filenames);
freq_band_names = {'delta', 'theta', 'alpha', ...
    'beta', 'gamma1', 'gamma2'};
freq_band_names = flip(freq_band_names);

num_seconds = 100;
for j=1:num_seconds %run 5 operations at a time because 12 
    %nodes are available, but we don't want to overload memory
    xyz_currs = containers.Map();
    numDigits=3;
    trial_num = sprintf('%0*d', numDigits, j);
    
   for b=1:numel(freq_band_names)
        band_filename = freq_band_filenames{b};
        band_name = freq_band_names{b};
        for s=1:numel(SubjectNames)        
            subject=SubjectNames{s};
            run_time = subject_band_runtimes([subject, band_name]);  
            thickness_in_m = subject_thicknesses(subject);
            support_areas = subject_support_areas(subject);

            %link|sub-0007/sub-0007_ses-01_task-rest_run-01_meg_notch_high_band_05/
            %results_MN_MEG_KERNEL_230801_1221.mat|sub-0007/
            %sub-0007_ses-01_task-rest_run-01_meg_notch_high_band_05/data_block002.mat

            data = in_bst_results(sprintf(['link|%s/%s_ses-01_task-rest_run-01_' ...
                'meg_notch_high_%s/results_MN_MEG_KERNEL_230801_%s.mat|' ...
                '%s/%s_ses-01_task-rest_run-01_meg_notch_high_%s/data_block%s.mat'], ...
                subject, subject, band_filename, run_time, subject, subject, ...
                band_filename, trial_num), 1);
        
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
            parsave(sprintf('%s%s_%s_%s_dip_p.mat', output_dir, subject, band_name, trial_num), dip_p);
            parsave(sprintf('%s%s_%s_%s_dip_i.mat', output_dir, subject, band_name, trial_num), dip_i);
               
        end
    end
end
%%%%%%%%%%HELPER FUNCTIONS*********************************
function resultMatrix = extractEveryOtherColumn(inputMatrix) %3D matrix 
    resultMatrix = inputMatrix(:, 1:2:end,:);
end