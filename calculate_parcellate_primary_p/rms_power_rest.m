%GOAL: Calculate RMS of 3D dipole currents across 120 secs

num_seconds = 120;
nrow=45006;
ncol=2400;
sub_fem_fwd_fnames={'sub-0002/sub-0002_ses-01_task-rest_run-01_meg_notch_high/results_MN_MEG_KERNEL_240118_1840.mat',
'sub-0004/sub-0004_ses-01_task-rest_run-01_meg_notch_high/results_MN_MEG_KERNEL_240306_1639.mat',
'sub-0006/sub-0006_ses-01_task-rest_run-01_meg_notch_high/results_MN_MEG_KERNEL_240226_1120.mat',
'sub-0007/sub-0007_ses-01_task-rest_run-01_meg_notch_high/results_MN_MEG_KERNEL_240306_1641.mat'};
SubjectNames = {'sub-0002', ...
    'sub-0004', 'sub-0006', 'sub-0007'};

for s=1:numel(SubjectNames) %run 5 operations at a time because 12 
    %nodes are available, but we don't want to overload memory
    sum_sq_idip=zeros(15002,3);
    subject=SubjectNames{s};
    fname=sub_fem_fwd_fnames{s};
    for j=1:num_seconds 
        numDigits=3;
        trial_num = sprintf('%0*d', numDigits, j);
        data = in_bst_results(sprintf(['link|%s|' ...
            '%s/%s_ses-01_task-rest_run-01_meg_notch_high/data_block%s_02.mat'], ...
            fname, subject, subject, trial_num), 1);
    
        amp = data.ImageGridAmp;

        curr_xyz = convertTo3DArray(amp, nrow, ncol);
        ss_curr=squeeze(sum(curr_xyz.^2,2));
        sum_sq_idip=sum_sq_idip+ss_curr;
        %disp(size(ss_curr));
        %disp(mean(ss_curr));
        %disp(size(mean(curr_xyz,2)));

    end
    ms_idip = sum_sq_idip./(ncol*num_seconds);
    
end

