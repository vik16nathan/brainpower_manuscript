%%%%%%%%%%%%%Calculate primary dipole power (as
%%%%%%%%%%%%%earlier)%%%%%%%%%%%%%%%%%%%%
%change directories to ./brainpower/rest
subject_list = {'sub-0002','sub-0004','sub-0006', 'sub-0007'};
%BEFORE in Brainstorm: Remap thickness and support areas to 15002 vertices 
% in CENTRAL SURFACE and not pial surface b/c FEM source estimation is
% w.r.t central surface (read DUNEuro documentation??). See
% preprocess_omega_fem.m (inspired by preprocess_rest.m)
nrow = 15002;
ncol = 3;
total_t = 1; %sampling frequency
num_v = 15002;
gm_res = 3; %grey matter resistivity in ohms*meters (from Bouda et al., 2021)
output_dir='/export02/data/vikramn/duneuro_malte/';
sub_dendr_r = {
      '/export02/data/vikramn/brainstorm3/resistances/sub-0002_region_res_cell_column.mat',
      '/export02/data/vikramn/brainstorm3/resistances/sub-0004_region_res_cell_column.mat',
      '/export02/data/vikramn/brainstorm3/resistances/sub-0006_region_res_cell_column.mat',
      '/export02/data/vikramn/brainstorm3/resistances/sub-0007_region_res_cell_column.mat'
    };
%Remap thickness and support areas to 15002 vertices in FEM head model
load('../../thickness_areas_omega_fem.mat'); %load cortical thickness 
for s=1:numel(subject_list)
    subject=string(subject_list(s));
    
    thickness_in_m = subject_thicknesses(subject);
    support_areas = subject_support_areas(subject);
    
    %ignore spherical head model for now
    %use three-layer FEM head model
    curr_xyz = load(sprintf("%s_idip_rms.csv", subject));
    %curr_xyz=load(sprintf("%s_idip_rms.csv",subject));
    %3D matrix is normally 15002 x t x 3 (where t is recording, and t/2400
    %is time in s)
    curr_xyz = reshape(curr_xyz, [nrow,1,ncol]); %3D matrix (2D matrices across time) for calculatePower
    %R = gm_res * thickness_in_m ./ support_areas;
    
    power_vs_t_all = calculatePower(curr_xyz, R, num_v, total_t);
    dip_p = power_vs_t_all{2};

    dip_i = vecnorm(curr_xyz, 2, 3);

    save(sprintf('%s%s_dip_p_rms.mat', output_dir, subject), "dip_p"); %in W
    %save(sprintf('%s%s_dip_i_0s.mat', output_dir, subject), "dip_i"); %in A.m

    %save(sprintf('%s_dip_p_rms.mat', subject), "dip_p"); %in W
    %save(sprintf('%s_dip_i_rms.mat', subject), "dip_i"); %in A.m
    
end
%%%%%%%%Correlate this with secondary power for sub-0002 (calculated for each
%%%%%%%%dipole)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%FIXED PIAL VS. CENTRAL SURFACE ISSUE FOR FEM VS DIP P!! 
%change slightly for RMS vs. 0 seconds
for s=1:numel(subject_list)
    subject=subject_list{s};
    fprintf("Subject %s\n",subject)
    %sub_pri_p=load(sprintf('%s_dip_p_rms.mat',subject)).dip_p;
    %sub_sec_p=load(sprintf("%s%s_dip_fem_p_rms.csv",output_dir,subject));
    sub_pri_p=load(sprintf('%s_dip_p_0s.mat',subject)).dip_p;
    sub_sec_p=load(sprintf("%s%s_fem_p_0s.csv",output_dir,subject));

    fprintf("Correlation between primary and secondary P (0s I_{dip}): %s\n", num2str(corr(sub_pri_p, sub_sec_p)));
    fprintf("Method 1 Summed P_{sec} for Sub %s: %s\n", subject,num2str(sum(sub_sec_p)));
    fprintf("**************************************************\n")
end


%%%%%%%%Compare the net amplitudes of primary vs. secondary (for all inds)%
for s=1:numel(subject_list)
    sub=string(subject_list(s));
    %subp = load(sprintf("%s_dip_p_rms.mat", sub));
    subp = load(sprintf("%s_dip_p_0s.mat",sub));
    fprintf("Total primary power for %s: %s\n", sub, num2str(sum(subp.dip_p)));
end
%%%%%%%%%%%%%OLD: Power consumption averaged over one second%%%%%%%%%%%%%%%
%dip_p_1s_avg=load('./brainpower/rest/sub-0002_dip_p_1s.mat'); %in W
%dip_i_1s_avg=load('./brainpower/rest/sub-0002_dip_i_1s.mat'); %in A.m

%mean_dip_p = mean(dip_p,2);
%mean_dip_i = mean(dip_i, 2);
%Load in FEM power 
%fem_p=csvread('/export02/data/vikramn/duneuro_malte/sub0002_fem_p.csv');
%fem_i=csvread('/export02/data/vikramn/duneuro_bst_inputs/sub0002_mean_idip_1s.csv');
%fem_i=vecnorm(fem_i, 2,2);

%Summary statistics
%sum(fem_p);
