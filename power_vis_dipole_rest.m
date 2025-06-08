%Need to visualize separately for each subject (or combine onto reference
%anatomy...)

%But average across time
num_v=15002;
SubjectNames = {'sub-0002', 'sub-0003', ...
    'sub-0004', 'sub-0006', 'sub-0007'};

%Repeat for a single individual
%SubjectNames = {'sub-0002'};
input_dir = '../resting_power/';
load('../thickness_areas_rest.mat');
% for s=1:numel(SubjectNames)
%     subject=SubjectNames{s};
%     dip_avg_sec_p = findAndAverageFiles(input_dir, subject, 'p');
%     dip_avg_sec_i = findAndAverageFiles(input_dir, subject, 'i');
%     save(sprintf('%s_dip_p_i_avg.mat', subject), 'dip_avg_sec_p', 'dip_avg_sec_i');
% end
%Should have saved in schaeffer_rest_avg_var_p_i.m

%Load input data - take averages across subjects
total_current = zeros(num_v, 1);
total_power = zeros(num_v, 1);
total_cort_thick = zeros(num_v, 1);
total_supp_areas = zeros(num_v, 1);

%Create a dictionary of avg. power/current for each subject
sub_avg_p_dict = containers.Map();
sub_avg_i_dict = containers.Map();
sub_ct_dict = containers.Map();
sub_sa_dict = containers.Map();

for s=1:numel(SubjectNames)
    subject=SubjectNames{s};
    dip_avg_file = load(sprintf('%s_dip_p_i_avg.mat', subject));
    dip_avg_sec_p = cell2mat(dip_avg_file.dip_avg_sec_p);
    dip_avg_p = mean(dip_avg_sec_p,2);
    sub_avg_p_dict(subject) = dip_avg_p;
    
    dip_avg_sec_i = cell2mat(dip_avg_file.dip_avg_sec_i);
    dip_avg_i = mean(dip_avg_sec_i,2);
    sub_avg_i_dict(subject) = dip_avg_i;

    %total_current = total_current + dip_avg_i;
    %total_power = total_power + dip_avg_p;

    thickness_in_m = subject_thicknesses(subject);
    sub_ct_dict(subject) = thickness_in_m;
    support_areas = subject_support_areas(subject);
    sub_sa_dict(subject) = support_areas;
    %total_cort_thick = total_cort_thick + thickness_in_m;
    %total_supp_areas = total_supp_areas + support_areas;
    
end

save('sub_p_i_ct_sa_dicts.mat', 'sub_avg_p_dict', 'sub_avg_i_dict', ...
    'sub_ct_dict', 'sub_sa_dict');

%Manually project average power and current across all 100s for each individual
% onto default anatomy

%Repeat for thickness and area

%Start with current - load in variables for min norm kernels for each
%subject called sub02, sub03, ... 
load('sub_p_i_ct_sa_dicts.mat');

%Repeat with power, support area, cortical thickness (change the
%dictionary)


mean_i = total_current/numel(SubjectNames);
mean_p = total_power/numel(SubjectNames);
mean_thick = total_cort_thick/numel(SubjectNames);
mean_sa = total_supp_areas/numel(SubjectNames);

%Fit a quadratic regression between power and current
%Find the residuals, plot the curve
x=mean_i;
y=mean_p;
coefficients = polyfit(x, y, 2); 

y_fit = polyval(coefficients, x);
resid = y-y_fit;
x_pf = linspace(min(x), 0.8*max(x));
y_pf = polyval(coefficients, x_pf);
figure;
set(gcf, 'Position', get(0, 'Screensize'));
scatter(x, y);
hold on;
plot(x_pf, y_pf, 'r');
r2 = calculate_r2(y, y_fit);
xlabel('Mean Current Dipole Across Individuals (A/m)')
ylabel('Mean Dipole Power Across Individuals (W)')
title('Average Power vs. Current for n=5 Resting Omega Individuals')
text(0.8*max(x), 0.8*max(y), ['R squared: ' num2str(r2)]);
%For resting state scan

data_table = table(mean_i, mean_p, mean_thick, mean_sa, resid, 'VariableNames', ...
    {'mean_i', 'mean_p', 'mean_thick', 'mean_sa', 'resid'});

filename = 'omega_avg_data.csv';
writetable(data_table, filename);
save('omega_mean_p_i_thick_area.mat', 'mean_i', 'mean_p', 'mean_thick', 'mean_sa');
%Parcellate the residuals, mean cortical thicknesses, and support areas
%using default anatomy... figure this out later

% Set the plot to be grid-like
grid on;
saveas(gcf, 'p_v_i_avg_omega.png');
% Release the hold on the plot
hold off;

%Visualize residuals (setting outlier threshold)
resid_out = filter_outliers(resid);
resid_filt = resid_out{1};
%Export dummy time series as time_series
num_freqs = 181;
time_series.TF = repmat(resid_filt, 1, 1, num_freqs);


%Correlate residual with cortical thickness
x=1 ./ mean_thick; %1/L relationship based on formula
y=resid;
coefficients = polyfit(x, y, 1); %switch to line of best fit
%Should be a 1/L relationship - inverse

y_fit = polyval(coefficients, x);
resid = y-y_fit;
x_pf = linspace(min(x), 0.8*max(x));
y_pf = polyval(coefficients, x_pf);
scatter(x, y);
hold on;
plot(x_pf, y_pf, 'r');
r2 = calculate_r2(y, y_fit);
xlabel('Inverse cortical thickness (1/m)');
ylabel('Residual P vs. I ');
title('P-I Deviation vs. Cortical Thickness');
text(0.8*max(x), 0.8*max(y), ['R squared: ' num2str(r2)]);
grid on;
saveas(gcf, 'resid_v_thick_avg.png');
hold off;

%Correlate residuals with support areas
x=1 ./ mean_sa;
coefficients = polyfit(x, y, 1); %switch to line of best fit
%Should be a 1/A relationship - inverse

y_fit = polyval(coefficients, x);
resid = y-y_fit;
x_pf = linspace(min(x), 0.8*max(x));
y_pf = polyval(coefficients, x_pf);
scatter(x, y);
hold on;
plot(x_pf, y_pf, 'r');
r2 = calculate_r2(y, y_fit);
xlabel('Inverse Support Area (1/m^2)')
ylabel('Residual P vs. I (no outliers)')
title('P-I Deviation vs. Support Area');
text(0.8*max(x), 0.8*max(y), ['R squared: ' num2str(r2)]);
grid on;
saveas(gcf, 'resid_v_area_avg.png')
hold off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Repeat without outliers for power
output = filter_outliers(mean_p);
y=output{1}; %Replace all outliers with max/min value
out_ind=output{2};
non_out_ind=output{3};
x = mean_i(cell2mat(non_out_ind));
y=y(cell2mat(non_out_ind));
fprintf('Number of outlier dipoles: %i', length(cell2mat(out_ind)));
coefficients = polyfit(x, y, 2); 

y_fit = polyval(coefficients, x);
resid = y-y_fit;
x_pf = linspace(min(x), 0.8*max(x));
y_pf = polyval(coefficients, x_pf);
scatter(x, y);
hold on;
plot(x_pf, y_pf, 'r');
r2 = calculate_r2(y, y_fit);
xlabel('Mean Current Dipole Across Individuals (A/m)')
ylabel('Mean Dipole Power Across Individuals (W)')
title('Average Power vs. Current for n=5 Resting (No Outlier Power)')
text(0.8*max(x), 0.8*max(y), ['R squared: ' num2str(r2)]);
grid on;
saveas(gcf, 'p_v_i_avg_omega_no_out.png');
hold off;

%Correlate residual with cortical thickness (no outliers)
thick_filt = mean_thick(cell2mat(non_out_ind));
x=1 ./ thick_filt; %1/L relationship based on formula
y=resid;
coefficients = polyfit(x, y, 1); %switch to line of best fit
%Should be a 1/L relationship - inverse

y_fit = polyval(coefficients, x);
resid = y-y_fit;
x_pf = linspace(min(x), 0.8*max(x));
y_pf = polyval(coefficients, x_pf);
scatter(x, y);
hold on;
plot(x_pf, y_pf, 'r');
r2 = calculate_r2(y, y_fit);
xlabel('Inverse cortical thickness (1/m)')
ylabel('Residual P vs. I (no outliers)')
title('P-I Deviation vs. Cortical Thickness');
text(0.8*max(x), 0.8*max(y), ['R squared: ' num2str(r2)]);
grid on;
saveas(gcf, 'resid_v_thick_no_out_avg.png');
hold off;

%Correlate residuals with support areas (no outliers)
area_filt = mean_sa(cell2mat(non_out_ind));
x=1 ./ area_filt;
coefficients = polyfit(x, y, 1); %switch to line of best fit
%Should be a 1/A relationship - inverse

y_fit = polyval(coefficients, x);
resid = y-y_fit;
x_pf = linspace(min(x), 0.8*max(x));
y_pf = polyval(coefficients, x_pf);
scatter(x, y);
hold on;
plot(x_pf, y_pf, 'r');
r2 = calculate_r2(y, y_fit);
xlabel('Inverse Support Area (1/m^2)')
ylabel('Residual P vs. I (no outliers)')
title('P-I Deviation vs. Support Area');
text(0.8*max(x), 0.8*max(y), ['R squared: ' num2str(r2)]);
grid on;
saveas(gcf, 'resid_v_area_no_out_avg.png');
hold off;

%Export a 'dummy' TS map for for power/current

%Parcellate everything based on default anatomy 

