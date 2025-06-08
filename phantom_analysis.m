%%% DIPOLE SCANNING %%%

%Elekta phantom analysis
amplitudes={elekta_dipole_scan.Dipole.Amplitude};
goodness={elekta_dipole_scan.Dipole.Goodness};
amp_vec=[];
num_dipoles=32;
for i=1:num_dipoles
  amp_vec = [amp_vec, norm(amplitudes{i})];
end

disp(mean(amp_vec)/10^-9);
disp(std(amp_vec/10^-9)/sqrt(length(amp_vec)));
disp(mean(cell2mat(goodness)));

%%% MIN NORM %%%

%CTF

ctf_mn_200 = in_bst_results(['link|PhantomCTF-ds/phantom_200uA_20150709_01/' ...
    ['results_MN_MEG_KERNEL_231018_1231.mat|PhantomCTF-ds/'] ...
    ['phantom_200uA_20150709_01/data_stim_average_231014_1930.mat']], 1);

%time is zero at column 43 after epoching
%t0_col = 43; 
tmin_col=62;
ctf_mn_200_amp = ctf_mn_200.ImageGridAmp(:,tmin_col);
ctf_mn_200_total_idip = calculateTotalIDip(ctf_mn_200_amp);
disp(ctf_mn_200_total_idip);

%Repeat for 20 microA
t0_col=169;
tmin_col=243;
ctf_mn_20 = in_bst_results(['link|PhantomCTF-ds/phantom_20uA_20150603_03/' ...
    'results_MN_MEG_KERNEL_231018_1231.mat|PhantomCTF-ds/' ...
    'phantom_20uA_20150603_03/data_stim_average_231014_1930.mat'], 1);
ctf_mn_20_amp = ctf_mn_20.ImageGridAmp(:,tmin_col);
ctf_mn_20_total_idip = calculateTotalIDip(ctf_mn_20_amp);
disp(ctf_mn_20_total_idip);

%ELEKTA (load diff protocol... don't run at same time as lines above)

elekta_idips = [];
t60_col = 8; %column corresponding to 60 ms "peak" in Elekta dipole activity
for i=1:num_dipoles
    dip_num = sprintf('%02d', i);
    elekta_dip = in_bst_results(sprintf(['link|Kojak/kojak_all_200nAm_pp_no_chpi_no_ms_raw/' ...
        'results_MN_MEG_GRAD_KERNEL_231020_1214.mat|' ...
        'Kojak/kojak_all_200nAm_pp_no_chpi_no_ms_raw/' ...
        'data_%s_average_231016_1441.mat'], dip_num), 1);
    elekta_amp = elekta_dip.ImageGridAmp(:, t60_col);
    elekta_total_idip = calculateTotalIDip(elekta_amp);
    elekta_idips = [elekta_idips, elekta_total_idip];

end

disp(mean(elekta_idips));
disp(std(elekta_idips)/sqrt(32));

%Find the time where the dipole is the lowest and see if 