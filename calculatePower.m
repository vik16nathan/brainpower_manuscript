
function outputCellArray=calculatePower(curr_xyz, thickness_in_m, support_areas, gm_res, ...
    num_v, total_t)
%Assume we're only doing cortical columns (no single-neuron divided current
%approach with neuron volume density - see power_calc.m)
    I_x = curr_xyz(:, :, 1) ./ thickness_in_m;
    I_y = curr_xyz(:, :, 2) ./ thickness_in_m;
    I_z = curr_xyz(:, :, 3) ./ thickness_in_m;
   
    R = gm_res * thickness_in_m ./ support_areas;
    
    correction_factor = ones(num_v, total_t); % Initialize with ones for options 1 and 2
    %(see power_calc.m)
    
    I_equiv = sqrt(I_x.^2 + I_y.^2 + I_z.^2);
    power_dip_mat_full = correction_factor .* I_equiv.^2 .* R;
    power_vs_t = sum(power_dip_mat_full, 1);
    outputCellArray={};
    outputCellArray{1} = power_vs_t;
    outputCellArray{2} = power_dip_mat_full;

end

  

