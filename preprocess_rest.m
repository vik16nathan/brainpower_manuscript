num_v = 15002;
SubjectNames = {'sub-0002', 'sub-0003', ...
    'sub-0004', 'sub-0006', 'sub-0007'};

%SubjectNames = {'sub-0003', ...
   %'sub-0004', 'sub-0006', 'sub-0007'};

subject_v_area_dicts = containers.Map();
subject_support_areas = containers.Map();
subject_thicknesses = containers.Map(); %after converting to m
%Load in PIAL 15k surface files as variables

subject_cat12_time = containers.Map();
subject_cat12_time('sub-0003') = '230727_1840';
subject_cat12_time('sub-0004') = '230727_2000';
subject_cat12_time('sub-0006') = '230727_2133';
subject_cat12_time('sub-0007') = '230727_2255';
subject_cat12_time('sub-0002') = '230728_1747'; %fill in after

for s=1:numel(SubjectNames)
    subject=SubjectNames{s};
    cat12_time = subject_cat12_time(subject);
    surface_file = in_bst_data(sprintf(['/export02/data/vikramn/TutorialOmega/' ...
        'anat/%s/tess_cortex_pial_low.mat'], subject));
    v_area_dict = create_tri_area_dict(surface_file.Faces, surface_file.Vertices);
    subject_v_area_dicts(subject) = v_area_dict;
    subject_support_areas(subject) = oneThirdSumAreas(num_v, v_area_dict);
    thickness = in_bst_data(sprintf(['/export02/data/vikramn/TutorialOmega/data/%s/' ...
        'CAT12/results_surface_thickness_%s_15002V.mat'], subject, cat12_time));
    thickness_in_m = thickness.ImageGridAmp(:,1) ./ 1000;
    subject_thicknesses(subject) = thickness_in_m;
   
end

save('thickness_areas_rest.mat', 'subject_support_areas','subject_v_area_dicts', 'subject_thicknesses');





