function averageValues = findAndAverageFiles(directory, subject, band, p_or_i)
    % Initialize an empty cell array to store the average values
    averageValues = {};
    band_string='';
    if band == ""
        band_string="";
    else
        band_string=sprintf('%s_',band);
    end
    % Get a list of all files in the specified directory
    stringToMatch = sprintf('%s_%s*_dip_%s.mat', ...
        subject, band_string, p_or_i);
    files = dir(fullfile(directory, stringToMatch));
    
    % Loop through each file
    for i = 1:length(files)
        % Load the MAT file
        data = load(fullfile(directory, files(i).name));
        varInfo = whos('-file', fullfile(directory,files(i).name));
        varName = varInfo.name;
        tsValues = data.(varName);
        % Assuming the variable containing the values is named 'values'
        % Calculate the average of the values in the file
        avgValue = mean(tsValues, 2);
        
        % Append the average value to the array
        averageValues = [averageValues, avgValue];
    end
end

