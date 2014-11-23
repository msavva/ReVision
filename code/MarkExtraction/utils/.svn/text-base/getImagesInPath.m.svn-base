function [chart_files, fullpath, outpath] = getImagesInPath(path)
    % Input/output A = filepaths
    chart_path = path;
    %chart_path = 'Testing/VisInterpCorpus';
    fullpath = fullfile(pwd, chart_path);
    outpath = fullfile(fullpath,'output');

    % Retrieve file names for test charts
    all_image_path = dir(fullpath);

    chart_files = [];
    for i=1:length(all_image_path)
        if (~isempty(strfind(all_image_path(i).name, 'jpg')) || ...
            ~isempty(strfind(all_image_path(i).name, 'png')))
            chart_files = [chart_files; {all_image_path(i).name}];
        end
    end
end