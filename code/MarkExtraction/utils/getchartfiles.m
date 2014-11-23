function files = getchartfiles(subset_file)
    chartFile = fopen(fullfile(pwd, subset_file),'r');
    files = textscan(chartFile, '%s');
    files = files{1};
    fclose(chartFile);
end