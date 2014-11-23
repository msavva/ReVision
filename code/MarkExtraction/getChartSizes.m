% Get the sizes of a set of chart images.

%% Retrieve the charts
clear all
close all

setupPath();
thresholds();

IMAGE_READ_DIRECTORY = 1;
IMAGE_READ_SUBSET = 2;

%%%%%%%%%%%%%%%%%%%%%%
% Flags

CURVE_OPTIMIZE = false;
IMAGE_READ_MODE = IMAGE_READ_DIRECTORY;

% [chart_files fullpath outpath] = getImagesInPath('Testing/UMDpies');
% textpath =
% '/Users/nick/Documents/Projects/VisionCharts/VisInterp/corpus/UMDImages/p
% ie';
% cachepath = fullfile(pwd, 'cache');
%

%% Use a text file that specifies the files
if IMAGE_READ_MODE == IMAGE_READ_SUBSET

    chart_subset_file = 'Testing/largepies/large_pies.txt';
    %chart_subset_file = 'Testing/UMDpies/umd_pies.txt';
    chart_path = 'Visualization/PieChart';
    %chart_path = 'UMDImages/pie';
    chart_files = getchartfiles(chart_subset_file);
    outpath = fullfile(pwd, 'Testing/largepies');
    %outpath = fullfile(pwd, 'Testing/UMDpies');

    for i=1:length(chart_files)
        chart_files{i} = strcat(chart_files{i}, '.png');
    end

elseif IMAGE_READ_MODE == IMAGE_READ_DIRECTORY

    %% Use a directory

    % Large corpus paths
    chart_path = 'Scraped/bing_pie chart/UIST11_assumptions';
    outpath = fullfile(pwd, 'Testing/bing_pie chart');

    fullpath = fullfile(corpus_root, chart_path);
    listing = dir(fullpath);
    chart_files = {};
    curInd = 1;
    for i=1:size(listing,1)
        if (~isempty(regexpi(listing(i).name, '.jpg')) || ...
            ~isempty(regexpi(listing(i).name, '.png')) || ...
            ~isempty(regexpi(listing(i).name, '.gif')))
            chart_files{curInd} = listing(i).name;
            curInd = curInd + 1;

        end
    end

end

fullpath = fullfile(corpus_root, chart_path);

%% Compute the sizes
sizes = zeros(length(chart_files),2);
for i=1:length(chart_files)
    curFile = chart_files{i};
    [curIm curImGray] = readimage(fullfile(fullpath, curFile));
    sizes(i,:) = size(curImGray);
end