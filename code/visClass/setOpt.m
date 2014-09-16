%% Load VL_FEAT
curDir = pwd;
cd C:\code\VisInterp\code\vlfeat-0.9.9\toolbox\
vl_setup('NOPREFIX');
cd(curDir);

%% Configuration
opt.rfSize = 6;
opt.Npatches = 50000;
opt.Ncentroids = 200;
opt.kmeansIterations = 50;
opt.whitening = true;
opt.Ntrain = 1;
opt.Ntest = 1;
opt.DIM=[128 128 1];
opt.randomizeImageOrder = false;
opt.preserveAR = true;
opt.binaryClassification = false;
opt.useTextFeatures = false;
opt.filterTextRegions = false;
opt.subdivisionLevels = 1;
opt.addIs3Dtag = false;

opt.baseDir = 'C:\code\VisInterp\images\';
opt.catDirs = {'bar', 'bar3D','lines', 'lines3D', 'pie', 'pie3D', 'scatter','scatter3D', 'surface3D'};
% opt.catDirs = {'bar','lines','pie', 'scatter'};
% opt.catDirs = {'bar','pie'};
% opt.baseDir = 'C:\code\VisInterp\imagesNew\';
% opt.catDirs = { 
%     'AreaGraph\', ...         %1
%     'BarGraph\', ...          %2
%     'BoxPlot\', ...           %3
%     'ColumnGraph\', ...       %4
%     'FlowChart\', ...         %5
%     'LineGraph\', ...         %6
%     'Map\', ...               %7
%     'NetworkDiagram\', ...    %8
%     'ParetoChart\', ...       %9
%     'PieChart\', ...          %10
%     'RadarPlot\', ...         %11
%     'ScatterGraph\', ...      %12
%     'Table\', ...             %13
%     'TreeDiagram\', ...       %14
%     'VennDiagram\'};          %15
% opt.catDirs = strcat(opt.catDirs,'jpg');