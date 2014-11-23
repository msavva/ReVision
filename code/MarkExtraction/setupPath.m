clear path;

visionChartRoot = '';  % Path to the top-level folder.
sharedCache = '';  % Path to cache for intermediate results
dirstr = '';

rootstr = [pathsep, visionChartRoot];
dirstr = [dirstr, rootstr, '/pieAlgorithms'];
dirstr = [dirstr, rootstr, '/barAlgorithms'];
dirstr = [dirstr, rootstr, '/utils'];
dirstr = [dirstr, rootstr, '/utils/bilateralfiltering'];
dirstr = [dirstr, rootstr, '/utils/StanfordMATLABTools'];

path(pathdef);
path([dirstr, pathsep, path]);

clear dirstr rootstr visionChartRoot;

% Path to the root of the corpus of images. See bar_finder_batch.m and
% pie_finder_batch.m.
corpus_root = '';
