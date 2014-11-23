%% Define magic numbers and thresholds used in the algorithms

%% Bilateral filter
global BFILTER_WIDTH BFILTER_DOMAINSIGMA BFILTER_RANGESIGMA
BFILTER_WIDTH = 2;
BFILTER_DOMAINSIGMA = 2;
BFILTER_RANGESIGMA = 50/255;

%% Edge detection

% Starting sigma
EDGE_SIGMA = 1.5;
% Starting edge strength
EDGE_STRENGTH = 2^8;

% Minimum edge strength
MIN_EDGE_STRENGTH = 0.01;

% Stop halving sigma when the number of edges reaches at least this
% proportion, or when the sigma is at the minimum sigma
MIN_EDGE_PROPORTION = 1/30;

% Flag to do non-max suppression in the edge detector
NON_MAX_SUPPRESSION = true;

%% Color detection
LAB_COLOR_THRESH = 20;

%% Curve tracing
% Opening element for the edgel map
curveTracing_closingElement = strel('square', 2);

% Maximum normal angle difference two neighboring pixels can have to be
% grouped together
curveTracing_maxDirDiff = pi/16;

%% Curve classification
lineAspectThreshold = 0.1;
textAspectThreshold = 0.5;

%% Pie fitting
%ELLIPSE_RANSAC_ITERATIONS = 500;
ELLIPSE_RANSAC_ITERATIONS = 100000;
%ELLIPSE_RANSAC_ITERATIONS = 10;

CIRCLE_RANSAC_ITERATIONS = 100000;

%% Pie line extraction
pie_center_tolerance = 15;
GOOD_LINE_THRESH = 2;
ELLIPSE_SUBDIVISIONS = 1000;
MIN_SEGMENT_LENGTH = 8;
MIN_VOTES = 4;

%% Bar extraction
MINIMUM_CANDIDATE_BARS = 2;

GOOD_BAR_THRESH = 3;

% Rectangle coverage
RECTANGLE_COVERAGE = 0.90;

% Width/height slack
DIM_SLACK = 3;

% Axis distance slack
AXIS_SLACK = 1;

% Rectangle color check travel distance
COLOR_CHECK_STEP = 5;

%% Text detection variables
text_centroid_distance = 10;
