function [data, ws, ps] = ufl_load(opt)
%% Load data
data = DataSet(opt);
ws = data.getWS();

%% Extract random patches
ps = PatchSet(opt, ws.trainX, ws.trainY);

%% Run K-means to obtain centroid patches
ps.kmeansCentroids();
% ps.kmeansCentroids2Tier();
ps.showCentroids();

end
