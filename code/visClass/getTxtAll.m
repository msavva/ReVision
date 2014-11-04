%% Set Range

clear;
setOpt;

rfs = 6;%5:1:20;
DIMS = 128;%[64 128 256];

rfd = combvec(rfs,DIMS);

for i=1:size(rfd,2)
    %% Get iteration variables
    rf = rfd(1,i);
    d = rfd(2,i);
    setName = strcat('rf',num2str(rf),'_d',num2str(d));
    opt.rfSize = rf;
    opt.DIM = [d d 1];
    
    %% Load DataSet
    data = DataSet(opt);
    ws = data.getWSFull();
    tic;
    %% Learn image features
    ps = PatchSet(opt, ws.X, ws.Y);
    ps.kmeansCentroids();
%     ps.showCentroids();
    
    % Compute feature vectors
    fs = ps.extractFeatures(ws.X,1);
    featureNames = strnum('centroid',4*opt.Ncentroids);
    
    if (opt.useTextFeatures)
        %% Get text features
        tws = data.getTextFeatureSetVectors(3);
        
        %% Normalize text features and concatenate with image feature vectors
        tXmean = mean(tws.X);
        tXsd = sqrt(var(tws.X)+0.01);
        tXs = bsxfun(@rdivide, bsxfun(@minus, tws.X, tXmean), tXsd);
        % tXs = tws.X;
        
        X = [fs.XCs(:,1:end-1) tXs];
        Y = ws.Y;
        featureNames = [featureNames; data.cat{1}.txt(1).featureVectorLabels(3)];
    else
        X = fs.XCs(:,1:end-1);
        Y = ws.Y;
    end
    
    %% Add category indicator binary variables
%     catIndicators = bsxfun(@eq,ws.Y,1:data.Ncategories);
%     catIndicatorNames = strcat('is_',opt.catDirs);
%     X = [X catIndicators];
%     featureNames = [featureNames; catIndicatorNames'];
    %% Add filename tags
%     cats = [data.cat{:}];
%     txts = [cats.txt]';
%     fullnames = char(txts.fname);
%     fnames = zeros(data.NtotalImages,1);
%     for j=1:data.NtotalImages
%         fnames(j) = filenameToIndex(fullnames(j,:));
%     end
%     X = [X fnames];
%     featureNames = [featureNames; 'filenameIndex'];

% Add is3D tag
if (opt.addIs3Dtag)
    featureNames = [featureNames; 'is3D'];
    X = [X ~cellfun('isempty',strfind(ws.C,'3D'))];
end

    %% Output
    if (opt.useTextFeatures)
        prefix = '../../feats_txt_';
    else 
        prefix = '../../feats_';
    end
    filename = strcat(prefix, setName, '.arff');
    data.featuresToARFF(filename, X, Y, setName, featureNames);
    fprintf('Saved to %s\n', filename);
end

totalTime = toc;