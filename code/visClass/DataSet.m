classdef DataSet < handle
    
    properties
        baseDir;                            % Base Corpus Directory
        catDirs;                            % Cell array of category subdirs
        DIM = [128 128 1];                  % Image dimensions
        Ntrain = 30;                        % Number of training images per category
        Ntest = 20;                         % Number of test images per category
        randomizeImageOrder = false;        % Whether to randomize image order
        binaryClassification = -1;          % <= 0 implies multi-class, otherwise index to positive category
        preserveAR = true;                  % Whether to preserve Aspect Ratio when rescaling input images
        useTextFeatures = false;            % Whether to use extracted text features
        filterTextRegions = false;          % Whether to filter out text regions from image data
        
        cat;                                % Cell array of categories
        maxCatSize = 0;                     % Maximum category size
        NtotalImages = 0;                   % Number of total images
        Ncategories = 0;                    % Number of categories
    end % Instance properties
    
    properties (Dependent = true)
        trainX;                             % [Ntrain*Ncategories X Npixels] training data
        trainY;                             % [Ntrain*Ncategories X 1] training category labels
        testX;                              % [Ntest*Ncategories X Npixels] testing data
        testY;                              % [Ntest*Ncategories X 1] testing category labels
    end % Derived properties
    
    properties (Hidden = true, Access = private)
        loaded = false          % Images loaded flag
        mask;                   % [nCategories X 2] random permutation mask cell array for pulling out train and test subset indices (1st and 2nd correspondingly)
    end
    

    methods
        % Save a labeled set of feature vectors into Weka format ARFF files
        function wekaOBJ = featuresToARFF(ds, filename, X, Y, setName, featureNames)
            addpath('./wekaInterface');
            javaaddpath('./wekaInterface/weka.jar');
            
            featureNames = [featureNames ; 'class'];
            data = [num2cell(X) ds.catDirs(Y)'];

            wekaOBJ = matlab2weka(setName, featureNames, data);
            saveARFF(filename, wekaOBJ);
        end

        % Constructor
        function d = DataSet(opt)
            if nargin > 0
                % Pull variables
                d.baseDir = opt.baseDir;
                d.catDirs = opt.catDirs;
                d.DIM = opt.DIM;
                d.Ntrain = opt.Ntrain;
                d.Ntest = opt.Ntest;
                d.randomizeImageOrder = opt.randomizeImageOrder;
                d.binaryClassification = opt.binaryClassification;
                d.preserveAR = opt.preserveAR;
                d.useTextFeatures = opt.useTextFeatures;
                d.filterTextRegions = opt.filterTextRegions;

                d.Ncategories = length(d.catDirs);
                
                % Initialize
                d.initialize();
            end

        end
        
        % Load image data
        function loadData(data)
            fprintf('Loading data...\n');
            tic;
            
            data.cat = cell(data.Ncategories,1);
            data.NtotalImages = 0;
            maxN = 0;
            imDims = [];
            
            for i=1:data.Ncategories
                catdir = data.catDirs{i};
                [catImages, catTxt, catDims] = loadImages(strcat(data.baseDir,catdir), data.DIM, data.preserveAR, data.useTextFeatures, data.filterTextRegions);
                Nimages = size(catImages,1);
                data.cat{i}.i = i;
                data.cat{i}.name = catdir;
                data.cat{i}.images = catImages;
                data.cat{i}.Nimages = Nimages;
                data.cat{i}.txt = catTxt;
                
                data.NtotalImages = Nimages + data.NtotalImages;
                if (Nimages > maxN)
                    maxN = Nimages;
                end
                imDims = [imDims; catDims];
                fprintf('%i/%i Image Categories Loaded\n',i, data.Ncategories);
            end
            data.maxCatSize = maxN;
            meanDim = mean(imDims,1);
            loadtime = toc;
            fprintf('%i images loaded in %.2d s. Average Dimensions: %u %u\n', data.NtotalImages, loadtime, round(meanDim));
            data.loaded = true;
        end
        
        % Return category name corresponding to index i
        function c = catName(d, i)
            c = d.cat{i}.name;
        end
        
        % Convenience function to return train-test workspace
        function ws = getWS(d)
            ws.trainX = d.trainX;
            ws.testX = d.testX;
            ws.trainY = d.trainY;
            ws.testY = d.testY;
        end
        
        % Convenience function to return entire dataset in a single workspace
        function ws = getWSFull(d)
            ws.X = [];
            ws.Y = [];
            ws.categoryNames = cell(d.Ncategories,1);
            for i=1:d.Ncategories
                ws.X = [ws.X; vertcat(d.cat{i}.images)];
                ws.Y = [ws.Y; repmat(d.cat{i}.i, d.cat{i}.Nimages, 1)];
                ws.categoryNames(i) = {d.cat{i}.name};
            end
            ws.C = ws.categoryNames(ws.Y);
        end
        
        % Get TextFeatureSet descriptor vectors at subdiv level l
        function ws = getTextFeatureSetVectors(d, subdivLevel)
            ws.X = [];
            ws.Y = [];
            ws.categoryNames = cell(d.Ncategories,1);
            for i=1:d.Ncategories
                ws.X = [ws.X; d.cat{i}.txt(:).featureVector(d.DIM, subdivLevel)];
                ws.Y = [ws.Y; repmat(d.cat{i}.i, d.cat{i}.Nimages, 1)];
                ws.categoryNames(i) = {d.cat{i}.name};
            end
            ws.C = ws.categoryNames(ws.Y);
        end
        
        % Computes a mask into data conditioned on current Ntrain, Ntest,
        % randomize state
        function newMask(d)
            if (d.Ntrain + d.Ntest > d.maxCatSize)
                error('Ntrain+Ntest > maxCatSize');
            end
            d.mask = cell(d.Ncategories,1);
            
            for i=1:d.Ncategories
                if (d.randomizeImageOrder)
                    r = randperm(d.cat{i}.Nimages);
                else
                    r = 1:(d.cat{i}.Nimages);
                end
                d.mask{i}.trainInd = r(1:d.Ntrain);
                d.mask{i}.testInd = r((d.Ntrain+1):(d.Ntrain+d.Ntest));
            end
        end
        
        function initialize(d)
            if ~d.loaded
               d.loadData();
            end
            
%             d.newMask();
        end

        % Getter methods
        function n = get.Ncategories(d)
            n = length(d.catDirs);
        end
        
        function trainX = get.trainX(d)
            trainX = [];
            for i=1:d.Ncategories
                trainX = [trainX; d.cat{i}.images(d.mask{i}.trainInd,:)];
            end
        end
        
        function testX = get.testX(d)
            testX = [];
            for i=1:d.Ncategories
                testX = [testX; d.cat{i}.images(d.mask{i}.testInd,:)];
            end
        end
        
        function trainY = get.trainY(d)
            trainY = [];
            for i=1:d.Ncategories
                trainY = [trainY; repmat(d.cat{i}.i,d.Ntrain,1)];
            end
            if (d.binaryClassification > 0)
               m = (trainY == d.binaryClassification);
               trainY(m) = 1;
               trainY(~m) = 2;
            end
        end
        
        function testY = get.testY(d)
            testY = [];
            for i=1:d.Ncategories
                testY = [testY; repmat(d.cat{i}.i,d.Ntest,1)];
            end
            if (d.binaryClassification > 0)
               m = (testY == d.binaryClassification);
               testY(m) = 1;
               testY(~m) = 2;
            end
        end
        
        % Setter methods
        function set.Ntrain(d, n)
           d.Ntrain = n;
        end
        
        function set.Ntest(d, n)
           d.Ntest = n;
        end
        
        function set.binaryClassification(d, catIndex)
            if (catIndex <= 0) % Set multi-class classification
                fprintf('Multi-class classification for %i classes.\n', d.Ncategories);
                d.binaryClassification = 0;
            else % Set binary classification on given category index
                if (catIndex > length(d.cat))
                    error('Binary Classification Category Index exceeds number of Categories in DataSet');
                else
                    fprintf('Binary classification on category: %i : %s', catIndex, d.catDirs{catIndex});
                    d.binaryClassification = catIndex;
                end
            end
        end
        
        function set.randomizeImageOrder(d, bool)
            d.randomizeImageOrder = bool;
        end
        
        % Disp method
        function disp(d)
            fprintf(1,'DataSet with %d Categories, Ntrain=%i, Ntest=%i\n',...
                d.Ncategories, d.Ntrain, d.Ntest);
        end % disp

    end % methods
    
end % classdef
