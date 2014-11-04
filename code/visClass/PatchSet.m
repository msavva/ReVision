classdef PatchSet < handle
    
    properties
        N = 50000;              % Number of patches
        Ncentroids = 150;       % Number of centroids
        rfSize = 6;             % Receptor Field Size (i.e. Patch Size)
        whitening = true;       % Whether to use whitening
        normContrast = true;    % Whether to normalize patches for contrast
        DIM = [128 128 1]       % Image Dimensions
        kmeansIterations = 100  % Iterations for kmeans centroid computation
        patches = [];           % [N X rfSize^2] Patch Data Matrix
        patchLabels = [];       % [N X 1] Vector assigning category indices to each patch
        M = [];                 % Patch Mean Matrix
        P = [];                 % Patch Alignment Matrix (Right-Multiplies for whitening)
        centroids = [];         % Patch Centroids (computed through k-means clustering)
        centroidFrequency = []; % Centroid Occurence Frequencies
        
        MIN_PATCH_VAR = 38;    % Minimum Patch Variance for accepting as potential centroid (empirically set to about 25% quartile of var)
    end % Instance properties
    
    methods(Static)
        % K-means using VL_feat routine
        function [C, Cfreq] = kmeansVL(X, k)
            X = X';
            [C, A] = vl_kmeans(X, k, 'algorithm', 'elkan', 'initialization', 'plusplus');
            % Also sort centroids in descending order of popularity
            count = hist(single(A),single(unique(A)));
            [Cfreq, I] = sort(count,'descend');
            Cfreq = 100 .* (Cfreq ./ size(X,2));
            C = C(:,I)';
        end
        
        function [A, M, P] = normalizeAndWhiten(A)
            % normalize for contrast
            pM = mean(A,2);
            pSqVar = sqrt(var(A,[],2)+10);
            A = bsxfun(@rdivide, bsxfun(@minus, A, pM), pSqVar);
            
            % whiten
            C = cov(A);
            M = mean(A);
            [V,D] = eig(C);
            P = V * diag(sqrt(1./(diag(D) + 0.1))) * V';
            A = bsxfun(@minus, A, M) * P;
        end
        
        function A = invertWhiteningAndNormalization(A, M, P)
            A = bsxfun(@plus, (A / P), M);
%             A = bsxfun(@plus, bsxfun(@times, A, mean(pSQVAR)), mean(pM));
        end

        function Q = subdivPooling(X,l)
            n = min(size(X));
            split = round(n/2);
            if (l==0)
                Q = squeeze(sum(sum(X,1),2));
            else
                Q = [PatchSet.subdivPooling(X(1:split     , 1:split    , :), l-1);
                     PatchSet.subdivPooling(X(split+1:end , 1:split    , :), l-1);
                     PatchSet.subdivPooling(X(1:split     , split+1:end, :), l-1);
                     PatchSet.subdivPooling(X(split+1:end , split+1:end, :), l-1)];
            end
        end
    end
    
    methods
        % Constructor
        function p = PatchSet(opt, data, labels)
            if nargin > 0               % Pull variables
                p.N = opt.Npatches;
                p.Ncentroids = opt.Ncentroids;
                p.kmeansIterations = opt.kmeansIterations;
                p.rfSize = opt.rfSize;
                p.whitening = opt.whitening;
                p.DIM = opt.DIM;
            end
            % Extract patches from data
            if nargin == 2               
                p.extractPatches(data);
            elseif nargin == 3
                p.extractPatches(data, labels);
            end
        end
        
        % Extract Random Patches from Data (if labels provided, save patch
        % label as well)
        function extractPatches(p, data, labels)
            fprintf('Extracting random patches from data...\n');
            tic;
            
            A = zeros(p.N, p.rfSize*p.rfSize*p.DIM(3));
            L = ones(p.N, 1);
            nx = p.DIM(1); ny = p.DIM(2); nc = p.DIM(3);
            rf = p.rfSize;
            i=1;
            trials=0;
            while (i <= p.N)
                if (mod(trials,10000) == 0)
                    fprintf('%d / %d patches accepted.\n', i, p.N);
                end
                
                r = random('unid', nx - rf + 1);
                c = random('unid', ny - rf + 1);
                index = mod(i-1,size(data,1))+1;
                patch = reshape(data(index, :), [nx ny nc]);
                patch = patch( r:r+rf-1 , c:c+rf-1 , : );
                if (var(patch(:)) > p.MIN_PATCH_VAR)
                    A(i,:) = patch(:)';
                    L(i) = labels(index);
                    i = i + 1;
                end
                trials = trials + 1;
            end
            
            p.patches = A;
            p.patchLabels = L;
            
            fprintf('%i patches extracted in %d\n', p.N, toc);

        end
        
        % Use current centroids to extract a feature set for the given data
        % matrix (subdivLevels defines # of block subdivisions at which to sum feature
        % vectors)
        function fs = extractFeatures(p, X, subdivLevels, trainedFS)
            fprintf('Extracting feature vectors using centroid PatchSet...\n');
            tic;
            
            % If no subdivLevels given assume single quartering
            if nargin < 3
                subdivLevels = 1;
            end
            
            % Get local copies of variables
            rf = p.rfSize;
            pDIM = p.DIM;
            pM = p.M;
            pP = p.P;
            pWhitening = p.whitening;
            pCentroids = p.centroids;
            cc = sum(pCentroids.^2, 2)';
            pNcentroids = p.Ncentroids;
            
            sz = pDIM(1)*pDIM(2);
            
            XC = zeros(size(X,1), (4^subdivLevels)*p.Ncentroids);

            parfor i=1:size(X,1)
                if (mod(i,1000) == 0)
                    fprintf('Extracting features: %d / %d\n', i, size(X,1));
                end
                
                % extract overlapping sub-patches into rows of 'patches'
                ps = im2col(reshape(X(i,1:sz),pDIM(1:2)), [rf rf])';
                
                % do preprocessing for each patch
                
                % normalize for contrast
                ps = bsxfun(@rdivide, bsxfun(@minus, ps, mean(ps,2)), sqrt(var(ps,[],2)+1));

                % whiten
                if (pWhitening)
                    ps = bsxfun(@minus, ps, pM) * pP;
                end
                
                % compute 'triangle' activation function
                xx = sum(ps.^2, 2);
                xc = ps * pCentroids';
                
                z = sqrt( bsxfun(@plus, cc, bsxfun(@minus, xx, 2*xc)) ); % distances
                [v,inds] = min(z,[],2);
                mu = mean(z, 2); % average distance to centroids for each patch
                ps = max(bsxfun(@minus, mu, z), 0);
                
                %-- NOTE: 1-of-K hard assignment overwrites activations for now --%
                off = 0:pNcentroids:(size(z,1)-1)*pNcentroids;
                ps = ps(:);
                ps(:) = 0;
                ps(off'+inds) = 1;
                ps = reshape(ps, size(z,2), size(z,1))';
                
                
                % reshape to numCentroids-channel image
                prows = pDIM(1)-rf+1;
                pcols = pDIM(2)-rf+1;
                ps = reshape(ps, prows, pcols, pNcentroids);

                % pool over quadrants and concatenate into feature vector
                XC(i,:) = PatchSet.subdivPooling(ps,subdivLevels)';
            end
            
            fprintf('%i feature vectors computed in %d\n', size(X,1), toc);
            
            % Save standardized features into FeatureSet struct
            fs.XC = XC;
            if (nargin > 3) % trained FeatureSet with mean and sd passed in so use that to normalize this FeatureSet
                fs.mean = trainedFS.mean;
                fs.sd = trainedFS.sd;
            else
                fs.mean = mean(XC);
                fs.sd = sqrt(var(XC)+0.01);
            end
            fs.XCs = bsxfun(@rdivide, bsxfun(@minus, XC, fs.mean), fs.sd);
            fs.XCs = [fs.XCs, ones(size(fs.XCs,1),1)];
        end
        
        % K-means for Centroid Computation (uses VL_feat subroutine)
        function kmeansCentroids(ps)
            [normedPatches, ps.M, ps.P] = PatchSet.normalizeAndWhiten(ps.patches);
            [ps.centroids, ps.centroidFrequency] = PatchSet.kmeansVL(normedPatches, ps.Ncentroids);
        end
        
        % Hierarchical K-means (two-tier using Ncentroids per category and
        % reducing to Ncentroids total)
        function kmeansCentroids2Tier(ps)
            nCategories = max(ps.patchLabels);
            C = [];
            Cfreq = [];
            for i=1:nCategories
                origPatches = ps.patches(ps.patchLabels == i , :);
                [normedPatches, Mi, Pi] = PatchSet.normalizeAndWhiten(origPatches);
                [Ci, Cfreqi] = PatchSet.kmeansVL(normedPatches, ps.Ncentroids);
                ps.M = Mi;
                ps.P = Pi;
                ps.centroids = Ci;
                ps.centroidFrequency = Cfreqi;
                ps.showCentroids;
                Ci = PatchSet.invertWhiteningAndNormalization(Ci, Mi, Pi);
                C = [C; Ci];
                Cfreq = [Cfreq; Cfreqi];
            end
            [normedC, ps.M, ps.P] = PatchSet.normalizeAndWhiten(C);
            [ps.centroids, ps.centroidFrequency] = PatchSet.kmeansVL(normedC, ps.Ncentroids);
        end
        
        % K-means for Centroid Computation
        function kmeansCentroidsOLD(ps, NITER)
            iterations = ps.kmeansIterations;
            if nargin > 1
                iterations = NITER;
            end
            
            X = ps.patches;
            k = ps.Ncentroids;
            x2 = sum(X.^2,2);
            C = randn(k,size(X,2))*0.1;%X(randsample(size(X,1), k), :);
            BATCH_SIZE=1000;
            
            
            for itr = 1:iterations
                fprintf('K-means iteration %d / %d\n', itr, iterations);
                
                c2 = 0.5*sum(C.^2,2);
                
                summation = zeros(k, size(X,2));
                counts = zeros(k, 1);
                
                loss =0;
                
                for i=1:BATCH_SIZE:size(X,1)
                    lastIndex=min(i+BATCH_SIZE-1, size(X,1));
                    m = lastIndex - i + 1;
                    
                    [val,labels] = max(bsxfun(@minus,C*X(i:lastIndex,:)',c2));
                    loss = loss + sum(x2(i:lastIndex) - val');
                    
                    S = sparse(1:m,labels,1,m,k,m); % labels as indicator matrix
                    summation = summation + S'*X(i:lastIndex,:);
                    counts = counts + sum(S,1)';
                end
                
                
                C = bsxfun(@rdivide, summation, counts);
                
                % just zap empty centroids so they don't introduce NaNs everywhere.
%                 badIndex = find(counts == 0);
%                 C(badIndex, :) = 0;
                C(counts == 0, :) = 0;
            end
            
            ps.centroids = C;
        end
        
        % Show a Patch within the PatchSet
        function showPatch(p, i)
            figure;
            imshow(reshape(p.patches(i,:) , p.rfSize , p.rfSize , p.DIM(3)));
        end
        
        % Show Centroids for this PatchSet
        function showCentroids(p, highlight)
            if nargin < 2
                highlight = [];
            end
            % Check for object array case and empty case
            nP = length(p);
            if nP>1
                for i=1:nP
                    showCentroids(p(i));
                end
                return;
            elseif isempty(p.centroids)
                error('Centroids have not been computed.');
            end
            
            highlight = mod(highlight, p.Ncentroids);
            figure;
            
            H = p.rfSize;
            W = H;
            
            No=size(p.centroids,2)/(H*W);
            assert(No == 3 || No == 1);  % color and gray images
            
            K=size(p.centroids,1);
            COLS=round(sqrt(K));
            ROWS=ceil(K / COLS);
            
            % Invert whitening and normalization
            C = PatchSet.invertWhiteningAndNormalization(p.centroids, p.M, p.P);
            C = (C .* 40) + 190; % Approximate contrast denormalization for visibility (empirical values for mean and sqvar)
            C(C < 0) = 0;
            C(C > 255) = 255;
            
            clf; hold on;
            image=ones(ROWS*(H+1), COLS*(W+1), No)*100;
            for i=1:size(p.centroids,1)
                r= floor((i-1) / COLS);
                c= mod(i-1, COLS);
                centr = reshape(C(i,1:W*H*No),H,W,No);
                if (any(highlight == i))
                    centr(1:2,1:2) = [0 255; 255 0];
                end
                image((r*(H+1)+1):((r+1)*(H+1))-1,(c*(W+1)+1):((c+1)*(W+1))-1,:) = centr;
            end
            
%             mn=-1.5;
%             mx=+1.5;
%             image = (image - mn) / (mx - mn);

            imshow(image, [0 255]);
            title('Centroid Patches');
            
            % Also plot centroid frequencies as histogram
            figure;
            bar(p.centroidFrequency);
            title('Centroid Occurence Frequencies');
            drawnow;
        end
        
        % Disp method
        function disp(p)
            n = length(p);
            if n>1
                for i=1:n
                    disp(p(i));
                end
                return;
            else
                fprintf('PatchSet with %d Patches of size %i^2\n', p.N, p.rfSize);
            end
        end % disp

    end % methods
    
end % classdef
