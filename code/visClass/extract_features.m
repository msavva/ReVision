function XC = extract_features(X, centroids, rfSize, CIFAR_DIM, M,P)
  assert(nargin == 4 || nargin == 6);
  whitening = (nargin == 6);
  numCentroids = size(centroids,1);
  
  % compute features for all training images
  XC = zeros(size(X,1), numCentroids*4);

  parfor i=1:size(X,1)
    if (mod(i,1000) == 0)
        fprintf('Extracting features: %d / %d\n', i, size(X,1));
    end
    
    % extract overlapping sub-patches into rows of 'patches'
    sz = CIFAR_DIM(1)*CIFAR_DIM(2);
    %--- Grayscale Hack ---%
    patches = im2col(reshape(X(i,1:sz),CIFAR_DIM(1:2)), [rfSize rfSize])';

    % do preprocessing for each patch
    
    % normalize for contrast
    patches = bsxfun(@rdivide, bsxfun(@minus, patches, mean(patches,2)), sqrt(var(patches,[],2)+10));
    % whiten
    if (whitening)
      patches = bsxfun(@minus, patches, M) * P;
    end
    
    % compute 'triangle' activation function
    xx = sum(patches.^2, 2);
    cc = sum(centroids.^2, 2)';
    xc = patches * centroids';
    
    z = sqrt( bsxfun(@plus, cc, bsxfun(@minus, xx, 2*xc)) ); % distances
    [v,inds] = min(z,[],2);
    mu = mean(z, 2); % average distance to centroids for each patch
    patches = max(bsxfun(@minus, mu, z), 0);
    % patches is now the data matrix of activations for each patch
    
    % reshape to numCentroids-channel image
    prows = CIFAR_DIM(1)-rfSize+1;
    pcols = CIFAR_DIM(2)-rfSize+1;
    patches = reshape(patches, prows, pcols, numCentroids);
    
    % pool over quadrants
    halfr = round(prows/2);
    halfc = round(pcols/2);
    q1 = sum(sum(patches(1:halfr, 1:halfc, :), 1),2);
    q2 = sum(sum(patches(halfr+1:end, 1:halfc, :), 1),2);
    q3 = sum(sum(patches(1:halfr, halfc+1:end, :), 1),2);
    q4 = sum(sum(patches(halfr+1:end, halfc+1:end, :), 1),2);
    
    % concatenate into feature vector
    XC(i,:) = [q1(:);q2(:);q3(:);q4(:)]';
  end

end