function [labels, testX] = binaryClassifier(testX, theta, centroids, rfSize, CIFAR_DIM, M, P, trainXC_mean, trainXC_sd)
%% Configuration
addpath minFunc;

fprintf('Extracting features...\n');
tic;
% compute testing features and standardize
testXC = extract_features(testX, centroids, rfSize, CIFAR_DIM, M,P);
testXCs = bsxfun(@rdivide, bsxfun(@minus, testXC, trainXC_mean), trainXC_sd);
testXCs = [testXCs, ones(size(testXCs,1),1)];
toc

fprintf('Classifying...\n');
[val,labels] = max(testXCs*theta, [], 2);

end
