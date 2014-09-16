function [val, labels, acc] = ufl_svmpredict(svm, X, Y)

[val,labels] = max(X*svm.theta, [], 2);
acc = 100 * (1 - sum(labels ~= Y) / length(Y));
svm.testAccuracy = [svm.testAccuracy ; acc];
fprintf('Prediction accuracy %f%%\n', acc);

end