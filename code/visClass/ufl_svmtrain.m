function [svm, val, labels] = ufl_svmtrain(X, Y)

svm.C = 1000;
svm.testAccuracy = [];
svm.theta = train_svm(X, Y, svm.C);

[val,labels] = max(X*svm.theta, [], 2);
svm.trainAccuracy = 100 * (1 - sum(labels ~= Y) / length(Y));
fprintf('Train accuracy %f%%\n', svm.trainAccuracy);

end