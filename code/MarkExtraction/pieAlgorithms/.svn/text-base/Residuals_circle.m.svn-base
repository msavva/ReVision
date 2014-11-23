function RSS = Residuals_circle(XY,params)
%
% Computing distance of a set of points to a circle.
%
% Input: XY(n,2) is the array of coordinates of n points x(i)=XY(i,1), y(i)=XY(i,2)
%        params is the vector of the circle parameters
%          params = [x_c, y_c, r]
%

    % Compute residual error
    XY(:,1) = (XY(:,1)-repmat(params(1),[size(XY,1) 1])).^2;
    XY(:,2) = (XY(:,2)-repmat(params(2),[size(XY,1) 1])).^2;
    RSS = abs(sqrt(XY(:,1)+XY(:,2))-params(3));

end   % Residuals_circle
