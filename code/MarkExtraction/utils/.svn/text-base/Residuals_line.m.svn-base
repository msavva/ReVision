%%
% Compute distance of the points XY to a line.
% Endpoints are (x0, x0+d). This function will compute the residuals for
% lines, semi-infinite lines, or line segments. fixPoint1 fixes x0, and
% fixPoint2 fixes x0+d.
%
% XY contains points in the form [x0 y0; x1 y1; ...]
%
% XY_proj contains the projections of the input points onto the line
% or line segment.
%
%
function [RSS, XY_proj, varargout] = Residuals_line(XY, x0, d, fixPoint1, fixPoint2)

    if ~exist('fixPoint1', 'var')
        fixPoint1 = true;
    end
    
    if ~exist('fixPoint2', 'var')
        fixPoint2 = true;
    end

    % Compute the normal
    n = [d(2) -d(1)];
    n = n/sqrt(n(1)^2+n(2)^2);    

    % Find the projection of the outliers onto the line, not considering the
    % endpoints.
    XY_proj = repmat(x0, [size(XY,1) 1])-XY;
    dotProd = XY_proj;
    dotProd(:,1) = dotProd(:,1)*n(1);
    dotProd(:,2) = dotProd(:,2)*n(2);
    dotProd = sum(dotProd,2);
    XY_proj = XY+repmat(n, [size(XY,1) 1]).*[dotProd dotProd];

    % Check for overrun of the endpoints
    aMat = XY_proj-repmat(x0, [size(XY,1) 1]);
    if d(2) == 0
        aMat = aMat*[1/d(1); 1];
    else
        aMat = aMat*[1; (1-d(1))/d(2)];
    end
    
    % Assign the closest endpoint as the projection
    lt = find(aMat < 0);
    gt = find(aMat > 1);
    if fixPoint1
        XY_proj(lt, :) = repmat(x0, [length(lt) 1]);
        aMat(lt, :) = repmat(0, [length(lt) 1]);
    end
    if fixPoint2
        XY_proj(gt, :) = repmat(x0+d, [length(gt) 1]);
        aMat(gt, :) = repmat(1, [length(gt) 1]);
    end
    clear lt gt dotProd

    RSS = sqrt(sum((XY-XY_proj).^2,2));

    varargout{1} = aMat;
end