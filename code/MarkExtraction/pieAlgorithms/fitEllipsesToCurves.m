function [allEllipses boundingBoxArea] = fitEllipsesToCurves(curves, imsize)
    % Fit an ellipse to each curve
    allEllipses = {};
    boundingBoxArea = [];
    for k=1:length(curves)
        x = [];
        y = [];
        for z=1:length(curves{k})
            [yind xind] = ind2sub(imsize, curves{k}(z));
            x(end+1) = xind;
            y(end+1) = yind;
        end
        boundingBoxArea(end+1) = (max(x)-min(x)+1)*(max(y)-min(y)+1);
        % Don't consider vertical or horizontal lines
        if( (max(x)-min(x) ~= 0) && (max(y)-min(y) ~= 0))
            ellipseParams = fitellipse(x',y');
            allEllipses{end+1} = ellipseParams;
        else
            allEllipses{end+1} = repmat(NaN, [1 5]);
        end
    end
end