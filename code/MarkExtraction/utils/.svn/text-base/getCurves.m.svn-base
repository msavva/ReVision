% Extract arcs from an edge image.
%
% INPUTS
%
% edgeIm - a binary image with edge points
% dirIm - an image containing values from [0, 2*Pi) with normal directions
% magIm - an image containing gradient magnitudes
function [curves, toosmall, dirIm_restricted] = getCurves(edgeIm, dirIm, magIm, maxDirDiff)

    DISPLAY_TIME = true;

    if ~exist('maxDirDiff','var')
        maxDirDiff = pi/16;
    end
    minCurveLength = 5;

    %dirIm(dirIm >= pi) = 2*pi - dirIm(dirIm >= pi);
    
    % Save the gradient direction values for pixels in the edgel image.
    % Set all non-values to NaN.
    dirIm_restricted = dirIm;
    dirIm_restricted(edgeIm == 0) = 0;
    dirIm_restricted(dirIm_restricted == 0) = NaN;
    dirIm_restricted = padarray(dirIm_restricted, [1 1], NaN);
    
    magIm_restricted = magIm;
    magIm_restricted(edgeIm == 0) = 0;
    magIm_restricted = padarray(magIm_restricted, [1 1], NaN);
    
    paddedImSize = size(dirIm_restricted);
    imSize = size(dirIm);
    
    curves = {};
    toosmall = {};
   
    total_time = 0;
    curPix = 1;
    while ~isempty(find(~isnan(dirIm_restricted)))
        tic
        
        if DISPLAY_TIME && (mod(curPix,1000) == 0 || curPix == 1)
            fprintf(2, 'Current pixel: %d\n', curPix);
            if mod(curPix, 1000) == 0
                fprintf(2, 'Time for last 1000 pixels: %6.4f seconds\n', total_time);
                total_time = 0;
            end        
        end
        
        % Start with a random non-zero point
        nonzeroPoints = find(~isnan(dirIm_restricted));
        curInd = nonzeroPoints(ceil(length(nonzeroPoints*rand(1,1))));
        %curInd = find(magIm_restricted == max(max(magIm_restricted)));

        curvePoints = stack();
        toConsider = stack(curInd);

        while ~isempty(toConsider)
            curInd = toConsider.pop();

            %sprintf('Value of dirIm_restricted(%d): %d',curInd, dirIm_restricted(curInd))
            [csi, csj] = ind2sub(paddedImSize, curInd);
            %sprintf('Value of dirIm_restricted(%d, %d): %d', csi, csj, dirIm_restricted(csi,csj))

            % Find differences in direction with the curve's neighbors

            % Extract the subimage consisting of a pixel and its immediate
            % neighbors
            neighborDirIm = dirIm_restricted((csi-1):min(size(dirIm_restricted,1),csi+1), ...
                                             (csj-1):min(size(dirIm_restricted,2),csj+1));

            % Compute the absolute directional differences
            neighborDirIm(~isnan(neighborDirIm)) = abs(neighborDirIm(~isnan(neighborDirIm)) - neighborDirIm(2,2));

%             neighborMagIm = magIm_restricted((csi-1):min(size(magIm_restricted,1),csi+1), ...
%                                              (csj-1):min(size(magIm_restricted,2),csj+1));
%             
%             neighborMagIm(~isnan(neighborMagIm)) = abs(neighborMagIm(~isnan(neighborMagIm)) - neighborMagIm(2,2));
                                         
            % Find nearby edge points by thresholding the directional
            % differences
            
            % Add the current point to the current curve point set.
            if ~isnan(dirIm_restricted(curInd))
                [cii, cij] = ind2sub(paddedImSize, curInd);
                curvePoints.push(sub2ind(imSize, cii-1, cij-1));
                dirIm_restricted(curInd) = NaN;
                magIm_restricted(curInd) = NaN;
            end
            
            % Add additional points to consider
            if(~isempty(find(~isnan(neighborDirIm))))
                %sprintf('adding new points')
%                 nextPoints = find(neighborDirIm <= maxDirDiff & ~isnan(neighborDirIm) & ...
%                                   neighborMagIm <= 5);
                nextPoints = find(neighborDirIm <= maxDirDiff & ~isnan(neighborDirIm));
                for k=1:length(nextPoints)
                    [smdi, smdj] = ind2sub(size(neighborDirIm),nextPoints(k));
                    toConsider.push(sub2ind(paddedImSize, csi+(smdi-2), csj+(smdj-2)));
                end
            end
        end

        % Only keep a curve if it's longer than a certain number of points
        if(length(curvePoints) > minCurveLength)
            curves{end+1} = curvePoints.data;
        else
            toosmall{end+1} = curvePoints.data;
        end
        
        total_time = total_time + toc;
        curPix = curPix + 1;
    end
end