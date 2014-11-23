%% RANSAC for circles
%
%
function [best_params, best_inliers, best_error, varargout] ...
    = circleRANSAC(edgelIm, numIterations, initialInliers, minInliers, numBestEllipses)

    DISPLAY_ON = false;
    INTERACTIVE = false;
    
    NO_ELLIPSE = 100;
    
    if ~exist('numIterations', 'var')
        numIterations = 5000;
    end

    if ~exist('initialInliers', 'var')
        initialInliers = 3;
    end
    
    if ~exist('minInliers', 'var')
        minInliers = 20;
    end
    
    if ~exist('numBestEllipses','var')
        numBestEllipses = 10;
    end
    
    % Array keeping track of the best ellipses currently found. Has the
    % form:
    %
    % [params error]
    %
    % This is a 6 x numBestEllipses array.
    best_ellipses = zeros(numBestEllipses,6);
    best_ellipses(:,6) = NO_ELLIPSE;
    
    minFitDistance = 5;
    iterations = 1;

    nonzeroEdgels = find(edgelIm ~= 0);
    [iInds jInds] = ind2sub(size(edgelIm), nonzeroEdgels);
    XY = [jInds iInds];

    best_params = zeros(1,5);
    best_inliers = [];
    best_initial_points = [];
    best_error = Inf;
    best_inlier_error = Inf;
    best_proportion_error = Inf;

    if DISPLAY_ON
        debugIm = zeros(size(edgelIm));
        cmap = [1 1 1;
                0 0 1;
                0 1 0;
                1 0 0];
    end

    total_time = 0;
    ellipse_fit_time = 0;
    ellipse_error_time = 0;
    ellipse_error_count = 0;
    ellipse_residual_time = 0;
    while iterations <= numIterations && best_error > 0.8
        
        tStart = tic;
        
        if mod(iterations, 100) == 0 || iterations == 1
            fprintf(2, 'Current iteration: %d\n', iterations);
            if mod(iterations, 100) == 0
                fprintf(2, 'Time for last 100 iterations: %6.4f seconds\n', total_time);
                fprintf(2, 'ellipse_residual_time: %6.4f | ellipse_fit_time: %6.4f | ellipse_error_time: %6.4f\n', ...
                           ellipse_residual_time, ellipse_fit_time, ellipse_error_time);
                total_time = 0;
                ellipse_error_time = 0;
                ellipse_error_count = 0;
                ellipse_residual_time = 0;
                ellipse_fit_time = 0;
            end
        end

        try
        % Choose the initial set of inliers
        dataPermutation = XY(randperm(length(XY)), :);
        new_inliers = dataPermutation(1:initialInliers, :);
        new_outliers = dataPermutation(initialInliers+1:end,:);

        if INTERACTIVE
            figure; imshow(edgelIm); hold on;
            for k=1:initialInliers
                tr = ginput(1)+0.5; %% Shift pixel borders to integer positions
                plot(tr(1), tr(2), '--ro', 'MarkerEdgeColor', 'b', 'MarkerFaceColor','b','MarkerSize',10);
                new_inliers(k,:) = tr;
            end            
            close;
        end
        
        tFit = tic;
        if (max(new_inliers(:,1)) - min(new_inliers(:,1)) == 0) || ...
           (max(new_inliers(:,2)) - min(new_inliers(:,2)) == 0)
            iterations = iterations + 1;
            continue;
        end
        
        % Fit the circle
        [xc,yc,r] = circfit(new_inliers(:,1),new_inliers(:,2));
        ellipseParams = [xc,yc,r,r,0];
        
        % Break if the radius is larger than half the width or height
        if r >= min(size(edgelIm))
            iterations = iterations + 1;
            continue;
        end
        
        tRes = tic;
        err = Residuals_circle(new_outliers, ellipseParams);
        ellipse_residual_time = ellipse_residual_time + toc(tRes);
        
        % Find additional points that fit the model
        additional_inlier_inds = find(err <= minFitDistance);
        additional_inliers = new_outliers(additional_inlier_inds, :);

        consensus_set = [new_inliers; additional_inliers];
        
        if(length(consensus_set) >= minInliers)
            %new_ellipse = fitellipse(consensus_set(:,1), consensus_set(:,2));
            %[dist_error new_proj] = Residuals_ellipse(consensus_set, new_ellipse);
            
            [xc,yc,r] = circfit(consensus_set(:,1), consensus_set(:,2));
            new_ellipse = [xc,yc,r,r,0];
            dist_error = Residuals_circle(consensus_set, new_ellipse);
            new_proj = Projections_circle(consensus_set, new_ellipse);
            
            % Penalize inliers that lie outside the ellipse.
            radmat = repmat([new_ellipse(1) new_ellipse(2)], [size(consensus_set,1) 1]);
            proj_dist_to_radius = sqrt(sum((new_proj-radmat).^2,2));
            set_dist_to_radius = sqrt(sum((consensus_set-radmat).^2,2));
            outside_ellipse = set_dist_to_radius > proj_dist_to_radius;
            dist_error(outside_ellipse) = dist_error(outside_ellipse)*1.2;
            dist_error = sum(dist_error)/size(consensus_set,1);
            
            tEllipseError = tic;
            new_proj = ceil(new_proj);
            in_ind = find(new_proj(:,1) >= 1 & new_proj(:,1) <= size(edgelIm,2) & ...
                          new_proj(:,2) >= 1 & new_proj(:,2) <= size(edgelIm,1));
            new_proj = new_proj(in_ind, :);
            proj_indices = sub2ind(size(edgelIm), new_proj(:,2), new_proj(:,1));
            
            % How much of the ellipse is covered by the edgel points
            good_ellipse_error = 1-min(1,length(unique(proj_indices))/(2*pi*new_ellipse(3)));
            ellipse_error_time = ellipse_error_time + toc(tEllipseError);
            ellipse_error_count = ellipse_error_count + 1;
            
            alpha = 0.5;
            beta = 3;
            gamma = 2;
            
            good_aspect_error = 1-min(new_ellipse(3),new_ellipse(4))/max(new_ellipse(3), new_ellipse(4));
            
            %total_error = dist_error*alpha + proportion_error*beta;
            %total_error = dist_error*alpha + good_ellipse_error*beta;
            total_error = dist_error*alpha + good_ellipse_error*beta + good_aspect_error*gamma;
            %total_error = dist_error;
            
            if(total_error < best_error)
                best_params = new_ellipse;
                best_inliers = consensus_set;
                best_initial_points = new_inliers;
                best_error = total_error;
                
                % Add the new best ellipse to the list.
                best_ellipses(2:size(best_ellipses,1),:) = best_ellipses(1:size(best_ellipses,1)-1,:);
                best_ellipses(1,:) = [best_params best_error];
                
                if DISPLAY_ON
                    fprintf(2, '\nBest so far! Iteration %d\n', iterations);
                end               
                
                fprintf(2,'# points in proj im: %d | ellipse circ: %6.4f | dist_error: %6.4f |\ngood_ellipse_error: %6.4f | good_aspect_error: %6.4f | total_error: %6.4f\n',...
                        length(unique(proj_indices)), ellipse_circumference(new_ellipse), dist_error*alpha, good_ellipse_error*beta, good_aspect_error*gamma, total_error);                
            end

            if DISPLAY_ON
                % Create inlier/outlier image
                debugIm = zeros(size(edgelIm));
                figure; imshow(edgelIm); colormap(cmap); hold on;
                %scatter(new_outliers(:,1), new_outliers(:,2), 'MarkerEdgeColor','r','MarkerFaceColor', 'red', 'SizeData', 1);
                %scatter(additional_inliers(:,1), additional_inliers(:,2), 'MarkerEdgeColor', 'g','MarkerFaceColor', 'g', 'SizeData', 100);
                %scatter(new_inliers(:,1), new_inliers(:,2), 'MarkerEdgeColor','b','MarkerFaceColor', 'b', 'SizeData', 100);
                scatter(new_proj(:,1), new_proj(:,2), 'MarkerEdgeColor', 'b', 'MarkerFaceColor','b','SizeData',40);
                drawEllipse(ellipseParams, 'blue');
                drawEllipse(new_ellipse, 'green');
                hold off;

                fprintf(2,'# points in proj im: %d | ellipse circ: %6.4f | dist_error: %6.4f |\ngood_ellipse_error: %6.4f | total_error: %6.4f\n',...
                        length(unique(proj_indices)), ellipse_circumference(new_ellipse), dist_error*alpha, good_ellipse_error*beta, total_error);

                fprintf(2,'Press any key to continue...\n'); pause;
                close all;  
            end
        end
        
        catch ME1
            fprintf(2,'Something went wrong! %s', ME1.message);
        end
        
        iterations = iterations + 1;
        
        total_time = total_time + toc(tStart);
    end
    
    best_ellipses = best_ellipses(best_ellipses(:,6) < NO_ELLIPSE,:);
    
    varargout{1} = best_initial_points;
    varargout{2} = best_ellipses;
end