%% RANSAC for rectangles
%
%  A rectangle consists of four line segments.
%
function [best_params, best_inliers, best_error] = rect_RANSAC(edgelIm, initialInliers, minInliers, numIterations)

    DISPLAY_ON = false;
    INTERACTIVE = false;

    if ~exist('initialInliers', 'var')
        initialInliers = 4;
    end
    
    if ~exist('minInliers', 'var')
        minInliers = 20;
    end
    
    if ~exist('numIterations', 'var')
        numIterations = 200;
    end
    
    minFitDistance = 5;

    nonzeroEdgels = find(edgelIm ~= 0);
    [iInds jInds] = ind2sub(size(edgelIm), nonzeroEdgels);
    XY = [jInds iInds];

    best_params = struct();
    best_inliers = [];
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

    iterations = 1;
    while iterations <= numIterations

        fprintf(2, 'Current iteration: %d\n', iterations);
        
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
        
        % Create the rectangle
        minX = min(new_inliers(:,1));
        maxX = max(new_inliers(:,1));
        minY = min(new_inliers(:,2));
        maxY = max(new_inliers(:,2));
        
        rect = [minX minY; ...
                minX maxY; ...
                maxX maxY; ...
                maxX minY; ...
                minX minY];
            
%         if DISPLAY_ON
%             figure; imshow(edgelIm); hold on;
%             plot(new_inliers(:,1), new_inliers(:,2), '--ro', 'MarkerSize', 10);
%             plot(rect(:,1), rect(:,2), 'Color', 'red', 'LineWidth', 2);
%             hold off;
%         end
        
        % Calculate the errors to each line segment
        line1 = [rect(1, :); rect(2,:) - rect(1,:)];
        [err1 XYproj1] = Residuals_line(new_outliers, line1(1,:), line1(2,:), true, true);
        line2 = [rect(2,:); rect(3,:)-rect(2,:)];
        [err2 XYproj2] = Residuals_line(new_outliers, line2(1,:), line2(2,:), true, true);
        line3 = [rect(3,:); rect(4,:)-rect(3,:)];
        [err3 XYproj3] = Residuals_line(new_outliers, line3(1,:), line3(2,:), true, true);
        line4 = [rect(4,:); rect(5,:)-rect(4,:)];        
        [err4 XYproj4] = Residuals_line(new_outliers, line4(1,:), line4(2,:), true, true);
        
        %err = [err1 err2];
        %XYproj = [XYproj1 XYproj2];
        err = [err1 err2 err3 err4];
        XYproj = [XYproj1 XYproj2 XYproj3 XYproj4];
        
        clear err1 err2 err3 err4 XYproj1 XYproj2 XYproj3 XYproj4
        
        [minerr minerrind] = min(err, [], 2);
        final_err = zeros(length(err), 1);
        final_proj = zeros(length(err), 2);
        for i=1:length(minerrind)
            final_err(i, :) = err(i, minerrind(i));
            final_proj(i, :) = XYproj(i, 2*(minerrind(i)-1)+1:2*(minerrind(i)-1)+2);
        end        
        
        % Find additional points that fit the model
        additional_inlier_inds = find(final_err <= minFitDistance);
        additional_inliers = new_outliers(additional_inlier_inds, :);

        consensus_set = [new_inliers; additional_inliers];
        
%            if DISPLAY_ON
%                 % Create inlier/outlier image
%                 debugIm = zeros(size(edgelIm));
%                 figure; imshow(edgelIm); colormap(cmap); hold on;
%                 %scatter(new_outliers(:,1), new_outliers(:,2), 'MarkerEdgeColor','r','MarkerFaceColor', 'red', 'SizeData', 1);
%                 %scatter(additional_inliers(:,1), additional_inliers(:,2), 'MarkerEdgeColor', 'g','MarkerFaceColor', 'g', 'SizeData', 100);
%                 %scatter(new_inliers(:,1), new_inliers(:,2), 'MarkerEdgeColor','b','MarkerFaceColor', 'b', 'SizeData', 100);
%                 scatter(final_proj(:,1), final_proj(:,2), 'MarkerEdgeColor', 'b', 'MarkerFaceColor','b','SizeData',30);
%            
%                 fprintf(2,'Press any key to continue...\n'); pause;
%                 close all;
%            end        
        
        if(length(consensus_set) >= minInliers)
            minX = min(consensus_set(:,1));
            maxX = max(consensus_set(:,1));
            minY = min(consensus_set(:,2));
            maxY = max(consensus_set(:,2));

            new_rect = [minX minY; ...
                    minX maxY; ...
                    maxX maxY; ...
                    maxX minY; ...
                    minX minY];            
            
            line1 = [rect(1, :); rect(2,:) - rect(1,:)];
            [err1 XYproj1] = Residuals_line(consensus_set, line1(1,:), line1(2,:), true, true);
            line2 = [rect(2,:); rect(3,:)-rect(2,:)];
            [err2 XYproj2] = Residuals_line(consensus_set, line2(1,:), line2(2,:), true, true);
            line3 = [rect(3,:); rect(4,:)-rect(3,:)];
            [err3 XYproj3] = Residuals_line(consensus_set, line3(1,:), line3(2,:), true, true);
            line4 = [rect(4,:); rect(5,:)-rect(4,:)];        
            [err4 XYproj4] = Residuals_line(consensus_set, line4(1,:), line4(2,:), true, true);

            %err = [err1 err2];
            %XYproj = [XYproj1 XYproj2];
            err = [err1 err2 err3 err4];
            XYproj = [XYproj1 XYproj2 XYproj3 XYproj4];

            clear err1 err2 err3 err4 XYproj1 XYproj2 XYproj3 XYproj4

            [minerr minerrind] = min(err, [], 2);
            final_err = zeros(length(err), 1);
            final_proj = zeros(length(err), 2);
            for i=1:length(minerrind)
                final_err(i, :) = err(i, minerrind(i));
                final_proj(i, :) = XYproj(i, 2*(minerrind(i)-1)+1:2*(minerrind(i)-1)+2);
            end 
            
            % Distance error
            dist_error = mean(final_err);
            
            % How much of the rectangle is covered
            % # unique pixels
            unique_pixels = length(unique(final_proj, 'rows'));
            % Rectangle perimeter
            perimeter = 2*(maxX-minX) + 2*(maxY-minY);
            good_rectangle_error = 1-double(unique_pixels)/double(perimeter);
            
            alpha = 0.5;
            beta = 5;
            
            %out_of_bounds_error = 1000;
            
            %total_error = dist_error*alpha + proportion_error*beta;
            total_error = dist_error*alpha + good_rectangle_error*beta;
            %total_error = dist_error;           
            
            if(total_error < best_error)
                best_params.rect = new_rect;
                
                best_inliers = consensus_set;
                best_error = total_error;
                

                fprintf(2,'dist_error: %6.4f |\ngood_rectangle_error: %6.4f | total_error: %6.4f\n',...
                        dist_error*alpha, good_rectangle_error*beta, total_error);                
                
                if DISPLAY_ON
                    fprintf(2, '\nBest so far!\n');
                end
                
           if DISPLAY_ON
                % Create inlier/outlier image
                debugIm = zeros(size(edgelIm));
                figure; imshow(edgelIm); colormap(cmap); hold on;
                %scatter(new_outliers(:,1), new_outliers(:,2), 'MarkerEdgeColor','r','MarkerFaceColor', 'red', 'SizeData', 1);
                %scatter(additional_inliers(:,1), additional_inliers(:,2), 'MarkerEdgeColor', 'g','MarkerFaceColor', 'g', 'SizeData', 100);
                %scatter(new_inliers(:,1), new_inliers(:,2), 'MarkerEdgeColor','b','MarkerFaceColor', 'b', 'SizeData', 100);
                scatter(final_proj(:,1), final_proj(:,2), 'MarkerEdgeColor', 'b', 'MarkerFaceColor','b','SizeData',30);
                
                hold off;

                fprintf(2,'dist_error: %6.4f |\ngood_rectangle_error: %6.4f | total_error: %6.4f\n',...
                        dist_error*alpha, good_rectangle_error*beta, total_error);
           
                fprintf(2,'Press any key to continue...\n'); pause;
                close all;
           end                
            end
        end
        
        close all;
        iterations = iterations + 1;
    end
end