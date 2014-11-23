classdef barFinder < handle

    properties (Constant)
        % Bar orientation flag
        VERTICAL_BARS = 1;
        HORIZONTAL_BARS = 2;
    end    
    
    properties
        % Initialization: chart image, cache, and current file
        im
        cache
        curFile_noext
        
        % Connected components properties
        components
        componentStats
        cc_time
        ccstats_time

        % Color check metadata
        colorcheck_time
        
        % Stores indices of all the rectangles.
        rectangles = [];

        % Stores indices of candidate bars.
        candidate_bars = [];
        candidate_widths = [];
        candidate_heights = [];

        % Stores the indices of bars that failed the color check or were too small
        other_bars = [];
        other_widths = [];
        other_heights = [];    

        % Stores the indices of rectangles that were too small, but passed the
        % color check.
        small_bars = [];
        small_widths = [];
        small_heights = [];
        
        % Stores the final bars
        final_bars = [];
        
        % Axis properties
        axis1 = struct('x1',NaN,'y1',NaN','x2',NaN,'y2',NaN);
        axis2 = struct('x1',NaN,'y1',NaN','x2',NaN,'y2',NaN);
        axis_peak
        
        % Orientation properties
        candidate_width
        candidate_height
        orientation
    end
    
    methods
    
        %% Constructor
        function obj = barFinder(im, curFile, cache)
            obj.im = im;
            obj.cache = cache;
            
            extInd = strfind(curFile, '.');
            obj.curFile_noext = curFile(1:extInd-1);
        end
        
        %% Computation methods
        % Find color-connected components and store statistics about those
        % components.
        function findConnectedComponents(obj)
        
            % Check for cached connected components
            cc_file = strcat(obj.curFile_noext,'_colorcc_nobg.mat');
            loadPath = obj.cache.getLoadPath(cc_file);

            % Retrieve the cached image
            if ~isempty(loadPath)
                fprintf(2,'Loading cached components...\n');
                load(loadPath);
            else
            % Generate the connected components if they aren't cached
                cctic = tic;
                % Connected components method
                fprintf(2,'Finding color-connected components...\n');
                [components, n] = conncompcolor(obj.im);
                cc_time = toc(cctic);
                save(obj.cache.getSavePath(cc_file, obj.cache.SHARED_CACHE), ...
                     'components','cc_time');
            end
            
            obj.components = components;
            obj.cc_time = cc_time;
        

            % Checked for cached cc stats
            ccstat_file = strcat(obj.curFile_noext,'_ccstats_',sprintf('%d', obj.cache.BFILTER_WIDTH), '_', ...
                 sprintf('%2.2f',obj.cache.BFILTER_DOMAINSIGMA),'_',sprintf('%2.2f',obj.cache.BFILTER_RANGESIGMA),'.mat');
            loadPath = obj.cache.getLoadPath(ccstat_file);

            % Retrieve cached files
            if ~isempty(loadPath)
                fprintf(2,'Loading cached cc stats...\n');
                load(loadPath);
                
            % Generate the image if it's not cached
            else
                % Compute bounding boxes and centroids for each connected component
                fprintf(2,'Finding connected component statistics...\n');
                ccstattic = tic;
                componentStats = regionprops(components, 'all');
                ccstats_time = toc(ccstattic);
                save(obj.cache.getSavePath(ccstat_file, obj.cache.LOCAL_CACHE), ...
                    'componentStats','ccstats_time');
            end
            
            obj.componentStats = componentStats;
            obj.ccstats_time = ccstats_time;
        end
        
        % Find the rectangles in the image. Classifies components as
        % rectangles if their area covers at least rectThreshold of their
        % bounding box. Stores rectangles that are below an area threshold
        % as small bars.
        function findRectangles(obj, rectThreshold, areaThreshold)
            
            % Hack to speed this function up
            cs = obj.componentStats;
            
            for zz=1:length(cs)
                % Check if the component is likely a rectangle by comparing the
                % bounding box area to the number of filled pixels. Discard any
                % rectangles with an area less than a threshold (25).
                curComp = cs(zz);
                boundingBoxArea = size(curComp.FilledImage, 1)* ...
                                  size(curComp.FilledImage, 2);
                onPoints = length(curComp.PixelIdxList);

                if (onPoints/boundingBoxArea) > rectThreshold

                    if boundingBoxArea <= areaThreshold

                        w = curComp.BoundingBox(3);
                        h = curComp.BoundingBox(4);                

                        obj.small_bars(end+1) = zz;
                        obj.small_widths(end+1) = round(w);
                        obj.small_heights(end+1) = round(h);
                    else
                        obj.rectangles(end+1) = zz;
                    end

                end
            end
        end

        % V2 COLOR CHECK
        % For the midpoints of each rectangle side
        % 1. Check color of pixel immediately adjacent, perpendicular to the
        % side. Call this color C.
        % 2. Walk perpendicularly out from the side until either
        %   (a) A pixel has a color significantly different from C (call this
        %   color D), or
        %   (b) We've walked out more than 5 pixels.
        function doColorCheck(obj, COLOR_CHECK_STEP, COLOR_THRESH)
            fprintf(2, 'Performing color check...\n');
            colorcheck_t = tic;
            
            local_cs = obj.componentStats;
            
            imr = obj.im(:,:,1);
            img = obj.im(:,:,2);
            imb = obj.im(:,:,3);     
            
            for zz=1:length(obj.rectangles)
                if mod(zz, 100) == 0
                    fprintf(2, 'Color check iteration %d\n', zz);
                end

                %%
                curComp = local_cs(obj.rectangles(zz));
                curColor = avgColor(obj.im(curComp.SubarrayIdx{:},:));
                % Ensure color ranges from 0 to 1
                if(max(curColor) > 1)
                    curColor = curColor/255;
                end

                ul = curComp.BoundingBox(1:2);
                w = curComp.BoundingBox(3);
                h = curComp.BoundingBox(4);

                adjColorsLocs = [ul(2)+(h/2) ul(1)-1 ;
                             ul(2)+(h/2) ul(1)+w+1 ;
                             ul(2)-1 ul(1)+(w/2) ;
                             ul(2)+h+1 ul(1)+(w/2)];
                %adjColorsLocs = round(adjColorsLocs);
                adjColorsLocs = ceil(adjColorsLocs);

                % Keep values in-bounds
                acl_col_1 = adjColorsLocs(:,1);
                acl_col_2 = adjColorsLocs(:,2);
                acl_col_1(acl_col_1 < 1) = 1;
                acl_col_1(acl_col_1 > size(imr,1)) = size(imr,1);
                acl_col_2(acl_col_2 < 1) = 1;
                acl_col_2(acl_col_2 > size(imr,2)) = size(imr,2);
                adjColorsLocs = [acl_col_1 acl_col_2];

                adjColors_inds = sub2ind(size(imr),adjColorsLocs(:,1),adjColorsLocs(:,2));
                adjColors = [imr(adjColors_inds) img(adjColors_inds) imb(adjColors_inds)];

                %%
                % Color check on each side
                passedColorCheck = true;
                for yy=1:size(adjColorsLocs,1)
                    
                    initialLocation = adjColorsLocs(yy,:);
                    initialColor = adjColors(yy,:);
                    testColor = [];

                    for qq=0:COLOR_CHECK_STEP-1

                        getColor = true;

                        testLocation = initialLocation;
                        % Left edge midpoint
                        if yy==1
                            % Move an additional unit left
                            testLocation(1,2) = testLocation(1,2)-qq;
                            % Out-of-bounds check
                            if testLocation(1,2) < 1
                                getColor = false;
                            end
                        % Right edge midpoint
                        elseif yy==2
                            % Move an additional unit right
                            testLocation(1,2) = testLocation(1,2)+qq;
                            % Out-of-bounds check
                            if testLocation(1,2) > size(imr,2);
                                getColor = false;
                            end
                        % Top edge midpoint
                        elseif yy==3
                            % Move an additional unit up
                            testLocation(1,1) = testLocation(1,1)-qq;
                            % Out-of-bounds check
                            if testLocation(1,1) < 1
                                getColor = false;
                            end
                        % Bottom edge midpoint
                        else
                            % Move an additional unit down
                            testLocation(1,1) = testLocation(1,1)+qq;
                            % Out-of-bounds check
                            if testLocation(1,1) > size(imr,1);
                                getColor = false;
                            end                    
                        end

                        if getColor
                            % Retrieve the color at the new location
                            testLocation_ind = sub2ind(size(imr), testLocation(1), testLocation(2));
                            testLocation_color = [imr(testLocation_ind) img(testLocation_ind) imb(testLocation_ind)];

                            % Compute MSE
                            colorMSE = sqrt(sum((testLocation_color-initialColor).^2,2));

                            % If the color is different than the immediately adjacent color, or if the maximum number
                            % of steps has been taken outside the rectangle, set
                            % the test color.
                            if colorMSE > COLOR_THRESH || qq == COLOR_CHECK_STEP-1
                                testColor = testLocation_color;
                                break;
                            end
                        end

                    end

                    if getColor
                        % Compare the test color with the average component color
                        colorCheckMSE = sqrt(sum((curColor-testColor).^2,2));
                        if colorCheckMSE < COLOR_THRESH
                            passedColorCheck = false;
                            break;
                        end
                    end

                end

                % Classify component as a candidate bar or non-candidate bar.
                if passedColorCheck
                    obj.candidate_bars(end+1) = obj.rectangles(zz);
                    obj.candidate_widths(end+1) = round(w);
                    obj.candidate_heights(end+1) = round(h);
                else
                    obj.other_bars(end+1) = obj.rectangles(zz);
                    obj.other_widths(end+1) = round(w);
                    obj.other_heights(end+1) = round(h);
                end
            end
            obj.colorcheck_time = toc(colorcheck_t);        
        end
        
        function findCandidateBars(obj, minDim, MINIMUM_CANDIDATE_BARS)
            obj.removeSmallBars(minDim);
        
            % If no candidate bars are found (e.g., all rectangles failed the color
            % check because they are the same color as the background), classify
            % rectangles that are large enough as candidate bars.
            if length(obj.candidate_bars) <= MINIMUM_CANDIDATE_BARS
                bigenough = obj.other_widths > 2 & obj.other_heights > 2;

        %         small_bars = [];
        %         small_widths = [];
        %         small_heights = [];
                obj.candidate_bars = obj.other_bars(bigenough);
                obj.candidate_widths = obj.other_widths(bigenough);
                obj.candidate_heights = obj.other_heights(bigenough);

                obj.other_bars = obj.other_bars(~bigenough);
                obj.other_widths = obj.other_widths(~bigenough);
                obj.other_heights = obj.other_heights(~bigenough); 
            end        
        end
        
        % Remove bars who have a dimension smaller than minDim
        % from the candidate bars collection.
        function removeSmallBars(obj, minDim)
            bigenough = obj.candidate_widths > minDim & obj.candidate_heights > minDim;

            obj.other_bars = cat(2,obj.other_bars,obj.candidate_bars(~bigenough));
            obj.other_widths = cat(2,obj.other_widths,obj.candidate_widths(~bigenough));
            obj.other_heights = cat(2,obj.other_heights,obj.candidate_heights(~bigenough));
        
            obj.small_bars = cat(2,obj.small_bars,obj.candidate_bars(~bigenough));
            obj.small_widths = cat(2,obj.small_widths,obj.candidate_widths(~bigenough));
            obj.small_heights = cat(2,obj.small_heights,obj.candidate_heights(~bigenough));
            obj.candidate_bars = obj.candidate_bars(bigenough);
            obj.candidate_widths = obj.candidate_widths(bigenough);
            obj.candidate_heights = obj.candidate_heights(bigenough);
        end
        
        % Determine orientation of the bar chart. Requires
        % candidate_heights and candidate_widths to be filled.
        function value = inferOrientation(obj)
            edges = 1:max(obj.candidate_widths);
            nwidths = histc(obj.candidate_widths, edges);
            obj.candidate_width = max(edges(nwidths == max(nwidths)));

            % height histogram
            edges = 1:max(obj.candidate_heights);
            nheights = histc(obj.candidate_heights, edges);
            obj.candidate_height = max(edges(nheights == max(nheights)));

            % Infer bar orientation
            obj.orientation = obj.VERTICAL_BARS;
            if length(unique(obj.candidate_widths)) > length(unique(obj.candidate_heights))
                obj.orientation = obj.HORIZONTAL_BARS;
            end
            
            value = obj.orientation;
        end
        
        function findBaselineAxis(obj, DIM_SLACK_PERC, AXIS_SLACK, COLOR_THRESH)
            
            local_cs = obj.componentStats;
            
            % Get the colors of the candidate bars
            candidate_bar_colors = zeros(length(obj.candidate_bars),3);
            for zz = 1:length(obj.candidate_bars)
                curComp = local_cs(obj.candidate_bars(zz));
                candidate_bar_colors(zz,:) = avgColor(obj.im(curComp.SubarrayIdx{:},:));
            end            
            
            if obj.orientation == obj.VERTICAL_BARS
                candidate_dims = obj.candidate_widths;
                small_dims = obj.small_widths;
                bar_dim = obj.candidate_width;
                baseline_ind = 2;
                perp_ind = 1;
            else
                candidate_dims = obj.candidate_heights;
                small_dims = obj.small_heights;
                bar_dim = obj.candidate_height;
                baseline_ind = 1;
                perp_ind = 2;
            end
            
            filtered_candidate_bars = [];
            
            DIM_SLACK = max(bar_dim*DIM_SLACK_PERC,5);
            
            % Keep bars whose width is near the bar width.
            dim_diff = abs(candidate_dims-mean(bar_dim));
            filtered_candidate_bars = obj.candidate_bars(dim_diff < DIM_SLACK);
            good_bars = cat(2,dim_diff < DIM_SLACK,zeros(1,length(obj.small_bars)));

            dim_vals = zeros(1,(length(filtered_candidate_bars)+1)*2);
            for zz=1:length(filtered_candidate_bars)
                curComp = local_cs(filtered_candidate_bars(zz));
                curComp_extrema = unique(curComp.Extrema, 'rows');
                dim_vals(zz*2+1:zz*2+2) = [min(curComp_extrema(:,baseline_ind)) max(curComp_extrema(:,baseline_ind))];
            end
            dim_vals = round(dim_vals);               

            edges = 1:max(dim_vals);
            ndims= histc(dim_vals, edges);
            obj.axis_peak = edges(ndims == max(ndims));

            % Vertical bar: Take leftmost peak in case of finding multiple peaks
            % Horizontal bar: Take bottommost peak in case of finding
            % multiple peaks
            if obj.orientation == obj.VERTICAL_BARS
                obj.axis_peak = obj.axis_peak(1);
                
                % Baseline axis
                obj.axis1.y1 = obj.axis_peak;
                obj.axis1.y2 = obj.axis_peak;   
                
                % Perpendicular axis
                obj.axis2.y1 = obj.axis_peak;
            else
                obj.axis_peak = obj.axis_peak(end);
                
                % Baseline axis
                obj.axis1.x1 = obj.axis_peak;
                obj.axis1.x2 = obj.axis_peak;
                
                % Perpendicular axis
                obj.axis2.x1 = obj.axis_peak;
            end

            %%
            % Only keep bars that are on the axis, are close to the candidate width, 
            % and that passed the color check

            % Stores the x- or y-extents of the candidate bars. We assume no stacked
            % charts, so we don't add any bars that are above or below an
            % existing candidate bar.
            %
            % [ (min extent) (max extent) (0 if candidate bar, distance to axis if small bar); ... ]
            candidate_bar_extents = zeros(length(obj.candidate_bars),3);

            for zz=1:length(obj.candidate_bars)
                curComp = local_cs(obj.candidate_bars(zz));
                curComp_extrema = round(unique(curComp.Extrema, 'rows'));

                candidate_bar_extents(zz,:) = [min(curComp_extrema(:,perp_ind)) max(curComp_extrema(:,perp_ind)) 0];

                close_to_axis = abs(min(curComp_extrema(:,baseline_ind))-obj.axis_peak) <= AXIS_SLACK || ...
                   abs(max(curComp_extrema(:,baseline_ind))-obj.axis_peak) <= AXIS_SLACK;
                good_dim = abs(candidate_dims(zz)-bar_dim) <= DIM_SLACK;

                if close_to_axis && good_dim
                    good_bars(zz) = good_bars(zz) & 1;
                else
                    good_bars(zz) = 0;
                end
            end

            % Add bars that were too small but have the same width as existing
            % bars AND passed the color check
            for zz=1:length(obj.small_bars)
                curComp = local_cs(obj.small_bars(zz));
                curComp_extrema = round(unique(curComp.Extrema, 'rows'));
                dist_to_axis = min(abs(min(curComp_extrema(:,baseline_ind))-obj.axis_peak), ...
                                   abs(max(curComp_extrema(:,baseline_ind))-obj.axis_peak));
                close_to_axis = dist_to_axis <= AXIS_SLACK;
                good_dim = abs(small_dims(zz)-bar_dim) <= DIM_SLACK;

                %%
                add_bar = false;
                if close_to_axis && good_dim

                    %%
                    % Check if vertical space that small bar occupies is already
                    % occupied by a candidate bar or another small bar that's
                    % closer to the axis.
                    cur_bar_extents = [min(curComp_extrema(:,perp_ind)) max(curComp_extrema(:,perp_ind)) dist_to_axis];
                    %%
                    overlapping_bar = false;
                    for yy=1:size(candidate_bar_extents,1)
                        curext = candidate_bar_extents(yy,:);
                        % There is already a bar in the space
                        if (cur_bar_extents(1) >= curext(1) && cur_bar_extents(1) <= curext(2)) || ...
                           (cur_bar_extents(2) >= curext(1) && cur_bar_extents(2) <= curext(2)) || ...
                           (cur_bar_extents(1) <= curext(1) && cur_bar_extents(2) >= curext(2))

                            % The bar is either an existing candidate bar or is a
                            % small bar that is closer to the axis
                            if curext(3) == 0 || curext(3) < dist_to_axis
                                fprintf(2, '%d\n', yy);
                                overlapping_bar = true;
                            end
                        end
                    end
                    if ~overlapping_bar
                        % Only add the bar if its color matches one of the
                        % existing bars
                        curColor = avgColor(obj.im(curComp.SubarrayIdx{:},:));
                        if max(curColor) > 1
                            curColor = curColor / 255;
                        end
                        curColorTest = repmat(curColor, [size(candidate_bar_colors,1) 1]);
                        testcolordiffs = sqrt(sum((curColorTest-candidate_bar_colors).^2,2));
                        if ~isempty(find(testcolordiffs < COLOR_THRESH))
                            candidate_bar_extents(end+1,:) = cur_bar_extents;
                            good_bars(zz+length(obj.candidate_bars)) = 1;
                            add_bar = true;                        
                        end
                    end
                end

                if ~add_bar
                    good_bars(zz+length(obj.candidate_bars)) = 0;
                end
            end

            all_bars = cat(2, obj.candidate_bars, obj.small_bars);
            obj.final_bars = all_bars(good_bars == 1);
        end
        
        %% Utility methods
        
        % Return a label image and colormap of the rectangles contained in
        % rectIndList.
        %
        % Args
        %   rectIndList A vector of rectangle indices in obj.componentStats
        %   colorMode One of 'average' or 'center'. Used for creating the
        %       output color map -- use either the average color of the
        %       component, or the center pixel value. Defaults to
        %       'average'.
        function [outIm colormap] = getImage(obj, rectIndList, colorMode)
            
            if ~exist('colorMode', 'var')
               colorMode = 'average';
            end
            local_cs = obj.componentStats;
            
            outIm = zeros(size(obj.im,1),size(obj.im,2));
            
            if (nargout == 2)
                colormap = zeros(length(rectIndList),3);
            end
            
            for i=1:length(rectIndList)
                curComp = local_cs(rectIndList(i));
                outIm(curComp.PixelIdxList) = i;
            
                % Output a color map along with the image
                if (nargout == 2)
                    % Use the average color of the component
                    if strcmp('average', colorMode)
                        curColor = avgColor(obj.im(curComp.SubarrayIdx{:},:));
                    % Use the color of the center pixel of the component
                    else
                        curIm = obj.im(curComp.SubarrayIdx{:}, :);
                        center = [max(1, floor(size(curIm, 1) / 2))
                                  max(1, floor(size(curIm, 2) / 2))];
                        curColor = [curIm(center(1), center(2), 1)
                                    curIm(center(1), center(2), 2)
                                    curIm(center(1), center(2), 3)];
                    end
                    % Ensure color ranges from 0 to 1
                    if(max(curColor) > 1)
                        curColor = curColor/255;
                    end
                    colormap(i,:) = curColor;
                end
            end
        end
        
        %% Getters and setters
        function value = get.axis_peak(obj)
            value = obj.axis_peak;
        end
        
        function value = get.orientation(obj)
            value = obj.orientation;
        end
        
        function value = get.componentStats(obj)
           value = obj.componentStats;
        end
    end
end