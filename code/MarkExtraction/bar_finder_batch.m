% Prototype to find bars in a bar chart. We explore a few methods to do
% this.

%% Bars via connected components
clear all
close all

setupPath();
thresholds();

IMAGE_READ_DIRECTORY = 1;
IMAGE_READ_SUBSET = 2;

%%%%%%%%%%%%%%%%%%%%%%
% Flags
IMAGE_READ_MODE = IMAGE_READ_DIRECTORY;

% [chart_files fullpath outpath] = getImagesInPath('Testing/UMDbars');
% textpath =
% '/Users/nick/Documents/Projects/VisionCharts/VisInterp/corpus/UMDImages/bar';
% cachepath = fullfile(pwd, 'cache');
%

%% Use a text file that specifies the files
if IMAGE_READ_MODE == IMAGE_READ_SUBSET

    chart_subset_file = 'Testing/largebars/large_bars.txt';
    %chart_subset_file = 'Testing/UMDbars/umd_bars.txt';
    chart_path = 'Visualization/BarChart';
    %chart_path = 'UMDImages/bar';
    chart_files = getchartfiles(chart_subset_file);
    outpath = fullfile(pwd, 'Testing/largebars');
    %outpath = fullfile(pwd, 'Testing/UMDbars');

    for i=1:length(chart_files)
        chart_files{i} = strcat(chart_files{i}, '.png');
    end

elseif IMAGE_READ_MODE == IMAGE_READ_DIRECTORY

    %% Use a directory

    % Large corpus paths
    %chart_path = 'Scraped/bing_bar chart/UIST11_assumptions';
    
    %chart_path = 'bing_bar chart/UIST11_assumptions'; % Large corpus path
    chart_path = 'thesisims'; % Charts for the textref project
    
    %outpath = fullfile(pwd, 'Testing/bing_bar chart');
    outpath = fullfile(pwd, 'Testing/output/thesisims');
    
    fullpath = fullfile(corpus_root, chart_path);
    listing = dir(fullpath);
    chart_files = {};
    curInd = 1;
    for i=1:size(listing,1)
        if (~isempty(regexpi(listing(i).name, '.jpg')) || ...
            ~isempty(regexpi(listing(i).name, '.png')) || ...
            ~isempty(regexpi(listing(i).name, '.gif')))
            chart_files{curInd} = listing(i).name;
            curInd = curInd + 1;

        end
    end

end
%%

fullpath = fullfile(corpus_root, chart_path);

cacheManager = CacheManager(BFILTER_WIDTH, BFILTER_DOMAINSIGMA, BFILTER_RANGESIGMA);
cacheManager.LocalCachePath = fullfile(pwd, 'cache');
cacheManager.SharedCachePath = sharedCache;
% Set to:
% "corpusbar" if using Stanford corpus
% "largebar" if using large bar corpus
% "bingbar" if using the Bing bar corpus
% "custom" if using the custom corpus
% empty if using UMD.
cacheManager.cacheprefix = 'custom';


% fullpath = fullfile(pwd, chart_path);
% outpath = fullfile(fullpath,'output');
% 
% % Retrieve file names for test charts
% all_image_path = dir(fullpath);
% 
% chart_files = [];
% for i=1:length(all_image_path)
%     if (~isempty(strfind(all_image_path(i).name, 'jpg')) || ...
%         ~isempty(strfind(all_image_path(i).name, 'png')))
%         chart_files = [chart_files; {all_image_path(i).name}];
%     end
% end

outString = 'file,proj_decision,proj_time,cc_decision,cc_time\n';
 
%chart_files = {'bar_59.png', 'bar_89.png', 'bar_91.png'}; % Vertical bars with successful data
%extraction
% Failed the color check, or bars were too short
%chart_files = {'bar_31.png', 'bar_30.png', 'bar_33.png',...
%'bar_34.png','bar_35.png','bar_103.png','bar_104.png','bar_120.png','bar_76.png','bar_73.png',...
%'bar_129.png','bar_81.png','bar_84.png','bar_19.png','bar_74.png'};
%chart_files = {'bar_104.png'};

% Initialize simple second derivative
[gFilt_raw gxFilt_raw gxxFilt_raw] = genGaussians(0.3);
[gFilt gxFilt gxxFilt] = genGaussians(3);

dataparams = fopen(fullfile(outpath,'dataparams.txt'), 'a');

fprintf(dataparams, 'file,min_value,scale\n');
%%
for i=1:length(chart_files)
%for i=106
    try
        
    total_tic = tic;
        
    fprintf(2,'Processing %s...\n',chart_files{i});

    % Read in the image
    curFile = chart_files{i};
    %curFile = 'bar_105.png';
    extInd = strfind(curFile, '.');
    curFile_noext = curFile(1:extInd-1);
    textlocFile = fullfile(fullpath, 'text', strcat(curFile(1:extInd-1),'.txt'));
    
    clear improc
    improc = imageProcessor(curFile, fullpath, cacheManager);
    improc.readImages();
    improc.filterImage();
    
    im = improc.filteredIm;
    imgray = improc.grayIm;
    
    handle = imshowlarge(improc.im);
    saveImage(handle, curFile, '_origim.png', outpath);
    close;
    
    handle = imshowlarge(im,true);
    saveImage(handle, curFile, '_filteredim.png',outpath);
    close;
    
    outString = strcat(outString, ',', curFile);
    
    fprintf(2, 'Reading text locations...\n');
    % Retrieve the text locations, if available
    % Remove the text if text locations exist
    if exist(textlocFile, 'file')
        [textmask textlocs] = gettextmask(textlocFile,im);
        clear textLocFile
    else
        % No text file exists, so don't mask anything
        textmask = ones(size(im,1),size(im,2));
        textlocs = [];
    end

    % Output text mask debug images if text mask exists
    if ~isempty(textlocs)
        handle = figure('visible','off'); imshow(textmask);
        fprintf(2,'Saving text mask...\n');
        saveImage(handle, curFile, '_textmask.png', outpath);            

        handle = figure('visible','off');
        imshow(im);
        hold on;
        for rectInd=1:length(textlocs.boundingboxes)
            currect = [textlocs.boundingboxes{rectInd} ; textlocs.boundingboxes{rectInd}(1,:)];
            plot(currect(:,2), currect(:,1), 'LineWidth', 1, 'Color', 'red');
        end
        hold off;
        fprintf(2,'Saving text mask overlaid on image...\n');
        saveImage(handle, curFile, '_textmaskoverlay.png', outpath);       
    end    
    
    clear bf
    bf = barFinder(improc.filteredIm, curFile, cacheManager);
    bf.findConnectedComponents();
    
    components = bf.components;
    componentStats = bf.componentStats;
    
    falseIm = bf.getImage(1:length(bf.componentStats));    
    fprintf(2,'Saving the color connected components image...\n');
    handle = imshowlarge(label2rgb(falseIm, 'colorcube'),true);
    saveImage(handle, curFile, '_colorcc.png',outpath);
    close;    
    
    fprintf(2,'Finding rectangles...\n');
    
    rectThreshold = RECTANGLE_COVERAGE;
    areaThreshold = 0;
    bf.findRectangles(rectThreshold, areaThreshold);
    
    % Stores indices of all the rectangles.
    rectangles = bf.rectangles;
    
    % Stores indices of candidate bars.
    candidate_bars = bf.candidate_bars;
    candidate_widths = bf.candidate_widths;
    candidate_heights = bf.candidate_heights;
    
    % Stores the indices of bars that failed the color check or were too small
    other_bars = bf.other_bars;
    other_widths = bf.other_widths;
    other_heights = bf.other_heights; 
    
    % Stores the indices of rectangles that were too small, but passed the
    % color check.
    small_bars = bf.small_bars;
    small_widths = bf.small_widths;
    small_heights = bf.small_heights;

    all_bars = bf.rectangles;
    
    fprintf(2,'Took %6.6f seconds.\n',bf.ccstats_time);
    
    %%
    % Check colors outside of rectangles
    COLOR_THRESH = 0.1;
    bf.doColorCheck(COLOR_CHECK_STEP, COLOR_THRESH);
    [falseIm cmap] = bf.getImage(bf.rectangles, 'center');
 
    fprintf(2, 'Color check time: %6.4f\n', bf.colorcheck_time);
    
    fprintf(2,'Saving the rectangle image...\n');
    %handle = imshowlarge(colorIm,true);
    handle = imshowlarge(label2rgb(falseIm, cmap, 'k'), true);
    saveImage(handle, curFile, '_rectangles.png',outpath);
    close;
    
    fprintf(2,'Saving the false rectangle image...\n');
    handle = imshowlarge(label2rgb(falseIm),true);
    saveImage(handle, curFile, '_falserectangles.png',outpath);
    close;
    
    %%
    
    % Throw out bars that are too thin or short
    %
    % If no candidate bars are found (e.g., all rectangles failed the color
    % check because they are the same color as the background), classify
    % rectangles that are large enough as candidate bars.
    bf.findCandidateBars(2, MINIMUM_CANDIDATE_BARS);
    
    candidate_bars = bf.candidate_bars;
    candidate_widths = bf.candidate_widths;
    candidate_heights = bf.candidate_heights;
    
    other_bars = bf.other_bars;
    other_widths = bf.other_widths;
    other_heights = bf.other_heights;
    

    %% Debug images
    [falseIm cmap] = bf.getImage(bf.candidate_bars);
    fprintf(2,'Saving the candidate bar image...\n');    
    %handle = imshowlarge(colorIm,true);
    if(isempty(cmap))
        cmap = colormap('jet');
    end
    handle = imshowlarge(label2rgb(falseIm, cmap, 'k'), false);
    saveImage(handle, curFile, '_candidateBars.png',outpath);
    close;
    %%
    
    fprintf(2,'Saving the false candidate bar image...\n');
    handle = imshowlarge(label2rgb(falseIm),true);
    saveImage(handle, curFile, '_falsecandidatebars.png',outpath);
    close;     
    
    [falseIm cmap] = bf.getImage(bf.small_bars);
    fprintf(2,'Saving the small bars image...\n');
    if(isempty(cmap))
        cmap = colormap('jet');
    end
    handle = imshowlarge(label2rgb(falseIm, cmap, 'w'),true);
    saveImage(handle, curFile, '_smallbars.png',outpath);
    close;   
    
    [falseIm cmap] = bf.getImage(bf.other_bars);
    fprintf(2,'Saving the rejected rectangles image...\n');    
    if(isempty(cmap))
        cmap = colormap('jet');
    end
    handle = imshowlarge(label2rgb(falseIm, cmap, 'w'),true);
    saveImage(handle, curFile, '_rejected.png',outpath);
    close;
    
    %% Gather width and height votes for each candidate bar    
    bar_orientation = bf.inferOrientation();
    candidate_width = bf.candidate_width;
    candidate_height = bf.candidate_height;
    
%     % width histogram
%     edges = 1:max(candidate_widths);
%     nwidths = histc(candidate_widths, edges);
%     candidate_width = max(edges(nwidths == max(nwidths)));
%     
%     % height histogram
%     edges = 1:max(candidate_heights);
%     nheights = histc(candidate_heights, edges);
%     candidate_height = max(edges(nheights == max(nheights)));
%     
%     % Infer bar orientation
%     bar_orientation = VERTICAL_BARS;
%     if length(unique(candidate_widths)) > length(unique(candidate_heights))
%         bar_orientation = HORIZONTAL_BARS;
%     end
    
    %% Gather colors of the candidate bars
    fprintf(2, 'Getting colors from candidate bars...\n');
    candidate_bar_colors = zeros(length(candidate_bars),3);
    for zz = 1:length(candidate_bars)
        curComp = componentStats(candidate_bars(zz));
        candidate_bar_colors(zz,:) = avgColor(im(curComp.SubarrayIdx{:},:));
    end
    
    %%
%     figure;
%     subplot(2,1,1);
%     plot(nwidths, 'ro');
%     subplot(2,1,2);
%     plot(nheights,'ro');
    
%%    
    fprintf(2,'Finding the axes...\n');
    
    bf.findBaselineAxis(0.10, AXIS_SLACK, COLOR_THRESH);
    
    final_bars = bf.final_bars;
    if bf.orientation == bf.VERTICAL_BARS
        xaxis_peak = bf.axis_peak;
    else
        yaxis_peak = bf.axis_peak;
    end
    
    [falseIm cmap] = bf.getImage(final_bars);
    fprintf(2,'Saving the final bar image...\n');     
    handle = imshowlarge(label2rgb(falseIm, cmap, 'w'),true);
    saveImage(handle, curFile, '_finalbars.png',outpath);
    close;
    
    marks = {};
    for zz=1:length(final_bars)
        curComp = componentStats(final_bars(zz));
        curColor = avgColor(im(curComp.SubarrayIdx{:},:));        
        new_mark = {};
        new_mark.stats = curComp;
        new_mark.color = curColor;
        marks{end+1} = new_mark;
    end

    %%
    % Find the non-baseline axis (y-axis if vertical bars, x-axis if
    % horizontal)
    if bar_orientation == bf.VERTICAL_BARS
        min_x = Inf;
        max_x = 0;
        % Find the leftmost and rightmost bar sides
        for zz= 1:length(final_bars)
            curComp = componentStats(final_bars(zz));
            curComp_extrema = unique(curComp.Extrema, 'rows');   
            if(min_x > min(curComp_extrema(:,1)))
                min_x = min(curComp_extrema(:,1));
            end
            if(max_x < max(curComp_extrema(:,1)))
                max_x = max(curComp_extrema(:,1));
            end
        end

        % Y-difference
        imlab_y_diff = imdiff(im, 'y');
        ily_inds = find(imlab_y_diff < mean(imlab_y_diff(:)));
        ily_t = imlab_y_diff;
        ily_t(ily_inds) = 0;
        if ~isnan(textmask)
            ily_t = ily_t & textmask;
        end
        yproj = sum(ily_t,1);
        
        mstd =  mean(yproj)+std(yproj);
        [yproj_peak_locs ypl0 ypl1 ypl2] = findsmoothedpeaks(yproj, gFilt_raw, gxFilt_raw, gxxFilt_raw,-1);
        yproj_peak_locs = yproj_peak_locs(ypl0(yproj_peak_locs) > mstd);
        
        disttomin = yproj_peak_locs - repmat(min_x, size(yproj_peak_locs));
        disttomax = repmat(max_x, size(yproj_peak_locs)) - yproj_peak_locs;
        [m d] = max(disttomin(disttomin < -AXIS_SLACK));
        left_axis = yproj_peak_locs(disttomin < -AXIS_SLACK);
        left_axis = left_axis(d);
        [m d] = max(disttomax(disttomax < -AXIS_SLACK));
        right_axis = yproj_peak_locs(disttomax < -AXIS_SLACK);
        right_axis = right_axis(d);
        
        bf.axis1.x1 = left_axis;
        bf.axis2.x1 = left_axis;
        bf.axis2.x2 = left_axis;
        
        %%
        axis_labels = {};
        bar_labels = {};
        text_labels = {};
        
        % Get text labels and attempt to parse them
        if ~isempty(textlocs)

            % Find labels to the left of the left-axis
            for zz=1:length(textlocs.boundingboxes)
                curbb = textlocs.boundingboxes{zz};
                max_col = max(curbb(:,2));
                if max_col <= left_axis
                    newlab = {};
                    newlab.text = textlocs.text{zz};
                    newlab.boundingbox = textlocs.boundingboxes{zz};
                    newlab.bbcenter = [(max(curbb(:,1))+min(curbb(:,1)))/2 (min(curbb(:,2))+max(curbb(:,2)))/2];
                    newlab.top = curbb(1,1);
                    newlab.left = curbb(1,2);
                    newlab.right = curbb(2,2);
                    newlab.bottom = curbb(3,1);
                    axis_labels{end+1} = newlab;
                end
            end

            % Find labels below the x-axis
            for zz=1:length(textlocs.boundingboxes)
                curbb = textlocs.boundingboxes{zz};
                % Top of the bounding box is below the x-axis
                if curbb(1,1) > xaxis_peak
                    newlab = {};
                    newlab.text = textlocs.text{zz};
                    newlab.boundingbox = textlocs.boundingboxes{zz};
                    bar_labels{end+1} = newlab;
                end
            end            
            
            %successful_parse = true;
            for zz=1:length(axis_labels)

                curtext = axis_labels{zz}.text;
                new_label = parseValueLabel(curtext);

                if ~isempty(new_label.value)
                    new_label.boundingbox = axis_labels{zz}.boundingbox;
                    new_label.top = axis_labels{zz}.top;
                    new_label.bottom = axis_labels{zz}.bottom;
                    new_label.center = (new_label.top+new_label.bottom)/2;
                    text_labels{end+1} = new_label;
                end
            end
        end

        % Calculate scale and origin
        scale = 1;
        min_value = 0;

        % Assign values and labels to marks
        bbset = cell(1,length(bar_labels));
        for zz=1:length(bar_labels)
            bbset{zz} = bar_labels{zz}.boundingbox;
        end
        
        tlset = cell(1,length(text_labels));
        for zz=1:length(text_labels)
            tlset{zz} = text_labels{zz}.boundingbox;
        end
        %%
        if ~isempty(text_labels)
            origin_y = xaxis_peak;

            % Project labels onto the axes
%             xaxis_line.x0 = [xaxis_peak left_axis];
%             xaxis_line.d = [0 size(im,2)-left_axis];
%             
%             tl_projections = zeros(length(text_labels),2);
%             for zz=1:length(text_labels)
%                 curbb = text_labels{zz}.boundingbox;
%                 centerpt = [(min(curbb(:,1))+max(curbb(:,1)))/2 (min(curbb(:,2))+max(curbb(:,2)))/2];
%                 tl_line.x0 = centerpt;
%                 tl_projections(zz,:) = lineintpoint(xaxis_line);
%             end
            
            % Assume linear scale
            scalearr = zeros(sum(1:length(text_labels)-1),1);
            ind = 1;
            for t1=1:length(text_labels)
                for t2=t1+1:length(text_labels)
                    scalearr(ind) = (abs(text_labels{t1}.value-text_labels{t2}.value) / abs(text_labels{t1}.center-text_labels{t2}.center));
                    ind = ind + 1;
                end
            end
            scale = median(scalearr);
            if scale == 0
                scale = 1;
            % Estimate the min-value
            else
                min_value_found = false;
                origin = [xaxis_peak left_axis];
                if ~isempty(tlset) && length(origin) == 2
                    % Estimate the min-value using the closest label
                    [clp min_dist] = closestLabelFromPoint(origin, tlset);

                    closest_origin_label = text_labels{clp};
                    if closest_origin_label.value == 0
                        min_value = 0;
                        min_value_found = true;
                    elseif min_dist < 10
                        min_value = closest_origin_label.value;
                        min_value_found = true;
                    end
                end
                
                if ~min_value_found
                    min_vals = zeros(length(text_labels),1);
                    for t1=1:length(text_labels)
                        min_vals(t1) = text_labels{t1}.value-scale*(origin_y-text_labels{t1}.center);
                    end
                    min_value = median(min_vals);
                end
            end
        end
        fprintf(dataparams,'%s,%2.2f,%2.2f\n',curFile,min_value,scale);      

        % Fill out mark structure and output data table
        for m=1:length(marks)
            % TODO: Deal with x-axes that are not the baseline
%             bartop = unique(marks{m}.stats.Extrema, 'rows');
%             bartop = min(bartop(:,2));
%             marks{m}.value = min_value
            
            height = marks{m}.stats.BoundingBox(4);
            marks{m}.value = min_value+scale*height;

            if ~isempty(textlocs) && ~isempty(bbset)
                labind = closestLabel(marks{m}.stats.BoundingBox,bbset);
                %marks{m}.label = textset{labind};
                marks{m}.label = bar_labels{labind}.text;
            else
                marks{m}.label = num2str(m);
            end

            marks{m}.id = m;
        end
        
        %%
        handle = imshowlarge(im,true); hold on;
        plot([1 size(im,2)], [xaxis_peak xaxis_peak], '-r');
        %plot([min_x min_x], [1 size(im,1)],'-r');
        %plot([max_x max_x], [1 size(im,1)],'-r');
        if ~isempty(left_axis)
            plot([left_axis left_axis], [1 size(im,1)],'-b');
        end
        if ~isempty(right_axis)
            plot([right_axis right_axis], [1 size(im,1)],'-b');
        end
        for zz=1:length(text_labels)
            estvalue = min_value+scale*(xaxis_peak-text_labels{zz}.center);
            estheight = xaxis_peak - (text_labels{zz}.value - min_value)/scale;
            scatter(left_axis, estheight, 'ro','filled');
        end
        hold off;
        fprintf(2, 'Saving axis marked image...\n');
        saveImage(handle, curFile, '_axes.png',outpath);
        close;
        %%
        
%         figure;
%         subplot(3,1,1);
%         hold on;
%         plot(ypl0);
%         scatter(yproj_peak_locs, ypl0(yproj_peak_locs));
%         hold off;
%         subplot(3,1,2);
%         plot(ypl1);
%         subplot(3,1,3);
%         plot(ypl2);
    
    % Horizontal bars
    else
        min_y = Inf;
        max_y = 0;
        % Find the leftmost and rightmost bar sides
        for zz= 1:length(final_bars)
            curComp = componentStats(final_bars(zz));
            curComp_extrema = unique(curComp.Extrema, 'rows');   
            if(min_y > min(curComp_extrema(:,2)))
                min_y = min(curComp_extrema(:,2));
            end
            if(max_y < max(curComp_extrema(:,2)))
                max_y = max(curComp_extrema(:,2));
            end
        end
        
        % X-difference
        imlab_x_diff = imdiff(im, 'x');
        ilx_inds = find(imlab_x_diff < mean(imlab_x_diff(:)));
        ilx_t = imlab_x_diff;
        ilx_t(ilx_inds) = 0;
        if ~isnan(textmask)
            ilx_t = ilx_t & textmask;
        end
        xproj = sum(ilx_t,2);
        
        %ily_t(ily_t ~= 0) = 1   
        
        % LAB projections
%         figure;
%         hold on;
%         plot(xproj);
%         plot([1 length(xproj)], [mean(xproj) mean(xproj)], '-r');
%         plot([1 length(xproj)], [mean(xproj)+std(xproj) mean(xproj)+std(xproj)], '-r');
%         hold off;
        
        mstd =  mean(xproj)+std(xproj);
        [xproj_peak_locs xpl0 xpl1 xpl2] = findsmoothedpeaks(xproj, gFilt_raw, gxFilt_raw, gxxFilt_raw,-1);
        xproj_peak_locs = xproj_peak_locs(xpl0(xproj_peak_locs) > mstd);
        
        disttomin = xproj_peak_locs - repmat(min_y, size(xproj_peak_locs));
        disttomax = repmat(max_y, size(xproj_peak_locs)) - xproj_peak_locs;
        [m d] = max(disttomin(disttomin < -AXIS_SLACK));
        top_axis = xproj_peak_locs(disttomin < -AXIS_SLACK);
        top_axis = top_axis(d);
        [m d] = max(disttomax(disttomax < -AXIS_SLACK));
        bottom_axis = xproj_peak_locs(disttomax < -AXIS_SLACK);
        bottom_axis = bottom_axis(d);
        
       %%
        axis_labels = {};
        bar_labels = {};
        text_labels = {};
        
        % Get text labels and attempt to parse them
        if ~isempty(textlocs)
            
            % Find labels below the bottom axis
            for zz=1:length(textlocs.boundingboxes)
                curbb = textlocs.boundingboxes{zz};
                % Top of the bounding box is below the x-axis
                if curbb(1,1) > bottom_axis
                    newlab = {};
                    newlab.text = textlocs.text{zz};
                    newlab.boundingbox = textlocs.boundingboxes{zz};
                    newlab.bbcenter = [(max(curbb(:,1))+min(curbb(:,1)))/2 (min(curbb(:,2))+max(curbb(:,2)))/2];
                    newlab.top = curbb(1,1);
                    newlab.left = curbb(1,2);
                    newlab.right = curbb(2,2);
                    newlab.bottom = curbb(3,1);
                    axis_labels{end+1} = newlab;
                end
            end

            % Find labels to the left of the y axis
            for zz=1:length(textlocs.boundingboxes)
                curbb = textlocs.boundingboxes{zz};
                % Top of the bounding box is below the x-axis
                if max(curbb(:,2)) <= yaxis_peak
                    newlab = {};
                    newlab.text = textlocs.text{zz};
                    newlab.boundingbox = textlocs.boundingboxes{zz};
                    bar_labels{end+1} = newlab;
                end
            end            

            %successful_parse = true;
            for zz=1:length(axis_labels)

                curtext = axis_labels{zz}.text;
                new_label = parseValueLabel(curtext);

                if ~isempty(new_label.value)
                    new_label.left = axis_labels{zz}.left;
                    new_label.right = axis_labels{zz}.right;
                    new_label.center = (new_label.left+new_label.right)/2;
                    new_label.boundingbox = axis_labels{zz}.boundingbox;
                    text_labels{end+1} = new_label;
                end
            end
        end
%%
        % Calculate scale and origin
        scale = 1;
        min_value = 0;

        % Assign values and labels to marks
        bbset = cell(1,length(bar_labels));
        for zz=1:length(bar_labels)
            bbset{zz} = bar_labels{zz}.boundingbox;
        end
        
        tlset = cell(1,length(text_labels));
        for zz=1:length(text_labels)
            tlset{zz} = text_labels{zz}.boundingbox;
        end        
        
        if ~isempty(text_labels)
            origin_x = yaxis_peak;
            
%             yaxis_line.x0 = [bottom_axis yaxis_peak];
%             yaxis_line.d = [size(im,1)-bottom_axis 0];
%             
%             tl_projections = zeros(length(text_labels),2);
%             for zz=1:length(text_labels)
%                 curbb = text_labels{zz}.boundingbox;
%                 centerpt = [(min(curbb(:,1))+max(curbb(:,1)))/2 (min(curbb(:,2))+max(curbb(:,2)))/2];
%                 tl_line.x0 = centerpt;
%                 tl_projections(zz,:) = lineintpoint(xaxis_line);
%             end            
            
            % Assume linear scale
            scalearr = zeros(sum(1:length(text_labels)-1),1);
            ind = 1;
            for t1=1:length(text_labels)
                for t2=t1+1:length(text_labels)
                    scalearr(ind) = (abs(text_labels{t1}.value-text_labels{t2}.value) / abs(text_labels{t1}.center-text_labels{t2}.center));
                    ind = ind + 1;
                end
            end
            scale = median(scalearr);
            if scale == 0
                scale = 1;
            % Estimate the min-value
            else
                min_value_found = false;
                origin = [bottom_axis yaxis_peak];
                if ~isempty(tlset) && length(origin) == 2
                    % Estimate the min-value using the closest label
                    [clp min_dist] = closestLabelFromPoint([bottom_axis yaxis_peak], tlset);
                    closest_origin_label = text_labels{clp};
                    if closest_origin_label.value == 0
                        min_value = 0;
                        min_value_found = true;
                    elseif min_dist < 10
                        min_value = closest_origin_label.value;
                        min_value_found = true;
                    end
                end
                
                if ~min_value_found
                    min_vals = zeros(length(text_labels),1);
                    for t1=1:length(text_labels)
                        min_vals(t1) = text_labels{t1}.value-scale*(text_labels{t1}.center-origin_x);
                    end
                    min_value = median(min_vals);
                end
            end
        end
        fprintf(dataparams,'%s,%2.2f,%2.2f\n',curFile,min_value,scale);

        % Fill out mark structure and output data table
        for m=1:length(marks)
            width = marks{m}.stats.BoundingBox(3);
            marks{m}.value = min_value+scale*width;
            if ~isempty(textlocs) && ~isempty(bbset)
                labind = closestLabel(marks{m}.stats.BoundingBox,bbset);
                %marks{m}.label = textset{labind};
                marks{m}.label = bar_labels{labind}.text;
            else
                marks{m}.label = num2str(m);
            end

            markColor = avgColor(bf.im(marks{m}.stats.SubarrayIdx{:},:));
            % Ensure color ranges from 0 to 1
            if(max(markColor) > 1)
                markColor = markColor/255;
            end            
            
            marks{m}.color = markColor;
            
            marks{m}.id = m;
        end
    %%
        
        handle = imshowlarge(im,true); hold on;
        plot([yaxis_peak yaxis_peak], [1 size(im,1)], '-r');
        %plot([1 size(im,2)], [min_y min_y],'-r');
        %plot([1 size(im,2)], [max_y max_y],'-r');
        if ~isempty(top_axis)
            plot([1 size(im,2)], [top_axis top_axis],'-b');
        end
        if ~isempty(bottom_axis)
            plot([1 size(im,2)], [bottom_axis bottom_axis],'-b');
        end
        for zz=1:length(text_labels)
            estwidth = (text_labels{zz}.value - min_value)/scale + yaxis_peak;
            scatter(estwidth, bottom_axis, 'ro','filled');
        end        
        hold off;
        fprintf(2, 'Saving axis marked image...\n');
        saveImage(handle, curFile, '_axes.png',outpath);
        close;
    end
    
%%    
    % Save the marks
    dataTableFile = fullfile(outpath, strcat(curFile_noext,'_data.txt'));
    datatable = fopen(dataTableFile,'w');
    fprintf(datatable, 'id\tlabel\tvalue\tleft\ttop\twidth\theight\torientation\tr\tg\tb\n'); 
    for m=1:length(marks)
        fprintf(datatable, '%d\t"%s"\t%2.2f\t%2.2f\t%2.2f\t%2.2f\t%2.2f\t%d\t%2.2f\t%2.2f\t%2.2f\n',...
                        marks{m}.id, marks{m}.label, marks{m}.value,...
                        min(marks{m}.stats.Extrema(:,1)),min(marks{m}.stats.Extrema(:,2)),...
                        max(marks{m}.stats.Extrema(:,1))-min(marks{m}.stats.Extrema(:,1)),...
                        max(marks{m}.stats.Extrema(:,2))-min(marks{m}.stats.Extrema(:,2)),...
                        bar_orientation, marks{m}.color(1), marks{m}.color(2), marks{m}.color(3));
    end
    fclose(datatable);
    
    % Save the axes
    axisFile_handle = fullfile(outpath, strcat(curFile_noext,'_axes.txt'));
    axisfile = fopen(axisFile_handle, 'w');
    fprintf(axisfile, 'x1\ty1\tx2\ty2\n');
    fprintf(axisfile, '%d\t%d\t%d\t%d\n', ...
                      bf.axis1.x1, bf.axis1.y1, bf.axis1.x2, bf.axis1.y2);
    fprintf(axisfile, '%d\t%d\t%d\t%d\n', ...
                      bf.axis2.x1, bf.axis2.y1, bf.axis2.x2, bf.axis2.y2);
    fclose(axisfile);
    close all;    

    total_time = toc(total_tic);
    fprintf(2, 'Total time: %6.4f\n', total_time);

    catch ME1
        %beep;
        fprintf(2,'Something went wrong! %s\n', ME1.message);
    end

end

fclose(dataparams);