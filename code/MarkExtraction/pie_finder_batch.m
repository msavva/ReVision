% Prototype to find slices in a pie chart.

%% Find pie using RANSAC (doesn't work with exploded slices)
clear all
close all

setupPath();
thresholds();

IMAGE_READ_DIRECTORY = 1;
IMAGE_READ_SUBSET = 2;

%%%%%%%%%%%%%%%%%%%%%%
% Flags

CURVE_OPTIMIZE = false;
IMAGE_READ_MODE = IMAGE_READ_DIRECTORY;

% [chart_files fullpath outpath] = getImagesInPath('Testing/UMDpies');
% textpath =
% '/Users/nick/Documents/Projects/VisionCharts/VisInterp/corpus/UMDImages/p
% ie';
% cachepath = fullfile(pwd, 'cache');
%

%% Use a text file that specifies the files
if IMAGE_READ_MODE == IMAGE_READ_SUBSET

    chart_subset_file = 'Testing/largepies/large_pies.txt';
    %chart_subset_file = 'Testing/UMDpies/umd_pies.txt';
    chart_path = 'Visualization/PieChart';
    %chart_path = 'UMDImages/pie';
    chart_files = getchartfiles(chart_subset_file);
    outpath = fullfile(pwd, 'Testing/largepies');
    %outpath = fullfile(pwd, 'Testing/UMDpies');

    for i=1:length(chart_files)
        chart_files{i} = strcat(chart_files{i}, '.png');
    end

elseif IMAGE_READ_MODE == IMAGE_READ_DIRECTORY

    %% Use a directory

    % Large corpus paths
    chart_path = 'bing_pie chart/UIST11_assumptions';
    outpath = fullfile(pwd, 'Testing/bing_pie chart');

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
% empty if using UMD.
cacheManager.cacheprefix = '';

errors = fopen(fullfile(outpath, 'ellipse_errors.txt'), 'a');
times = fopen(fullfile(outpath, 'times.txt'), 'a');
edgeparams = fopen(fullfile(outpath,'edgeparams.txt'),'a');
usededgemap = fopen(fullfile(outpath,'usededgemap.txt'),'a');
ellipseparams = fopen(fullfile(outpath,'ellipseparams.txt'),'a');

fprintf(times, 'file,edgel_time,curve_tracing_time,ellipse_RANSAC_time,line_RANSAC_time\n');
fprintf(errors,'file,best_error\n');
fprintf(edgeparams, 'file,imsigma,imstrength,improportion,cropsigma,cropstrength,cropproportion\n');
fprintf(usededgemap, 'file,usededegmap\n');
fprintf(ellipseparams,'file,x,y,minor,major,orient\n');

% Initialize smoothing filters
sigma = MIN_SEGMENT_LENGTH;
[gFilt gxFilt gxxFilt] = genGaussians(sigma);
sigma_2 = 3;
[gFilt2 gxFilt2 gxxFilt2] = genGaussians(sigma_2);

% Initialize simple second derivative
[gFilt_raw gxFilt_raw gxxFilt_raw] = genGaussians(0.3);

normalize = false;

randomseed = 0;
rand('twister',randomseed);
% s = RandStream('mt19937ar','Seed', randomseed);
% RandStream.setDefaultStream(s)

% Specify a custom set of charts
%chart_files = {'pie_126.png','pie_127.png'};
%chart_files = {'00156.jpg', '00179.gif', '00014.png', '00131.Gif', '00232.png', '00073.jpg', '00084.png', ...
%               '00099.JPG', '00240.PNG', '00108.gif', '00120.gif'};


%%
for i=1:length(chart_files)
%for i=107
    total_tic = tic;

    try
        
    % Read in the image_files
    curFile = chart_files{i};
    %curFile = 'pie_103.png';
    extInd = strfind(curFile, '.');
    curFile_noext = curFile(1:extInd-1);
    
    textlocFile = fullfile(fullpath, 'text', strcat(curFile_noext,'.txt'));
    
    markPropertiesFile = fullfile(outpath, strcat(curFile_noext,'_marks.txt'));
    
    fprintf(2,'Processing %s...\n',curFile);
    fprintf(errors, '%s,',curFile);
    fprintf(times, '%s,',curFile);
    fprintf(edgeparams,'%s,',curFile);
    fprintf(usededgemap,'%s,',curFile);

%%    
    % Handle both 3-channel and indexed images. Don't resize!
    curIm = getimage(fullfile(fullpath,curFile));   
    
    % Flag for color image
    COLOR = size(curIm,3) > 1;
    
    % Convert image range to 0-1
    if(max(max(max(curIm))) > 1)
        curIm = double(curIm)/double(255);
    else
        curIm = double(curIm);
    end
    
    % Check for filtered image
    filteredIm_file = strcat(curFile_noext,'_filtered_',sprintf('%d', BFILTER_WIDTH), '_', ...
         sprintf('%2.2f',BFILTER_DOMAINSIGMA),'_',sprintf('%2.2f',BFILTER_RANGESIGMA),'.mat');    
    
    % Check the cache
    loadPath = cacheManager.getLoadPath(filteredIm_file);
    % In cache
    if ~isempty(loadPath)
        fprintf(2,'Loading cached image...\n');
        load(loadPath);
    % Not in cache: filter the image
    else
        fprintf(2,'Applying bilateral filter...\n');
        curIm = bfilter2(curIm, BFILTER_WIDTH, [BFILTER_DOMAINSIGMA BFILTER_RANGESIGMA]);
        save(cacheManager.getSavePath(filteredIm_file, cacheManager.LOCAL_CACHE), 'curIm', 'BFILTER_WIDTH', ...
                                            'BFILTER_DOMAINSIGMA', 'BFILTER_RANGESIGMA');
    end
    
    if COLOR
        curImGray = rgb2gray(curIm);
    else
        curImGray = curIm;
    end
    
    handle = figure('visible','off'); imshow(curIm);
    fprintf(2,'Saving filtered image...\n');
    saveImage(handle, curFile, '_filteredim.png', outpath);
    close
    
    % Retrieve the text locations, if available
    % Remove the text if text locations exist
    if exist(textlocFile, 'file')
        [textmask textlocs] = gettextmask(textlocFile,curIm);
        clear textLocFile
    else
        % No text file exists, so don't mask anything
        textmask = ones(size(curIm,1),size(curIm,2));
        textlocs = [];
    end    
    
    % Output text mask debug images if text mask exists
    if ~isempty(textlocs)
        handle = figure('visible','off'); imshow(textmask);
        fprintf(2,'Saving text mask...\n');
        saveImage(handle, curFile, '_textmask.png', outpath);            

        handle = figure('visible','off');
        imshow(curIm);
        hold on;
        for rectInd=1:length(textlocs.boundingboxes)
            currect = [textlocs.boundingboxes{rectInd} ; textlocs.boundingboxes{rectInd}(1,:)];
            plot(currect(:,2), currect(:,1), 'LineWidth', 1, 'Color', 'red');
        end
        hold off;
        fprintf(2,'Saving text mask overlaid on image...\n');
        saveImage(handle, curFile, '_textmaskoverlay.png', outpath);       
    end    
    
%%
    % Find the Canny edgels
    % Check for cached image
    adaptiveEdgel_file = strcat(curFile_noext,'_adaptiveedgels_',sprintf('%2.2f', EDGE_SIGMA), '_', ...
         sprintf('%2.2f',EDGE_STRENGTH),'_',sprintf('%2.4f',MIN_EDGE_PROPORTION),'.mat');
     
    % Check the cache
    loadPath = cacheManager.getLoadPath(adaptiveEdgel_file);
     
    % Generate the image if it's not cached
    if isempty(loadPath)
        fprintf(2,'Finding Canny edgels...\n'); 
        etTic = tic;
        [edgelIm, gradIm, nrmGradIm, dirIm, unsuppressedIm, cursigma,curstrength,proportion] = getAdaptiveEdgeMap(curImGray, textmask);
        edgel_time = toc(etTic);
        save(cacheManager.getSavePath(adaptiveEdgel_file, cacheManager.LOCAL_CACHE),...
            'edgelIm', 'dirIm', 'gradIm', 'nrmGradIm', 'unsuppressedIm', 'cursigma', 'curstrength', 'proportion','edgel_time');
    % Retrieve the cached image
    else
        fprintf(2,'Loading cached edge map...\n');
        load(loadPath);
    end 
    
    fprintf(times,'%6.4f,',edgel_time);
    fprintf(edgeparams,'%6.6f,%6.6f,%6.6f,',cursigma,curstrength,proportion);
%%   
    if(~isempty(find(edgelIm)))

        handle = figure('visible','off'); imshow(edgelIm);
        fprintf(2,'Saving text-removed edgel image...\n');
        saveImage(handle, curFile, '_notextedgel.png', outpath);
        close
        
        handle = figure('visible','off'); imshow(unsuppressedIm);
        fprintf(2,'Saving text-removed unsuppressed edgel image...\n');
        saveImage(handle, curFile, '_notextunsuppressededgel.png', outpath);
        close
        
        if CURVE_OPTIMIZE
        
            % Close the edge map before curve tracing to attempt to keep contours
            % joined
            edgelClosed = imclose(edgelIm, curveTracing_closingElement);            
            
            % Check for cached curves
            curve_file = strcat(curFile_noext,'_curves_',sprintf('%2.2f', EDGE_SIGMA), '_', ...
                 sprintf('%2.2f',EDGE_STRENGTH),'_',sprintf('%2.4f',MIN_EDGE_PROPORTION),'.mat');
            
            loadPath = cacheManager.getLoadPath(curve_file);
             
            % Generate the curves if they aren't cached
            if isempty(loadPath)
                fprintf(2, 'Curve tracing..\n');
                curvetic = tic;
                %[curves, dim] = getCurves(edgelIm, dirIm.*pi, nrmGradIm,curveTracing_maxDirDiff);
                [curves, toosmall, dim] = getCurves(edgelClosed, dirIm.*pi, nrmGradIm, curveTracing_maxDirDiff); 
                curve_tracing_time = toc(curvetic);
                save(cacheManager.getSavePath(curve_file, cacheManager.LOCAL_CACHE),...
                    'curves', 'toosmall', 'dim', 'curve_tracing_time');
            % Retrieve the cached image
            else
                fprintf(2,'Loading cached curves...\n');
                load(loadPath);
            end
            fprintf(times, '%6.4f,',curve_tracing_time);

            %%%%%%
            % Create an image showing each curve and ellipses fit to that curve
            %
            fprintf(2,'Fitting ellipses...\n');
            tic
            [allEllipses boundingBoxArea] = fitEllipsesToCurves(curves, size(edgelIm));
            toc

            %%%%    
            % Classify curves into curves, text, and lines
            %
            fprintf(2,'Classifying curves...\n');
            tic
            % Discard ellipses that aren't related
            lineInds = [];
            textInds = [];
            curveInds = [];
            areaThreshold = (size(edgelIm,1)*size(edgelIm,2))/500;
            for k=1:length(allEllipses)
                cx = allEllipses{k}(1);
                cy = allEllipses{k}(2);
                rx = allEllipses{k}(3);
                ry = allEllipses{k}(4);

                % If a radius is NaN, indicates a straight line
                if(isnan(rx))
                    lineInds(end+1) = k;
                % If one ellipse axis is significantly longer than another, the curve
                % is likely a line.
                elseif(min(rx,ry)/max(rx,ry) < lineAspectThreshold)
                    lineInds(end+1) = k;
                % Eliminate text by using small bounding boxes and circular
                % aspect
                % ratios
                elseif(boundingBoxArea(k) < areaThreshold && min(rx,ry)/max(rx,ry) > textAspectThreshold)
                    textInds(end+1) = k;
                else
                    curveInds(end+1) = k;
                end
            end
            toc

            fprintf(2,'Creating curve classification image...\n');
            otherColor = 4;
            curveColor = 3;
            lineColor = 2;
            textColor = 1;
            cltIm = zeros(size(edgelIm));
            for k=1:length(lineInds)
                cltIm(curves{lineInds(k)}) = lineColor;
            end
            for k=1:length(textInds)
                cltIm(curves{textInds(k)}) = textColor;
            end
            for k=1:length(curveInds)
                cltIm(curves{curveInds(k)}) = curveColor;
            end
            for k=1:length(toosmall)
                for j=1:length(toosmall{k})
                    cltIm(toosmall{k}(j)) = otherColor;
                end
            end
            curveTypeMap = zeros(4,3);
            curveTypeMap(1,:) = [1 1 1];
            curveTypeMap(2,:) = [1 0 0];
            curveTypeMap(3,:) = [0 1 0];
            curveTypeMap(4,:) = [0 0 1];
            curveTypeMap(5,:) = [1 0 1];
            handle = figure('visible','off'); imshow(cltIm/max(max(cltIm))); colormap(curveTypeMap);

            fprintf(2,'Saving curve image...\n');
            saveImage(handle, curFile, '_curves.png', outpath);
            close
        
        end
%%
        % Check for cached circle
        fittedCircle_file = strcat(curFile_noext,'_fittedCircle_',sprintf('%d',randomseed),'_',...
                sprintf('%d',CIRCLE_RANSAC_ITERATIONS),sprintf('%2.2f', EDGE_SIGMA), '_', ...
                sprintf('%2.2f',EDGE_STRENGTH),'_',sprintf('%2.4f',MIN_EDGE_PROPORTION),'.mat');
        
        loadPath = cacheManager.getLoadPath(fittedCircle_file);
            
        if isempty(loadPath)
            fprintf(2,'Use RANSAC to find the best circle (%d iterations)...\n',CIRCLE_RANSAC_ITERATIONS);
            ransac_tic = tic;
            
            if CURVE_OPTIMIZE
                curveImage = zeros(size(edgelIm));
                for k=1:length(curveInds)
                    curveImage(curves{curveInds(k)}) = 1;
                end
                [circ_best_params circ_best_inliers circ_best_error circ_best_initial_points best_circles] = circleRANSAC(curveImage,CIRCLE_RANSAC_ITERATIONS);
            else
                curveImage = edgelIm;
                [circ_best_params circ_best_inliers circ_best_error circ_best_initial_points best_circles] = circleRANSAC(curveImage,CIRCLE_RANSAC_ITERATIONS);
            end

            ransac_time = toc(ransac_tic);
            save(cacheManager.getSavePath(fittedCircle_file, cacheManager.SHARED_CACHE),...
                'circ_best_params', 'circ_best_inliers', 'circ_best_error', ...
                'circ_best_initial_points', 'best_circles', 'ransac_time');        
        else
            fprintf(2,'Loading cached circle...\n');
            load(loadPath);            
        end
        
        % If the best circle's error is still above a threshold, try and
        % fit an ellipse.
        if circ_best_error > 0.8
            % Check for cached ellipse
            fittedEllipse_file = strcat(curFile_noext,'_fittedEllipse_',sprintf('%d',randomseed),'_',...
                    sprintf('%d',ELLIPSE_RANSAC_ITERATIONS),sprintf('%2.2f', EDGE_SIGMA), '_', ...
                    sprintf('%2.2f',EDGE_STRENGTH),'_',sprintf('%2.4f',MIN_EDGE_PROPORTION),'.mat');
                
            loadPath = cacheManager.getLoadPath(fittedEllipse_file);
            % Generate the ellipse if it hasn't been cached
            if isempty(loadPath)
                fprintf(2,'Use RANSAC to find the best ellipse...\n');
                ransac_tic = tic;

                if CURVE_OPTIMIZE
                    curveImage = zeros(size(edgelIm));
                    for k=1:length(curveInds)
                        curveImage(curves{curveInds(k)}) = 1;
                    end
                    [best_params best_inliers best_error best_initial_points best_ellipses] = ellipseRANSAC(curveImage,ELLIPSE_RANSAC_ITERATIONS);
                else
                    curveImage = edgelIm;
                    [best_params best_inliers best_error best_initial_points best_ellipses] = ellipseRANSAC(curveImage,ELLIPSE_RANSAC_ITERATIONS);
                end

                ransac_time = ransac_time+toc(ransac_tic);
                save(cacheManager.getSavePath(fittedEllispe_file, cacheManager.SHARED_CACHE),...
                    'best_params', 'best_inliers', 'best_error', 'best_initial_points', 'best_ellipses', 'ransac_time');
            % Retrieve the cached image
            else
                fprintf(2,'Loading cached ellipse...\n');
                load(loadPath);
            end

            if best_error > circ_best_error
                best_params = circ_best_params;
                best_inliers = circ_best_inliers;
                best_error = circ_best_error;
                best_initial_points = circ_best_initial_points;
            end
        else
            best_params = circ_best_params;
            best_inliers = circ_best_inliers;
            best_error = circ_best_error;
            best_initial_points = circ_best_initial_points;
            clear best_ellipses;
        end

        fprintf(ellipseparams,'%s,%6.2f,%6.2f,%6.2f,%6.2f,%6.2f\n',curFile,best_params(1),best_params(2),best_params(3),best_params(4),best_params(5));
        fprintf(times, '%6.4f,',ransac_time);
        
        % Output the errors of the ellipse
        errorString = sprintf('%6.4f\n',best_error);
        fprintf(2, errorString);
        fprintf(errors, errorString);
        handle = figure('visible','off'); imshow(curImGray); hold on;
        
        % Draw the top 10 circles (if they exist)
        if exist('best_circles','var')
            for circind=2:size(best_circles)
                drawEllipse(best_circles(circind,1:5),'green',1);
            end
        end
        
        % Draw the top 10 ellipses (if they exist)
        if exist('best_ellipses','var')
            for ellipseind=2:size(best_ellipses)
                drawEllipse(best_ellipses(ellispseind,1:5),'blue',1);
            end        
        end
        
        % Draw the best ellipse and the inliers associated with that
        % ellipse
        drawEllipse(best_params,'red',3);
        scatter(best_initial_points(:,1), best_initial_points(:,2), 'MarkerEdgeColor', 'g','MarkerFaceColor', 'g', 'SizeData', 25);
        hold off;
        fprintf(2,'Saving best ellipse from all ellipses...\n');
        saveImage(handle, curFile, '_bestellipse.png', outpath);
        close all
%%
        
        %%%%%%%% After cropping to the ellipse, circle-walk

        %%%%%%%%% Use the walked circle to identify borders/transition
        %%%%%%%%% points, and mark out lines
        %ellipseradiidividers = [1.05 1.1 1.15 1.25 1.5 1.75 2 2.25 2.5 2.75 3];
        
        fprintf(2, 'Walking ellipses...\n');
        
        % Ellipse radii to trace
        ellipseradiidividers = 1.05:0.01:5;
        
        % Output the ellipse traces
        traceimage = zeros(length(ellipseradiidividers), ELLIPSE_SUBDIVISIONS,3);
        traceimage_debug = zeros(length(ellipseradiidividers), ELLIPSE_SUBDIVISIONS,3);
        trace_textmask = ones(length(ellipseradiidividers),ELLIPSE_SUBDIVISIONS);
        
        ellipse_peaks_votes = zeros(1, ELLIPSE_SUBDIVISIONS);
        pixeldiff_im = zeros(length(ellipseradiidividers), ELLIPSE_SUBDIVISIONS);
        
        for z=1:length(ellipseradiidividers)
            ellipse_to_trace = best_params;
            ellipse_to_trace(3:4) = ellipse_to_trace(3:4)/ellipseradiidividers(z);
            [et tt] = ellipsetrace(curIm, ellipse_to_trace, textmask, ELLIPSE_SUBDIVISIONS,true);
            if COLOR
                traceimage(z,:,:) = et;
                traceimage_debug(z,:,:) = et;
                et = colorspace('rgb->lab',et);
            else
                traceimage(z,:,:) = repmat(et,[1 1 3]);
                traceimage_debug(z,:,:) = repmat(et,[1 1 3]);
            end
            trace_textmask(z,:) = trace_textmask(z,:) & tt;
                
            pixel_diff = adjacentcolordiff(et);
            pixel_diff = pixel_diff.*tt;
            
            pixeldiff_im(z,:) = pixel_diff;
            
%            [ellipse_trace_peaks etp_0 etp_1 etp_2] = findsmoothedpeaks(pixel_diff, gFilt, gxFilt, gxxFilt);
%            [ellipse_trace_peaks etp_0 etp_1 etp_2] = findsmoothedpeaks(pixel_diff, gFilt_raw, gxFilt_raw, gxxFilt_raw,0);
            ellipse_trace_peaks = find(pixel_diff > mean(pixel_diff));

%             if z==50
%                 et0 = pixel_diff;
%                 et0_0 = etp_0;
%                 et0_1 = etp_1;
%                 et0_2 = etp_2;
%                 et0_p = ellipse_trace_peaks;
%             end
            
            
%%%% FIND PEAKS BY GROUPING            
%             ellipse_trace_peaks = find(pixel_diff > LAB_COLOR_THRESH);
%             % Group together peaks that are closer than MIN_SEGMENT_LENGTH. Given the
%             % values of the peaks in ascending order, traverse the array backwards and
%             % store peaks that are close together in peaks_in_segment. Once two peaks
%             % that are far apart are found, replace the current set of close peaks with
%             % their average. z == 2 covers the border case where the first two peaks
%             % are close together.
% 
%             MIN_SEGMENT_LENGTH = min(10,500/max(ellipse_to_trace(3:4)));
%             peaks_in_segment = [];
%             for q=length(ellipse_trace_peaks):-1:2
%                 peaks_in_segment(end+1) = q;
%                 if ellipse_trace_peaks(q)-ellipse_trace_peaks(q-1) < MIN_SEGMENT_LENGTH
%                     if q == 2
%                         peaks_in_segment(end+1) = 1;
%                     end
%                 else
%                     ellipse_trace_peaks(peaks_in_segment) = (max(ellipse_trace_peaks(peaks_in_segment))+min(ellipse_trace_peaks(peaks_in_segment)))/2;
%                     peaks_in_segment = [];
%                 end
%             end
%             ellipse_trace_peaks(peaks_in_segment) = (max(ellipse_trace_peaks(peaks_in_segment))+min(ellipse_trace_peaks(peaks_in_segment)))/2;
%             ellipse_trace_peaks = unique(ellipse_trace_peaks);
            
            ellipse_peaks_votes(round(ellipse_trace_peaks)) = ellipse_peaks_votes(round(ellipse_trace_peaks))+1;
            
            for zz=1:length(ellipse_trace_peaks)
                traceimage_debug(z, round(ellipse_trace_peaks(zz)), :) = [1 0 0];
            end
        end
        
        peakvotes = 1-ellipse_peaks_votes/max(ellipse_peaks_votes);
        if COLOR
            peakvotes = repmat(peakvotes, [1 1 3]);
        end
%         for zz=1:10
%             traceimage_debug(end+1,:,:) = peakvotes;
%         end
 
%% DEBUG OUTPUT
        fprintf(2, 'Saving pixel difference image...\n');
        handle = imshowlarge(pixeldiff_im/max(max(pixeldiff_im)),true);
        saveImage(handle, curFile, '_pixeldiffim.png', outpath);
        close;
        
        pdi_l = pixeldiff_im;
        pdi_l = pdi_l > 10;
        fprintf(2, 'Saving peak threshold image...\n');
        handle = imshowlarge(pdi_l,true);
        saveImage(handle, curFile, '_peakthresh.png', outpath);
        close;

        %%
%%%% BIN PEAK VOTES
%         BIN_SIZE = 3;
%         epv = ellipse_peaks_votes;
%         for zz=1:BIN_SIZE-mod(length(ellipse_peaks_votes),BIN_SIZE)
%             epv(end+1) = 0;
%         end
%         epv = reshape(epv, [BIN_SIZE length(epv)/BIN_SIZE ]);
%         epv = sum(epv,1);
%         epv = repmat(epv, [BIN_SIZE 1]);
%         epv = reshape(epv, [1 size(epv,2)*BIN_SIZE]);
%         epv = epv(1:length(ellipse_peaks_votes));
%         ellipse_peaks_votes = epv;
%% Check for votes at either end of the trace, within the smoothing width
%% (3 at the moment)
        circular_offset = 0;
        if ~isempty(find(ellipse_peaks_votes(1:6) > 0)) || ...
                ~isempty(find(ellipse_peaks_votes(length(ellipse_peaks_votes)-6:length(ellipse_peaks_votes))) > 0)
            % Find the first row of six zeros, and offset the trace so that
            % the row of six zeros straddles the end and beginning of the trace.
            for oo=1:length(ellipse_peaks_votes)-6
                if sum(ellipse_peaks_votes(oo:oo+6)) == 0
                     circular_offset = oo+3;
                     break;
                end
            end
        end

        %circular_offest = 0;
        ellipse_peaks_votes = circshift(ellipse_peaks_votes,[0 -circular_offset]);
        %
        [new_peaks sv_0 sv_1 sv_2] = findsmoothedpeaks(ellipse_peaks_votes, gFilt2, gxFilt2, gxxFilt2,0);
        epv_std = std(ellipse_peaks_votes);
        epv_m = mean(ellipse_peaks_votes);
        std_up = epv_m+epv_std;
        new_peaks = new_peaks(ellipse_peaks_votes(new_peaks) > std_up);
        
%% DEBUG OUTPUT
        fprintf(2, 'Saving vote plot...\n');
        handle = figure('visible','off'); plot(ellipse_peaks_votes);
        hold on;
            plot([1 length(ellipse_peaks_votes)], [epv_m epv_m], '-r');
            plot([1 length(ellipse_peaks_votes)], [epv_m+epv_std epv_m+epv_std], '-r');
        hold off;
        saveImage(handle, curFile, '_votes.png', outpath);
        close

        fprintf(2, 'Saving smoothed vote plot...\n');
        handle = figure('visible','off');
        subplot(3,1,1);
            hold on;
            plot(sv_0);
            % Mean
            plot([1 length(ellipse_peaks_votes)], [epv_m epv_m], '-r');
            plot([1 length(ellipse_peaks_votes)], [epv_m+epv_std epv_m+epv_std], '-r');
            hold off;
        subplot(3,1,2);
            plot(sv_1);
        subplot(3,1,3);
            plot(sv_2);
        saveImage(handle, curFile, '_smoothedvotes.png', outpath);
        close    
       
%         if COLOR
%             sv_0 = repmat(sv_0, [1 1 3]);
%         end
%         for zz=1:10
%             traceimage_debug(end+1,:,:) = sv_0;
%         end

        fprintf(2, 'Saving trace debug image...\n');
        handle = imshowlarge(traceimage_debug,true);
        saveImage(handle, curFile, '_tracedebug.png', outpath);
        close
        
        fprintf(2, 'Saving trace image...\n');
        handle = imshowlarge(traceimage,true);
        saveImage(handle, curFile, '_trace.png', outpath);
        close

        fprintf(2, 'Saving trace text mask...\n');
        handle = imshowlarge(trace_textmask,true);
        saveImage(handle, curFile, '_tracetextmask.png', outpath);
        close
        %%
%%%%% Find lines in the edge map of the ellipse trace image and use that
%%%%% instead of votes

%         if COLOR
%             trace_edgels = cannyEdgels(rgb2gray(traceimage_debug), 1, 0.03);
%         else
%             trace_edgels = cannyEdgels(traceimage_debug, 1, 0.03);
%         end
%         labeledIm = tracelines(trace_edgels);
%         imshowlarge(trace_edgels);
%         
%         fprintf(2, 'Saving labeled edges...\n');
%         handle = imshowlarge(label2rgb(labeledIm));
%         saveImage(handle, curFile, '_labelededges.png', outpath);
%         close
% 
%         % Find long-spanning labels
%         labels = unique(labeledIm);
%         % Ignore the background (0) label
%         labels = labels(labels > 0);
%         numEllipses = size(labeledIm,1);
%         lineIm = labeledIm;
% 
%         for q=1:length(labels)
%             labelInds = find(labeledIm == labels(q));
%             [labelI labelJ] = ind2sub(size(labeledIm), labelInds);
%             if (max(labelI)-min(labelI)+1) < numEllipses-5
%                 lineIm(lineIm == labels(q)) = 0;
%             end
%         end
%         fprintf(2, 'Saving labeled lines...\n');
%         handle = imshowlarge(label2rgb(lineIm));
%         saveImage(handle, curFile, '_labeledlines.png', outpath);
%         close
%         
%         labels = unique(lineIm);
%         labels = labels(labels > 0);
%         topLocs = zeros(1, length(labels));
%         bottomLocs = zeros(1, length(labels));
%         avgLocs = zeros(1, length(labels));
%         for q=1:length(labels)
%             [i j] = find(lineIm == labels(q));
%             temp = j(i == min(i));
%             topLocs(q) = temp(1);
%             temp = j(i == max(i));
%             bottomLocs(q) = temp(1);
% 
%             avgLocs(q) = topLocs(q);
%         end
%         
%         new_peaks = avgLocs;

%%%%% Group vote peaks
%         MIN_VOTES = 50;
%         new_peaks = find(ellipse_peaks_votes > MIN_VOTES);
%         % Group together peaks that are closer than MIN_SEGMENT_LENGTH. Given the
%         % values of the peaks in ascending order, traverse the array backwards and
%         % store peaks that are close together in peaks_in_segment. Once two peaks
%         % that are far apart are found, replace the current set of close peaks with
%         % their average. z == 2 covers the border case where the first two peaks
%         % are close together.
%         peaks_in_segment = [];
%         MIN_SEGMENT_LENGTH = 8;
%         for q=length(new_peaks):-1:2
%             peaks_in_segment(end+1) = q;
%             if new_peaks(q)-new_peaks(q-1) < MIN_SEGMENT_LENGTH
%                 if q == 2
%                     peaks_in_segment(end+1) = 1;
%                 end
%             else
%                 new_peaks(peaks_in_segment) = (max(new_peaks(peaks_in_segment))+min(new_peaks(peaks_in_segment)))/2;
%                 peaks_in_segment = [];
%             end
%         end
%         new_peaks(peaks_in_segment) = (max(new_peaks(peaks_in_segment))+min(new_peaks(peaks_in_segment)))/2;
%         new_peaks = unique(new_peaks);
        
        new_peaks = new_peaks + circular_offset;
        new_peaks(new_peaks > ELLIPSE_SUBDIVISIONS) = new_peaks(new_peaks > ELLIPSE_SUBDIVISIONS) - ELLIPSE_SUBDIVISIONS;
        angles = 2*pi*new_peaks/ELLIPSE_SUBDIVISIONS;
        epts = getEllipsePoints(best_params, angles);

        handle = imshowlarge(curIm,true); hold on;
        %drawEllipse(best_params, 'red');
        for z=1:size(epts,1)
            lineparams = {};
            %lineparams.x0 = [ellipse_to_trace(1) ellipse_to_trace(2)];
            lineparams.x0 = [best_params(1) best_params(2)];
            lineparams.d = epts(z,:) - lineparams.x0;
            drawLine(lineparams);
            %scatter(epts(:,1), epts(:,2), 'MarkerEdgeColor', 'g','MarkerFaceColor', 'g', 'SizeData', 25);
        end
        hold off;
        fprintf(2,'Saving found lines..\n');
        saveImage(handle, curFile, '_bordertracedlines.png', outpath);
        close all
        
        %% Output angles and slice colors for each sector

        fprintf(2, 'Locating points on the pie closest to text...\n');
        
        
        if exist('textlocs','var') && ~isempty(textlocs)
            % Determine text proximity to ellipse
            bbflat = zeros(length(textlocs.boundingboxes)*4,2);
            for z=1:length(textlocs.boundingboxes)
                bbflat((4*(z-1)+1):(4*z),:) = textlocs.boundingboxes{z};
            end
            bbflat = round(bbflat+0.5);
            % Switch columns
            bbflat(:, [1,2]) = bbflat(:, [2,1]);

            [bberr bbproj] = Residuals_ellipse_vector1(bbflat, best_params);
            bberr_reshape = reshape(bberr, 4, []);
            [minerr minind] = min(bberr_reshape);
            for z=1:length(minind)
                minind(z) = minind(z)+4*(z-1);
            end
            associated_points = bbflat(minind, :);
            associated_proj = bbproj(minind, :);
            closest_proj = bbproj(minind, :);
            center = repmat([best_params(1) best_params(2)], [size(closest_proj,1) 1]);
            cp_c = closest_proj-center;
            proj_angles = atan2(cp_c(:,2), cp_c(:,1));
            proj_angles(proj_angles < 0) = 2*pi+proj_angles(proj_angles < 0);
            proj_angles = proj_angles - repmat(best_params(5), [length(proj_angles) 1]);
            proj_angles(proj_angles > 2*pi) = proj_angles(proj_angles > 2*pi)-2*pi;

            [proj_angles order] = sort(proj_angles);
            minerr = minerr(order);
            text_marks = textlocs.text;
            text_marks = text_marks(order);
        else
            proj_angles = [];
        end
        
        fprintf(2, 'Extract slices...\n');
        tic
        
        [angles, order] = sort(angles);
        epts = epts(order, :);
        
        centerpt = [best_params(1) best_params(2)];
        slices = struct([]);
        
        if ~isempty(textlocs)
            curlabels = textlocs.text(:)';
            curbb = textlocs.boundingboxes;
        end

        dataTableFile = fullfile(outpath, strcat(curFile_noext,'_data.txt'));
        datatable = fopen(dataTableFile,'w');
        fprintf(datatable, 'id\tlabel\tvalue\tbegin\tend\n');        

        for z=1:length(angles)
            % Store the angle and the bounding lines
            slices(z).line1.x0 = centerpt;
            slices(z).line2.x0 = centerpt;
            slices(z).id = z;
            
            if z < length(angles)
                slices(z).begAngle = angles(z);
                slices(z).endAngle = angles(z+1);
                slices(z).line1.d = epts(z,:) - centerpt;
                slices(z).line2.d = epts(z+1,:) - centerpt;
                slices(z).percent = 100*(angles(z+1)-angles(z))/(2*pi);
            else
                slices(z).begAngle = angles(z);
                slices(z).endAngle = angles(1)+2*pi;
                slices(z).line1.d = epts(z,:) - centerpt;
                slices(z).line2.d = epts(1,:) - centerpt;
                slices(z).percent = 100*(2*pi-angles(z)+angles(1))/(2*pi);
            end
            
            % Use the closest text label
            if exist('textlocs','var') && ~isempty(textlocs) && ~isempty(curlabels)
                % Approximate the ellipse arc as a set of lines
                begAngle = slices(z).begAngle;
                endAngle = slices(z).endAngle;
                
                sampleAngles = begAngle:2*pi*(5/360):endAngle;
                if sampleAngles(end) ~= endAngle
                    sampleAngles(end+1) = endAngle;
                end
                samplePts = getEllipsePoints(best_params, sampleAngles);
                samplePts = samplePts(:,2:-1:1);
                sampleLines = {};
                for zz=1:size(samplePts,1)-1
                    newline = {};
                    newline.x0 = samplePts(zz,:);
                    newline.d = samplePts(zz+1,:)-samplePts(zz,:);
                    sampleLines{end+1} = newline;
                end

                % Find the closest text label
                tind = closestLabelFromLines(sampleLines,curbb);
                slices(z).label = curlabels{tind};

                curbb = cat(2,curbb(1:tind-1),curbb(tind+1:end));
                curlabels = cat(2,curlabels(1:tind-1),curlabels(tind+1:end));
            else
                slices(z).label = sprintf('%d',slices(z).id);
            end
        end

        for z=1:length(slices)
            fprintf(datatable, '%d\t"%s"\t%2.2f\t%2.3f\t%2.3f\n', slices(z).id,slices(z).label,slices(z).percent,...
                                slices(z).begAngle, slices(z).endAngle);
        end
        
        fclose(datatable);
        
        toc
    else
        fprintf(2,'No edges found: not enough gradient amplitude?\n');
    end
    
    catch ME1
        %beep;
        fprintf(2,'Something went wrong! %s', ME1.message);
    end
    
    total_time = toc(total_tic);
    fprintf(2, 'total_time = %6.4f\n', total_time);
end

fclose(errors);
fclose(times);
fclose(edgeparams);
fclose(usededgemap);
fclose(ellipseparams);