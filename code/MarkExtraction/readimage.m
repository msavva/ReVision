%%
% Reads in an color or grayscale image and returns the original and
% grayscale version of the images. maxDim is an optional parameter that
% specifies the maximum dimension of the image.
function [origim grayim] = readimage(imFile, maxDim)

    origim = getimage(imFile);
    
    %if max(max(max(origim))) == 255
    %    origim = double(origim)/255;
    %end
    
    if size(origim,3) > 1
        grayim = rgb2gray(origim);
    else   
        grayim = origim;
    end
    if exist('maxDim', 'var')
        if max(size(grayim)) > maxDim
            fprintf(2,'Resizing image...\n');
            scale = maxDim/max(size(grayim));
            origGrayMax = max(max(grayim));
            origColorMax = max(max(max(origim)));
            
            grayim = imresize(grayim, scale, 'nearest');
            grayim = min(grayim, origGrayMax);
            grayim = max(grayim, 0);
            
            origim = imresize(origim, scale, 'nearest');
            origim = min(origim, origColorMax);
            origim = max(origim, 0);
        end
    end
end