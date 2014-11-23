function color = avgColor(im, pixelList)
   
    if ~exist('pixelList', 'var')
        pixelList = 1:(size(im,1)*size(im,2));
    end

    rPlane = im(:,:,1);
    gPlane = im(:,:,2);
    bPlane = im(:,:,3);
    
    color = zeros(1,3);
    color(1) = sum(sum(rPlane(pixelList)))/length(pixelList);
    color(2) = sum(sum(gPlane(pixelList)))/length(pixelList);
    color(3) = sum(sum(bPlane(pixelList)))/length(pixelList);
end