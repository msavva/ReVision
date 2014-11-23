function [mask textlocs] = gettextmask(textlocFile,im)
    DISPLAY_ON = false;

    try
        [imsize textlocs] = readTextFile(textlocFile);
        rects = [textlocs.x*imsize(1) textlocs.y*imsize(2) ...
                 textlocs.w*imsize(1) textlocs.h*imsize(2) -textlocs.theta];   
             
        if DISPLAY_ON
            h = figure; imshow(im); hold on;
        end

        % Hack to fix differences in image size from the text extractor to
        % the actual image size
        im_offset = [size(im,2)-imsize(1); size(im,1)-imsize(2)];
        imsize = imsize+im_offset;
        
        boundingboxes = struct([]);
        mask = zeros(imsize(2), imsize(1));
        for i=1:size(rects,1)
            
            % Hack to fix differences in image size from the text extractor to
            % the actual image size
            rects(i,2) = rects(i,2)+im_offset(2)/2;
            rects(i,1) = rects(i,1)+im_offset(1)/2;
            
            currect = [rects(i,2) rects(i,1); ...
                       rects(i,2) rects(i,1)+rects(i,3); ...
                       rects(i,2)+rects(i,4) rects(i,1)+rects(i,3); ...
                       rects(i,2)+rects(i,4) rects(i,1); ...
                       rects(i,2) rects(i,1)];
            trans = repmat([(rects(i,2)+rects(i,4)/2) (rects(i,1)+rects(i,3)/2)], [size(currect,1), 1]);
            rot = [cos(rects(i,5)) -sin(rects(i,5)); ...
                  sin(rects(i,5)) cos(rects(i,5))];
            currect = (rot*(currect-trans)')'+trans;
            
            boundingboxes{end+1} = currect(1:4,:);
            
            mask = mask | poly2mask(currect(:,2), currect(:,1), imsize(2), imsize(1));
            if DISPLAY_ON
                plot(currect(:,2), currect(:,1), 'LineWidth', 2, 'Color', 'red'); pause;
            end
        end
        textlocs.boundingboxes = boundingboxes;

        if DISPLAY_ON
            hold off;
            close(h);
        end

        mask = 1-mask;
    catch ME1
        fprintf(2, '%s\n', ME1.message);
        mask = NaN;
        textlocs = [];
    end
end