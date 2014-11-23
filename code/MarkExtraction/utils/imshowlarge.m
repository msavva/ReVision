function h = imshowlarge(im,hidden)
    if ~exist('hidden','var')
        hidden = false;
    end
    if hidden
        hidden = 'off';
    else
        hidden = 'on';
    end
    figure('Position', [0 0 1440 900],'visible',hidden);
    h = imshow(im, 'InitialMagnification', 'fit');
end