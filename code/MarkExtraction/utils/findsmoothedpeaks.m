function [peaks varargout] = findsmoothedpeaks(pixel_diff, gFilt, gxFilt, gxxFilt, sderiv_thresh, normalize)

    if ~exist('sderiv_thresh','var')
        sderiv_thresh = -0.01;
    end

    if ~exist('normalize','var')
        normalize = false;
    end
    
    if sderiv_thresh > 0
        sderiv_thresh = -sderiv_thresh;
    end
    
    DISPLAY_ON = false;
    
    lf = length(gxxFilt);
    pd_0 = conv(pixel_diff, gFilt);
    pd_0 = pd_0(floor(lf/2)+1:length(pd_0)-floor(lf/2));
    
    if normalize
        pd_0 = pd_0/max(pd_0);
    end
    
    pd_1 = conv(pixel_diff, gxFilt);
    pd_1 = pd_1(floor(lf/2)+1:length(pd_1)-floor(lf/2));
    pd_2 = conv(pixel_diff, gxxFilt);
    pd_2 = pd_2(floor(lf/2)+1:length(pd_2)-floor(lf/2));

    %% Find zero crossings
    peaks = [];
    for zz=2:length(pd_1)
        if (pd_1(zz-1) >=0 && pd_1(zz) <= 0) && ...
           (pd_2(zz) < sderiv_thresh)
            peaks(end+1) = zz;
        end
    end

    if DISPLAY_ON
        figure;
        subplot(4,1,1); plot(pixel_diff);
        subplot(4,1,2); plot(pd_0);
        subplot(4,1,3); plot(pd_1);
        subplot(4,1,4); plot(pd_2);
        
        subplot(4,1,1); hold on; scatter(zc, zeros(length(zc), 1), 8);
        subplot(4,1,3); hold on; scatter(zc, zeros(length(zc), 1), 8);
    end
    
    varargout{1} = pd_0; varargout{2} = pd_1; varargout{3} = pd_2;

end