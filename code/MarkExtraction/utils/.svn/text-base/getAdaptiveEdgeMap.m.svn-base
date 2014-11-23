function [edgelIm, gradIm, nrmGradIm, dirIm, unsuppressedIm, varargout] = getAdaptiveEdgeMap(imGray, textmask, edgesigma, edgestrength, mep, suppress)
    if ~exist('EDGE_SIGMA') || ~exist('EDGE_STRENGTH') || ~exist('MIN_EDGE_PROPORTION')
        thresholds();
    end

    if ~exist('edgesigma', 'var')
        edgesigma = EDGE_SIGMA;
    end
    
    if ~exist('edgestrength', 'var')
        edgestrength = EDGE_STRENGTH;
    end
    
    if ~exist('mep', 'var')
        mep = MIN_EDGE_PROPORTION;
    end
    
    if ~exist('suppress', 'var')
        suppress = NON_MAX_SUPPRESSION;
    end
    
    cursigma = edgesigma;
    curstrength = edgestrength;
    edgelIm = zeros(size(imGray));
    while((sum(sum(edgelIm))/(size(edgelIm,1)*size(edgelIm,2))) < mep && ...
          curstrength > MIN_EDGE_STRENGTH)
        fprintf(2, 'curstrength: %6.6f | curproportion: %6.6f\n', curstrength, ...
                   (sum(sum(edgelIm))/(size(edgelIm,1)*size(edgelIm,2))));
        curstrength = curstrength / 2;
        [edgelIm, gradIm, nrmGradIm, dirIm, unsuppressedIm] = cannyEdgels(imGray,cursigma,curstrength,suppress);
        edgelIm = edgelIm & textmask;
        unsuppressedIm = unsuppressedIm & textmask;
    end
%     
    fprintf(2, 'cursigma: %6.6f | curstrength: %6.6f | proportion: %6.6f\n', ...
            cursigma, curstrength, (sum(sum(edgelIm))/(size(edgelIm,1)*size(edgelIm,2))));
    varargout{1} = cursigma; varargout{2} = curstrength;
    varargout{3} = (sum(sum(edgelIm))/(size(edgelIm,1)*size(edgelIm,2)));
end