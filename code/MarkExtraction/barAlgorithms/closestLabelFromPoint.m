function [textlocind varargout] = closestLabelFromPoint(point, textloc_boundingboxes)    
    % Break the textloc bounding boxes into lines
    textlines = cell(length(textloc_boundingboxes),1);
    for tl=1:length(textloc_boundingboxes)
        bb = textloc_boundingboxes{tl};
        newlines = cell(4,1);
        
        for line=1:size(bb,1)
            if line < size(bb,1)
                newlines{line}.x0 = bb(line,:);
                newlines{line}.d = bb(line+1,:)-bb(line,:);
            else
                newlines{line}.x0 = bb(line,:);
                newlines{line}.d = bb(1,:)-bb(line,:);
            end
        end
        
        textlines{tl} = newlines;
        
    end
    
    min_dist = Inf;
    for tl=1:length(textlines)
        for tll=1:length(textlines{tl})
            d = pointlinedist(point,textlines{tl}{tll});
            if d < min_dist
                min_dist = d;
                textlocind = tl;
            end
        end
    end
    
    varargout{1} = min_dist;
end