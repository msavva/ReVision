function textlocind = closestLabel(boundingbox, textloc_boundingboxes)
    
    textlocind = -1;

    lines = cell(4,1);
    y = boundingbox(1); x=boundingbox(2);
    h=boundingbox(3); w=boundingbox(4);
    %x = boundingbox(1); y=boundingbox(2);
    %w=boundingbox(3); h=boundingbox(4);
    
    % Break the bounding box into four lines
    lines{1}.x0 = [x y];
    lines{1}.d = [x+w y]-lines{1}.x0;
    lines{2}.x0 = [x+w y];
    lines{2}.d = [x+w y+h]-lines{2}.x0;
    lines{3}.x0 = [x+w y+h];
    lines{3}.d = [x y+h]-lines{3}.x0;
    lines{4}.x0 = [x y+h];
    lines{4}.d = [x y]-lines{4}.x0;
    
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
        for bbl=1:length(lines)
            for tll=1:length(textlines{tl})
                d = linedist(lines{bbl},textlines{tl}{tll});
                if d < min_dist
                    min_dist = d;
                    textlocind = tl;
                end
            end
        end
    end
end