function seginds = findsegments(r)
    rdiff = r(2:length(r))- r(1:length(r)-1);
    seginds = [];
    startind = 1;
    endind = -1;
    for q=1:length(rdiff)
        if rdiff(q) > 0
            startind = q+1;
        elseif rdiff(q) < 0 
            endind = q;
            seginds(end+1, :) = [startind endind];
            startind = -1;
            endind = -1;
        end
    end
    if endind == -1 && startind ~= -1
        if ~isempty(find(r(startind:end)))
            seginds(end+1, :) = [startind length(r)];
        end
    end
end