function [labeledIm junctures] = tracelines(ti)
    labeledIm = uint8(ti);
    labeledIm(labeledIm == 1) = -1;
    
    junctures = {};
    
    currow_segs = findsegments(ti(1,:));
    curlabel = 1;
    for zz=1:size(currow_segs,1)
        labeledIm(1, currow_segs(zz,1):currow_segs(zz,2)) = curlabel;
        curlabel = curlabel+1;
    end
    
    
    % Label each subsequent row given the 
    for r=2:size(ti,1)
        % Label the row given the previous row. Pixels are assigned the
        % previous row's labels.
        prevrow_segs = currow_segs;
        currow_segs = findsegments(ti(r,:));
        
%         if r > 40
%             r
%             prevrow_segs
%             currow_segs
%         end
        
        for zz=1:size(currow_segs,1)
            adj_segs = [];
            
            % Check adjacency to this segment with all segments in the
            % previous row.
            for qq=1:size(prevrow_segs,1)
                if(adjacentsegs(currow_segs(zz,:), prevrow_segs(qq, :)))
                    adj_segs(end+1, :) = prevrow_segs(qq,:);
                    break;
                end
            end
            while qq < size(prevrow_segs,1)
                qq = qq + 1;
                if(adjacentsegs(currow_segs(zz,:), prevrow_segs(qq,:)))
                    adj_segs(end+1, :) = prevrow_segs(qq,:);
                else
                    break;
                end
            end
            
%             if r > 40
%                 zz
%                 adj_segs
%                 pause;
%             end
            
            % This segment is adjacent to segments in the previous row
            if ~isempty(adj_segs)
                labels = labeledIm(r-1,adj_segs(:,1));
                labeledIm(r, currow_segs(zz,1):currow_segs(zz,2)) = labels(1);
                % More than one section was adjacent. Add these labels to a
                % juncture, and use the first label.                
                if size(adj_segs,1) > 1
                    junctures{end+1} = labels(:)';
                end
            else
                labeledIm(r, currow_segs(zz,1):currow_segs(zz,2)) = curlabel;
                curlabel = curlabel+1;
            end
        end
    end
end