function [d varargout] = linedist(l1, l2)
    x1 = l1.x0(1); y1 = l1.x0(2);
    t = l1.x0+l1.d;
    x2 = t(1); y2 = t(2);
    
    x3 = l2.x0(1); y3 = l2.x0(2);
    t = l2.x0+l2.d;
    x4 = t(1); y4 = t(2);
    
    denom = det([x1-x2 y1-y2; x3-x4 y3-y4]);
    
    % Parallel lines
    if denom < eps
        l1_pt = [x1 y1];
        l2_pt = [x3 y3];
    else
        numa = det([x1 y1; x2 y2]);
        numb = det([x3 y3; x4 y4]);
        int_x = det([numa x1-x2; numb x3-x4]) / denom;
        int_y = det([numa y1-y2; numb y3-y4]) / denom;

        a1 = [int_x int_y]-[x1 y1];
        if l1.d(2) == 0
            a1= a1*[1/l1.d(1); 1];
        else
            a1= a1*[1; (1-l1.d(1))/l1.d(2)];
        end

        a2 = [int_x int_y]-[x3 y3];
        if l2.d(2) == 0
            a2= a2*[1/l2.d(1); 1];
        else
            a2= a2*[1; (1-l2.d(1))/l2.d(2)];
        end

        l1_pt = [int_x int_y];
        if a1 < 0
            l1_pt = [x1 y1];
        elseif a1 > 1
            l1_pt = [x2 y2];
        end

        l2_pt = [int_x int_y];
        if a2 < 0
            l2_pt = [x3 y3];
        elseif a2 > 1
            l2_pt = [x4 y4];
        end
    end

    d = norm(l1_pt-l2_pt,2);
    
    varargout(1) = {l1_pt}; varargout(2) = {l2_pt};
end