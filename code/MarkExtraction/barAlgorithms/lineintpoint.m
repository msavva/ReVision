function intpoint = lineintpoint(l1, l2)
    x1 = l1.x0(1); y1 = l1.x0(2);
    t = l1.x0+l1.d;
    x2 = t(1); y2 = t(2);
    
    x3 = l2.x0(1); y3 = l2.x0(2);
    t = l2.x0+l2.d;
    x4 = t(1); y4 = t(2);
    
    denom = det([x1-x2 y1-y2; x3-x4 y3-y4]);
    
    % Parallel lines
    if denom < eps
        intpoint = [NaN NaN];
    else
        numa = det([x1 y1; x2 y2]);
        numb = det([x3 y3; x4 y4]);
        int_x = det([numa x1-x2; numb x3-x4]) / denom;
        int_y = det([numa y1-y2; numb y3-y4]) / denom;
        
        intpoint = [int_x int_y];
    end
end