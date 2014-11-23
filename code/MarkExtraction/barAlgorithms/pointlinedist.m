function [d varargout] = pointlinedist(p,l1)
    x1 = l1.x0(1); y1 = l1.x0(2);
    t = l1.x0+l1.d;
    x2 = t(1); y2 = t(2);

    % Project point p onto the line
    l1_norm = [-l1.d(2) l1.d(1)];
    l1_norm = l1_norm/norm(l1_norm,2);
    pline.x0 = p;
    pline.d = l1_norm;
    
    x3 = pline.x0(1); y3 = pline.x0(2);
    t = pline.x0+pline.d;
    x4 = t(1); y4 = t(2);
    
    denom = det([x1-x2 y1-y2; x3-x4 y3-y4]);    
  
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
    
    l1_pt = [int_x int_y];
    if a1 < 0
        l1_pt = [x1 y1];
    elseif a1 > 1
        l1_pt = [x2 y2];
    end

    d = norm(l1_pt-p,2);
    
    varargout(1) = {l1_pt};
end