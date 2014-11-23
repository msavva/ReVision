function points = getEllipsePoints(params, angle)

    orientation = params(5);

    xvals = params(1)+params(3)*cos(angle)*cos(orientation) - ...
           params(4)*sin(angle)*sin(orientation);
    yvals = params(2)+params(3)*cos(angle)*sin(orientation) + ...
           params(4)*sin(angle)*cos(orientation);
       
    points = [xvals(:) yvals(:)];

end