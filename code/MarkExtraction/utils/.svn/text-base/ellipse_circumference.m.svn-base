%%
% Use Ramanujan's approximation to compute the ellipse circumference

function circ = ellipse_circumference(params)
    % a = major axis radius
    % b = minor axis radius
    a = max(params(3), params(4));
    b = min(params(3), params(4));
    
    x = ((a-b)/(a+b));
    
    circ = pi*(a+b)*(1+((3*x^2)/(10+sqrt(4-3*x^2))));
end