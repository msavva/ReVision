function drawEllipse(params, color, varargin)
   npoints = 100;
   phi = 0:2*pi/npoints:2*pi;
   %phi = 0:pi/4/npoint:pi/4;
   
   width = 1;
   if(~isempty(varargin))
       width = varargin{1};
   end
   
   % Orientation
   orientation = params(5);
   
   xvals = params(1)+params(3)*cos(phi)*cos(orientation) - ...
           params(4)*sin(phi)*sin(orientation);
   yvals = params(2)+params(3)*cos(phi)*sin(orientation) + ...
           params(4)*sin(phi)*cos(orientation);
   
   plot(xvals,yvals,'Color',color,'LineWidth',width); hold on;
   plot(params(1),params(2),'--ro','LineWidth',2,'MarkerEdgeColor',color, ...
                                   'MarkerFaceColor',color,'MarkerSize',3);
end