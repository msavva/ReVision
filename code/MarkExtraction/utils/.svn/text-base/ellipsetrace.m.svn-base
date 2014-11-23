% Output a trace of pixels by walking the image at the specified ellipse
%
% TODO: optimize using union-find or similar data structure
function [out_trace text_trace varargout] = ellipsetrace(I, params, textmask, npoints, interp, disp)

   if ~exist('disp','var')
        disp = false;
   end
   
   if ~exist('npoints','var')
        npoints = 1000;
   end
   
   if ~exist('interp','var')
       interp = true;
   end

   npoints = npoints-1;
   
   phi = 0:2*pi/npoints:2*pi;
   %phi = 0:pi/8/npoints:pi/8;
   
   % Orientation
   orientation = params(5);
   
   xvals = params(1)+params(3)*cos(phi)*cos(orientation) - ...
           params(4)*sin(phi)*sin(orientation);
   yvals = params(2)+params(3)*cos(phi)*sin(orientation) + ...
           params(4)*sin(phi)*cos(orientation);

%    inTextBox = zeros(length(xvals),1);
%    if ~isempty(textlocs)
%        for textboxInd = 1:length(textlocs.boundingboxes)
%            inTextBox = inTextBox | inpolygon(xvals, yvals, textlocs.boundingboxes{textboxInd}(:,1), ...
%                                                            textlocs.boundingboxes{textboxInd}(:,2));
%        end
%    end
%    
%    textTraceInds = sub2ind(size(I), yvals(inTextBox), xvals(inTextBox));

   discreteXVals = floor(xvals);
   discreteXVals = max(discreteXVals, 1);
   discreteXVals = min(discreteXVals, size(I, 2));
   
   discreteYVals = floor(yvals);
   discreteYVals = max(discreteYVals, 1);
   discreteYVals = min(discreteYVals, size(I, 1));
   inds = sub2ind(size(I), discreteYVals, discreteXVals);
   
   text_trace = [];
   if ~isempty(textmask)
    text_trace = textmask(inds);
   end

   if ~interp
       if size(I,3) > 1
           I1 = I(:,:,1);
           I2 = I(:,:,2);
           I3 = I(:,:,3);
           out_trace = zeros(1, npoints+1, 3);
           out_trace(:,:,1) = I1(inds);
           out_trace(:,:,2) = I2(inds);
           out_trace(:,:,3) = I3(inds);
           clear I1 I2 I3
       else
           out_trace = zeros(1, npoints+1);
           out_trace(:,:,1) = I(inds);
       end
   else
        % Interpolate
        if size(I,3) > 1
           I1 = I(:,:,1);
           I2 = I(:,:,2);
           I3 = I(:,:,3);
           I1int = interp2(I1, xvals, yvals);
           I2int = interp2(I2, xvals, yvals);
           I3int = interp2(I3, xvals, yvals);
           out_trace = zeros(1, npoints+1, 3);
           out_trace(:,:,1) = I1int;
           out_trace(:,:,2) = I2int;
           out_trace(:,:,3) = I3int;
           clear I1 I2 I3
        else
            out_trace(:,:,1) = interp2(I, xvals, yvals);
        end
   end
       

   if disp
       size(I)
      figure; imshow(I); hold on;
      plot(xvals,yvals,'Color','red','LineWidth',1); hold on;
        plot(params(1),params(2),'--ro','LineWidth',2,'MarkerEdgeColor','red', ...
                                   'MarkerFaceColor','red','MarkerSize',5);
      hold off;
   end
   
   varargout{1} = xvals; varargout{2} = yvals;
end