function [edgelIm, gradIm, nrmGradIm, dirIm, enoughGradAmp] = cannyEdgels(im, sigma, minStrength, doNonMaxSuppression)
% [edgelIm nrmIm dirIm] = cannyEdgels(IM, SIGMA, MINSTRENGTH)
%
% Compute the Canny edgels using the 1st derivatives of a Gaussian filter.
% This does not perform hysteresis thresholding.
%
% PARAMS:  
%  IM the 2D gray level input image
%  SIGMA (optional, default = 1.0) the std dev of the Gaussian filter.
%  MINSTRENGTH (optional, default = 1.0)  The threshold
%    value for the minimum gradient magnitude.  
%
% OUTPUT:  
%  edgelIm  binary image marking edgel locations.
%  nrmGradIm  image of norm of gradient image (NOT thresholded)
%  dirIm    direction of gradient image (in terms of theta/pi), so
%           directions are periodic in the interval (0, 2]. 
%           dirIm is NaN at pixels for which nrmGradIm < minStrength.  

% ADJ 9/01.

%%%%%%%%%%%%%%% Check parameters  %%%%%%%%%%%%%%%%%%%%%%%

%%% Fill in default parameter values, and correct bogus parameter values.
if ~exist('sigma', 'var')  %% sigma is the std dev for the Gaussian filters
  sigma = 1.0;
elseif (sigma < 0.5/3)
  fprintf(2, 'User specified sigma = %e is too small.\n', sigma); 
  sigma = 0.5000001/3;
  fprintf(2, 'Using sigma = %f instead\n', sigma);
end

if ~exist('minStrength', 'var')  
  minStrength = 1; 
elseif minStrength<=0.0
  fprintf(2,'Invalid edge strength tolerance %e\n', minStrength); 
  minStrength = eps;
  fprintf(2,'Using edge strength tolerance = %e\n', minStrength);
end

if ~exist('doNonMaxSuppression','var')
    doNonMaxSuppression = true;
end

%%%%%%%%%%%%%%% Build Filter Kernels %%%%%%%%%%%%%%%%%%%%%%%

%%% Build Gaussian filter masks, along with derivatives.
%%% The second derivative is not used yet.
sigmaSqr = sigma*sigma;
gFiltSize = 2 * round(3.0 * sigma) + 1;
x = [1:gFiltSize] - round((gFiltSize+1)/2);
gFilt = exp(- x .* x / (2.0*sigmaSqr));
gFilt = gFilt/ sum(gFilt(:));
gxFilt = - (x / sigmaSqr) .* gFilt;
gxxFilt = ((x / sigmaSqr).^2 - 1.0/sigmaSqr) .* gFilt;

%%%%%%%%%%%%%%% Do separable convolutions %%%%%%%%%%%%%%%%%%%

gradIm = zeros([size(im), 2]);
gradIm(:,:,1) = rconv2sep(im, gxFilt, gFilt);
gradIm(:,:,2) = rconv2sep(im, gFilt, gxFilt);

%%%%%%%%%%%%%%% Gradient Magnitude Image %%%%%%%%%%%%%%%%%%%
 
nrmGradIm = sum(gradIm .* gradIm, 3).^0.5;

%%%%%%%%%%%%%%% Enforce Magnitude Threshold %%%%%%%%%%%%%%%%%%%

dirIm = zeros(size(im));
enoughGradAmp = (nrmGradIm >= minStrength);

%%%%%%%%%%%%%%% Compute Direction Image %%%%%%%%%%%%%%%%%%%
dxIm = gradIm(:,:,1);
dyIm = gradIm(:,:,2);
dirIm(enoughGradAmp) = atan2(dyIm(enoughGradAmp), dxIm(enoughGradAmp))/pi;
%% Note this value of dirIm is periodic in the range [-1, 1]

%% Branch cut at 0 degrees (edges look like "dark | bright")
cut = 0;
while (sum(sum(dirIm > cut+2.0)) > 0)
  dirIm(dirIm>cut+2.0) = dirIm(dirIm>cut+2.0) - 2.0;
end
while (sum(sum(dirIm <= cut)) > 0)
  dirIm(dirIm<=cut) = dirIm(dirIm<=cut) + 2.0;
end
%%% dirIm now in the range (0, 2.0]

%% Insert NaN in regions for which the gradient is too small.
dirIm(~enoughGradAmp) = NaN;

if doNonMaxSuppression
    edgelIm = zeros(size(im));
    
    %%%%%%%%%%%%%%% Setup for  Non-max Edgel Suppression %%%%%%%%%%%%%%%%
    %% Set alpha to be the distance from the origin to the
    %% point that is one pixel to the right, and pi/8 degrees off horizontal.
    theta = pi/8;
    alpha = 0.5/sin(theta);
    
    %% We scale the edge normal by this amount, and then round it.
    %% This causes the rounded value to be the closest direction
    %% vector on the discrete grid for a 3x3 neighbourhood.
    %%

    %%%%%% Matrix version of non-max suppression implementation

    %% Make the norm of the gradient artificially large at
    %% pixels that are below the minStrength threshold.
    maxGrad = max(nrmGradIm(:));
    nrm = nrmGradIm .* enoughGradAmp + (~enoughGradAmp) * (10.0 * maxGrad);

    %% Compute the shifts in the x and y directions.
    qx = round(alpha * dxIm ./ nrm);
    qy = round(alpha * dyIm ./ nrm);

    %% Interpolate the norm of the gradient image at the
    %% pixels given by the shifts.
    [x y] = meshgrid(1:size(im,2), 1:size(im,1));
    ampNear = interp2(x, y, nrmGradIm, x+qx, y+qy, 'nearest');
    ampNear(isnan(ampNear)) = 0.0;

    %% Do non-max supression in this shift direction.
    ok = enoughGradAmp;
    ok = ok & (ampNear < nrmGradIm);

    %% Do non-max supression in the opposite shift direction.
    ampNear = interp2(x, y, nrmGradIm, x-qx, y-qy, 'nearest');
    ampNear(isnan(ampNear)) = 0.0;
    edgelIm = ok & (ampNear < nrmGradIm);
else
    edgelIm = enoughGradAmp;
end
%% End: cannyEdgels
