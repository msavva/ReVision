function [gFilt gxFilt gxxFilt] = genGaussians(sigma)
    sigmaSqr = sigma*sigma;
    gFiltSize = 2 * round(3.0 * sigma) + 1;
    x = [1:gFiltSize] - round((gFiltSize+1)/2);
    gFilt = exp(- x .* x / (2.0*sigmaSqr));
    gFilt = gFilt/ sum(gFilt(:));
    gxFilt = - (x / sigmaSqr) .* gFilt;
    gxxFilt = ((x / sigmaSqr).^2 - 1.0/sigmaSqr) .* gFilt;
end