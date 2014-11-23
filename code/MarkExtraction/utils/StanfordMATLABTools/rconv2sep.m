function result = rconv2sep(im,rowfilt,colfilt)
% RCONV2SEP: Separable, convolution using reflecting edge handler.
% 
%      result=rconv2sep(im,rowfilt,colfilt)
%
%      im - input image.
%      rowfilt - 1d filter applied to the rows
%      colfilt - 1d filter applied to the cols
%
% Example: foo=rconv2sep(im,[1 4 6 4 1],[-1 0 1]);
%
% DJH '96

rowfilt=rowfilt(:)';
colfilt=colfilt(:);

im = double(im);
rowfilt = double(rowfilt);
colfilt = double(colfilt);

tmp = upConv(im,rowfilt,'reflect1');
result = upConv(tmp,colfilt,'reflect1');
return;

%%% Debug
im=mkImpulse(7);
filter = [1,2,4,2,1];
filter=filter/sum(filter);

res1=rconv2sep(im,filter,filter);
res2=rconv2(im,filter'*filter);
mse(res1,res2)
