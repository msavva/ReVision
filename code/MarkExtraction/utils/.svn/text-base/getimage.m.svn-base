function Img = getimage(Filename)
%GETIMAGE  Read an image given a filename
%   V = GETIMAGE(FILENAME) where FILENAME is an image file.  The image is
%   returned either as an MxN double matrix for a grayscale image or as an
%   MxNx3 double matrix for a color image, with elements in [0,1].

% Pascal Getreuer 2008-2009

% Read the file
[Img,Map,Alpha] = imread(Filename);
%Img = double(Img);

if ~isempty(Map)    % Convert indexed image to RGB
    Img = Img + 1;
    Img = reshape(cat(3,Map(Img,1),Map(Img,2),Map(Img,3)),size(Img,1),size(Img,2),3);
else
    Img = double(Img)/255;  % Rescale to [0,1]
    if max(max(max(Img))) <=  1/255
        Img = Img*255;
    end
end
