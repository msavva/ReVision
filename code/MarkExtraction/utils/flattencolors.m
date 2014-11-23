%%
% Converts a three-plane color image into a 2D image, where each row is one
% of the planes, or a grayscale image into a row vector.

function flatcolor = flattencolors(color)
    if size(color,3) > 1
        flatcolor = zeros(3, size(color,1)*size(color,2));
        fc_1 = color(:,:,1);
        fc_2 = color(:,:,2);
        fc_3 = color(:,:,3);
        flatcolor(1,:) = fc_1(:);
        flatcolor(2,:) = fc_2(:);
        flatcolor(3,:) = fc_3(:);
    else
        flatcolor = color(:)';
    end
end