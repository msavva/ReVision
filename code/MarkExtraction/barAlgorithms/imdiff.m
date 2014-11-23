function diff = imdiff(im, orient)

    flatim = zeros(3*size(im,1), size(im,2));
    flatim(1:3:end-2,:) = im(:,:,1);
    flatim(2:3:end-1,:) = im(:,:,2);
    flatim(3:3:end,:) = im(:,:,3);

    % x-gradient (for horizontal bars)
    if orient == 'x'
        flatimshift = zeros(size(flatim));
        flatimshift(1:end-3,:) = flatim(4:end,:);
        flatimshift(end-2:end,:) = flatim(1:3,:);
    elseif orient == 'y'
        flatimshift = zeros(size(flatim));
        flatimshift(:,1:end-1) = flatim(:,2:end);
        flatimshift(:,end) = flatim(:,1);
    else
        fprintf(2, 'Error! orient must be either "x" or "y"\n');
        return;
    end
    
    flatdiff = (flatimshift-flatim).^2;
    diff = zeros(size(im,1),size(im,2),3);
    diff(:,:,1) = flatdiff(1:3:end-2,:);
    diff(:,:,2) = flatdiff(2:3:end-1,:);
    diff(:,:,3) = flatdiff(3:3:end,:);
    diff = sqrt(sum(diff,3));    
    
end