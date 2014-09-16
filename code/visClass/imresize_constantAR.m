function B = imresize_constantAR(A, dim)
    if max(max(A)) > 1
        padColor = 255.0;
    else
        padColor = 1.0;
    end
    targX = dim(1);
    targY = dim(2);
    origX = size(A,1);
    origY = size(A,2);
    origAR = origX / origY;

    B = uint8(repmat(padColor, dim));

    if (origAR > 1)     % Wide
        T = imresize(A, [targX NaN]);
        tempY = size(T,2);
        offsetY = floor((targY - tempY) / 2);
        B(:,offsetY+(1:tempY)) = T;
    else                % Tall
        T = imresize(A, [NaN targY]);
        tempX = size(T,1);
        offsetX = floor((targX - tempX) / 2);
        B(offsetX+(1:tempX),:,:) = T;
    end

end