function diff = adjacentcolordiff(et)
    % Color image
    if size(et, 3) > 1
        etcol = zeros(3, length(et));
        etcol(1,:) = et(:,:,1);
        etcol(2,:) = et(:,:,2);
        etcol(3,:) = et(:,:,3);

        et_shift = zeros(size(et));
        et_shift(1,1:length(et)-1,:) = et(1,2:length(et),:);
        et_shift(1,length(et),:) = et(1,1,:);

        et_shift_col = zeros(3, length(et_shift));
        et_shift_col(1,:) = et_shift(:,:,1);
        et_shift_col(2,:) = et_shift(:,:,2);
        et_shift_col(3,:) = et_shift(:,:,3);

        %et_shift_col = sqrt(sum(et_shift_col.^2, 1));

        diff = sqrt(sum((et_shift_col - etcol).^2,1));
    % Grayscale image
    else
        et_shift = zeros(size(et));
        et_shift(1,1:length(et)-1) = et(1,2:length(et));
        et_shift(1,length(et)) = et(1,1);
        
        diff = sqrt(sum((et_shift-et).^2,1));
    end
end