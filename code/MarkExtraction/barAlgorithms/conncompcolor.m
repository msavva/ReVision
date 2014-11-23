% Return the labeled connected components for an RGB image
%
% TODO: optimize using union-find or similar data structure
function [labelIm, N, varargout] = conncompcolor(I)
%     if ~exist('bgcolor', 'var')
%         bgcolor = [255 255 255];
%     end

    delta = 10;

    if(max(I(:)) < 1+eps)
        I = I * 255;
    end
    
    % Initialize output variables
    labelIm = repmat(0,size(I(:,:,1)));
    curLabel = 0;
    
%     bgcolor = reshape(bgcolor, [1 1 3]);
%     bgcolorIm = repmat(bgcolor, size(I(:,:,1)));
%     mseIm = mse(I, bgcolorIm);
%     bgIm = zeros(size(mseIm));
%     bgIm(mseIm < delta) = 1;
    
    I = double(I);
    %bgIm = double(bgIm);
    
    equivalenceMatrix = {};
    
    tFirst = tic;
    % First pass: record equivalences and assign labels
    % i,j is row,col number
    
    tenrow_time = 0;
    for i=1:size(I,1)
        if mod(i, 25) == 0
            fprintf(2, '%d rows complete. Time: %6.6f seconds\n', i, tenrow_time);
        end
        tTenrow = tic;
        for j=1:size(I,2)
            
            %sprintf('i: %i, j: %i, bgIm: %i', i, j, bgIm(i,j))
            
            % Not a background color pixel
            %if(bgIm(i,j) == 0)
                
                %sprintf('i: %i, j: %i', i, j)
                
                neighbors = getNeighbors(i, j);
                smallestLabel = Inf;
                labels = [];
                
                for k=1:length(neighbors)
                    curInds = neighbors{k};
                    if(sum(curInds < 1) == 0 && curInds(2) <= size(I,2))
                        %sprintf('curneighbor - i, j: %i %i', curInds(1), curInds(2))
                        curNeighbor = I(curInds(1),curInds(2),:);
                        %if(bgIm(curInds(1),curInds(2)) == 0 && mse(curNeighbor, I(i,j,:)) < delta)
                        if(mse(curNeighbor, I(i,j,:)) < delta)
                            labels = [labels labelIm(curInds(1), curInds(2))];
                            if(labelIm(curInds(1),curInds(2)) < smallestLabel)
                                smallestLabel = labelIm(curInds(1),curInds(2));
                            end
                        end
                    end
                end
                
                labels = unique(labels);
                
                % Update equivalence matrix if there is more than one label
                if length(labels) > 1
                    assigned = false;
                    for z=1:length(equivalenceMatrix)
                        if(find(equivalenceMatrix{z} == min(labels)))
                            % Concatenate two rows together
                            equivalenceMatrix{z} = unique(cat(2,equivalenceMatrix{z},labels));
                            assigned = true;
                            break;
                        end
                    end
                    if(~assigned)
                        equivalenceMatrix{end+1} = labels;
                    end
                end
                
                if(smallestLabel ~= Inf)
                    %disp('smallestLabel not inf');
                    labelIm(i,j) = smallestLabel;
                else
                    %disp('smallestLabel is inf');
                    curLabel = curLabel + 1;
                    labelIm(i,j) = curLabel;
                end
            %end
        end
        tenrow_time = tenrow_time + toc(tTenrow);
    end
    firstpass_time = toc(tFirst);
    fprintf(2, 'first pass: %6.6f seconds\n', firstpass_time);
    
    % Second pass: reassign labels based on the equivalence matrix
    eTime = tic;
    for k=1:length(equivalenceMatrix)
        curEq= equivalenceMatrix{k};
        for j=1:length(curEq)
            if(~isempty(find(labelIm == curEq(j))))
                labelIm(labelIm == curEq(j)) = min(curEq);
            end
        end
    end
    equivalence_time = toc(eTime);
    fprintf(2, 'time to assign equivalence labels: %6.6f seconds\n', equivalence_time);
    
    % Assign labels in order
    oTime = tic;
    uniqueLabels = unique(labelIm);
    for k=1:length(uniqueLabels)
        labelIm(labelIm == uniqueLabels(k)) = k-1;
    end
    order_time = toc(oTime);
    fprintf(2, 'time to reassign labels in order: %6.6f seconds\n', order_time);
    
    N = max(max(labelIm));
    
    varargout{1} = firstpass_time; varargout{2} = equivalence_time;
end

function neighbors = getNeighbors(i,j)
    neighbors = { [i,j-1]; [i-1,j-1]; [i-1,j]; [i-1, j+1] };
end

function out = mse(I1, I2)
    out = sqrt(sum((I1-I2).^2,3));
end