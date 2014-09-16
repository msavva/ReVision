classdef TextFeatureSet < handle
    
    properties
        fname;                  % filename
        imageDim;               % Image dimensions
        boxes;                  % box data in normalized [0,1] coordinates
        boxesPx;                % box data in pixel dimensions
        labels;                 % box labels
        text;                   % box text
        mask;                   % Mask with text regions
    end % Instance properties
    
    methods(Static)
        % Draw rectangle rotated by theta around its center
        function XY = rect(p)
            x=p(1); y=p(2); w=p(3); h=p(4); t=p(5);
            xv=[x x+w x+w x x];
            yv=[y y y+h y+h y];
            cx = x + w/2;
            cy = y + h/2;
            T1 = eye(3);
            T2 = eye(3);
            T1(1:2,3) = [cx; cy];
            T2(1:2,3) = -[cx; cy];
            R = [cos(t) -sin(t) 0;sin(t) cos(t) 0; 0 0 1];
            M = T1*R*T2;
            XY = M*[xv; yv; ones(1,5)];
            XY = bsxfun(@rdivide, XY, XY(3,:));
            XY = XY(1:2,:);
        end
        function ret = percentageNumericalStrings(str)
            if nargin == 0
                ret = 0;
            else
                res = regexpi(str,'[0-9]+');
                ret = sum(~cellfun(@isempty,res)) / numel(str);
            end

        end
        function superimposeGrid(im, L)
            imshow(imsc(im));
            hold on;
            M = size(im,1);
            N = size(im,2);
            span = round(M / 2^L);
            for k = [1:span:M M]
                x = [1 N];
                y = [k k];
                plot(x,y,'Color','w','LineStyle','-');
                plot(x,y,'Color','k','LineStyle',':');
            end
            
            for k = [1:span:N N]
                x = [k k];
                y = [1 M];
                plot(x,y,'Color','w','LineStyle','-');
                plot(x,y,'Color','k','LineStyle',':');
            end
        end
    end
    
    methods
        % Constructor
        function tfs = TextFeatureSet(filename)
            % Load from file
            if nargin == 1
                tfs.fname = filename;
                % Read header
                fid = fopen(filename, 'r');
                dim = textscan(fid, '%u %u',1,'HeaderLines',1,'Delimiter',',');
                tfs.imageDim = double([dim{:}]);
                fclose(fid);
                
                % Read text boxes
                fid = fopen(filename, 'r');
                lines = textscan(fid, '%n %n %n %n %n %q %q','HeaderLines',3,'Delimiter',',');
                fclose(fid);
                
                % Save data
                B = [lines{1:5}];
                C = lines{6};
                T = lines{7};
                tfs.boxes = B;
                B(:,[1 3]) = B(:,[1 3]) .* tfs.imageDim(1);
                B(:,[2 4]) = B(:,[2 4]) .* tfs.imageDim(2);
                tfs.boxesPx = B;
                tfs.labels = C;
                tfs.text = T;
                
                % Compute mask
                tfs.mask = zeros(tfs.imageDim(2), tfs.imageDim(1));
                for i=1:size(B,1)
                    XY = TextFeatureSet.rect(B(i,:));
                    tfs.mask = tfs.mask | poly2mask(XY(1,:), XY(2,:), tfs.imageDim(2), tfs.imageDim(1));
                end
            end
        end
        
        % Construct a feature vector using text box properties at resolution
        % DIM, with 2^l subdivisions along each axis
        function Q = featureVector(tfs, dim, l)
            n = length(tfs);
            if n>1
                Q = [];
                for i=1:n
                    Q = [Q; featureVector(tfs(i),dim, l)'];
                end
                return;
            elseif n==0
                Q = [];
                return;
            else
                histX = 0.05:.1:.95;
                histXD = linspace(0,sqrt(2),10);
                histXT = -pi:pi/6:pi;
                NNumStr = TextFeatureSet.percentageNumericalStrings(tfs.text);
                Whist = hist(tfs.boxes(:,3), histX);                               Whist = Whist / norm(Whist);
                Hhist = hist(tfs.boxes(:,4), histX);                               Hhist = Hhist / norm(Hhist);
                Thist = hist(tfs.boxes(:,5), histXT);                              Thist = Thist / norm(Thist);
                Cx = tfs.boxes(:,1)+(tfs.boxes(:,3)./2);
                Cy = tfs.boxes(:,2)+(tfs.boxes(:,4)./2);
                PXhist = hist(Cx, histX);                                           PXhist = PXhist / norm(PXhist);
                PYhist = hist(Cy, histX);                                           PYhist = PYhist / norm(PYhist);
                [POhist, PDhist] = tfs.pairwiseOrientationsAndDistances([Cx Cy], histXT, histXD);
                
                span = round(max(dim) / (2^l));
                txtMask = ~imresize_constantAR(uint8(~tfs.mask), dim);
                MaskVec = blockproc(txtMask, [span span], @(b) sum(sum((b.data))) ./ numel(b.data));
                MaskVec = MaskVec(:);
%                 MaskVec = MaskVec - mean(MaskVec);
%                 MaskVec = MaskVec / sqrt(var(MaskVec));
%                 Q = MaskVec;
                Q = [NNumStr MaskVec' Whist Hhist Thist PXhist PYhist PDhist POhist]';
            end
        end
        % Return feature vector labels
        function Qlabels = featureVectorLabels(tfs, l)
                Qlabels = ['NNumStr'; strnum('txtMaskVec',2^(2*l)); strnum('Whist',10); ...
                           strnum('Hhist',10); strnum('Thist',13); strnum('PXhist',10); ...
                           strnum('PYhist',10); strnum('PDhist',10); strnum('POhist',13)];
        end

        % Works in normalized image coordinates for now so orientations are a
        % bit skewed
        function [histTN, histDN] = pairwiseOrientationsAndDistances(tfs, C, histXT, histXD)
            n = size(C,1);
            T = [];
            D = [];
            for i=1:n
                for j=(i+1):n
                    a = C(i,:);
                    b = C(j,:);
                    d = b - a;
                    dd = norm(d) + 0.01;
                    D = [D; dd];
                    tt = acos(d(1) ./ dd);
                    if isnan(tt)
                        tt  = 0;
                    end
                    T = [T; tt];
                end
            end
            histT = hist(T, histXT);
            histTN = histT / norm(histT);
            histTN(isnan(histTN)) = 0;
            histD = hist(D, histXD);
            histDN = histD / norm(histD);
            histDN(isnan(histDN)) = 0;
        end
        
        function im2 = filterOutTextAreasSmooth(tfs, im)
            try
                im2 = roifill(im, tfs.mask);
            catch E
                disp(tfs.fname);
                disp(E);
            end
        end
        function im2 = filterOutTextAreasSolidBlack(tfs, im)
            im2 = bsxfun(@times, im, cast(~tfs.mask,class(im)));
        end
        % Disp method
        function disp(p)
            n = length(p);
            if n>1
                for i=1:n
                    disp(p(i));
                end
                return;
            else
                fprintf('TextFeatureSet with %d text blocks.\n', length(p.text));
            end
        end % disp

    end % methods
    
end % classdef
