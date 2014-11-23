classdef imageProcessor < handle

    properties
        impath
        curFile
        curFile_noext
        im
        filteredIm
        grayIm
        cache
        metadata = struct;
        COLOR
    end
    
    methods
        %% Constructor
        function obj = imageProcessor(curFile, fullpath, cache)
            extInd = strfind(curFile, '.');
            
            obj.cache = cache;
            
            obj.curFile = curFile;
            obj.curFile_noext = obj.curFile(1:extInd-1);
            obj.impath = fullfile(fullpath, obj.curFile);
        end

        function out = readImages(obj)
            fprintf(2, strcat(obj.impath,'\n'));
            [im grayIm] = readimage(obj.impath);
            obj.im = im;
            obj.grayIm = grayIm;
            
            out = obj;
        end
        
        function obj = filterImage(obj)
            BFILTER_WIDTH = obj.cache.BFILTER_WIDTH;
            BFILTER_DOMAINSIGMA = obj.cache.BFILTER_DOMAINSIGMA;
            BFILTER_RANGESIGMA = obj.cache.BFILTER_RANGESIGMA;            
            
            % Check the cache
            filteredIm_file = strcat(obj.curFile_noext,'_filtered_',sprintf('%d', BFILTER_WIDTH), '_', ...
                 sprintf('%2.2f',BFILTER_DOMAINSIGMA),'_',sprintf('%2.2f',BFILTER_RANGESIGMA),'.mat');

            % Check the cache
            loadPath = obj.cache.getLoadPath(filteredIm_file);            
            % In cache
            if ~isempty(loadPath)
                load(loadPath);
            % Not in cache: filter the image
            else
                fprintf(2,'Applying bilateral filter...\n');
                bfiltert = tic;
                im = bfilter2(obj.im, BFILTER_WIDTH, [BFILTER_DOMAINSIGMA BFILTER_RANGESIGMA]);
                bfilter_time = toc(bfiltert);
                
                fprintf(2, 'Saving fitered image...\n');
                save(obj.cache.getSavePath(filteredIm_file, obj.cache.LOCAL_CACHE), ...
                    'im','BFILTER_WIDTH','BFILTER_DOMAINSIGMA','BFILTER_RANGESIGMA','bfilter_time');
            end
            
            obj.filteredIm = im;
            obj.metadata.('BFILTER_WIDTH') = BFILTER_WIDTH;
            obj.metadata.('BFILTER_DOMAINSIGMA') = BFILTER_DOMAINSIGMA;
            obj.metadata.('BFILTER_RANGESIGMA') = BFILTER_RANGESIGMA;
            obj.metadata.('bfilter_time') = bfilter_time;
        end
        
        %% Getters and setters
        function value = get.im(obj)
            value = obj.im;
        end
        function obj = set.im(obj, value)
            COLOR = size(value,3) > 1;
            if ~COLOR
                value = repmat(value, [1 1 3]);
            end

            obj.im = value;
        end
        
        function value = get.grayIm(obj)
            value = obj.grayIm;
        end
%         function obj = set.grayIm(obj, value)
%             obj.grayIm = value;
%         end
        
        function value = get.filteredIm(obj)
            value = obj.filteredIm;
        end
%         function obj = set.filteredIm(obj, value)
%             obj.filteredIm = value;
%         end
    end

end