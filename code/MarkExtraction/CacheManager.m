classdef CacheManager < handle
    
    properties (Constant)
        LOCAL_CACHE = 1;
        SHARED_CACHE = 2;
    end
    
    properties
        LocalCachePath   % Local cache path
        SharedCachePath % Cache path in dropbox
        cacheprefix % (Optional) prefix specifying the corpus we're using
    
        BFILTER_WIDTH
        BFILTER_DOMAINSIGMA
        BFILTER_RANGESIGMA
    end
    
    methods
        
        function obj = CacheManager(BFILTER_WIDTH, BFILTER_DOMAINSIGMA, BFILTER_RANGESIGMA)
            obj.BFILTER_WIDTH = BFILTER_WIDTH;
            obj.BFILTER_DOMAINSIGMA = BFILTER_DOMAINSIGMA;
            obj.BFILTER_RANGESIGMA = BFILTER_RANGESIGMA;
        end
        
        function value = get.LocalCachePath(obj)
            value = obj.LocalCachePath;
        end
        function set.LocalCachePath(obj,val)
            if ~exist(val, 'dir')
                error('Specified path is not a directory');
            else
                obj.LocalCachePath = val;
            end
        end
        function value = get.SharedCachePath(obj)
            value = obj.SharedCachePath;
        end
        function set.SharedCachePath(obj,val)
            if ~exist(val, 'dir')
                error('Specified path is not a directory');
            else
                obj.SharedCachePath = val;
            end
        end

        function fullpath = getLoadPath(obj,filename)
            loadFileName = strcat(obj.cacheprefix, filename);
            
            fullpath='';
            
            % Check local cache
            if exist(fullfile(obj.LocalCachePath, loadFileName),'file')
                fprintf(2, 'In local cache\n');
                fullpath = fullfile(obj.LocalCachePath, loadFileName);
                %load(fullfile(obj.LocalCachePath, loadFileName));
            % Check shared cache
            elseif exist(fullfile(obj.SharedCachePath, loadFileName),'file')
                fprintf(2, 'In shared cache\n');
                fullpath = fullfile(obj.SharedCachePath, loadFileName);
                %fprintf(2, strcat('Shared cache: ', loadMessage));
                %load(fullfile(obj.SharedCachePath, loadFileName));
            end
        end
        
        % Saves the elements into
        function fullpath = getSavePath(obj,filename,cache)
            saveFileName = strcat(obj.cacheprefix,filename);
            if cache == obj.LOCAL_CACHE
                fullpath = fullfile(obj.LocalCachePath,saveFileName);
            elseif cache == obj.SHARED_CACHE
                fullpath = fullfile(obj.SharedCachePath,saveFileName);
            else
                error('Cache location not specified');
            end
        end
    end
end