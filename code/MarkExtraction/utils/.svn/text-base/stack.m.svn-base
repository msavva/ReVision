% Implement a simple stack. Inherits from handle so object properties
% can be changed via methods.
classdef stack < handle
% Public properties    
    properties
        % Data is stored as a row vector
        data = [];
    end
    
% Methods
    methods
        function obj = stack(a)
            if(nargin > 0)
                % Store data as a row vector
                obj.data = a(:).';
            end
        end
        
        function obj = push(obj, val)
            obj.data = cat(2,obj.data,val);
        end
       
        function val = pop(obj)
            val = obj.data(end);
            obj.data(end) = [];
        end
        
        function l = length(obj)
            l = length(obj.data);
        end

        function e = isempty(obj)
            e = isempty(obj.data);
        end
        
        function clear(obj)
            obj.data = [];
        end
    end
end