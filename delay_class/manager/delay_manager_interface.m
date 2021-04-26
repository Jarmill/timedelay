classdef delay_manager_interface < manager_interface
    %DELAY_MANAGER_INTERFACE A generic class to pose LMIs in measures for
    %time-delay systems.    
    
    properties
        Property1
    end
    
    methods
        function obj = delay_manager_interface(loc_supp)
            %DELAY_MANAGER_INTERFACE Construct an instance of this class
            %   Detailed explanation goes here
            obj.Property1 = inputArg1 + inputArg2;
        end
        
        function outputArg = method1(obj,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
        end
    end
end

