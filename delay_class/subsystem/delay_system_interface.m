classdef delay_system_interface < subsystem_interface
    %DELAY_SYSTEM_INTERFACE A container for time delay subsystem methods
    
    
    properties
        
    end
    
    methods
        function obj = delay_system_interface(delay_supp, f, sys_id, meas_handle)
            %DELAY_SYSTEM_INTERFACE Construct an instance of this class
            %   Detailed explanation goes here

            obj@subsystem_interface(delay_supp, f, sys_id, [], meas_handle);


        end
        
        function cons = abscont_box(obj,d)
            %ABSCONT_BOX Summary of this method goes here
            %   Detailed explanation goes here
            cons = [];
        end
    end
end

