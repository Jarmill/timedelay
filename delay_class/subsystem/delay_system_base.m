classdef delay_system_base < delay_system_interface
    %DELAY_SYSTEM_BASE A subsystem for time-delay system containing the
    %joint occupation measure in the absence of uncertainty or input
    
    properties
        varnames = {'t', 'x', 'x_lag'};        
    end
    
    methods
        function obj = delay_system_base(delay_supp, f)
            %DELAY_SYSTEM_BASE Construct an instance of this class
            %   Detailed explanation goes here

            
            obj@delay_system_interface(delay_supp, f, [], @meas_joint_base);
        end
        
        function Ay = cons_liou(obj,d)
            %CONS_LIOU Liouville Equation contribution of Lie derivative
            
            Ay = obj.meas_occ.mom_lie(d, obj.vars, obj.f);
        end
        
        function mom_out = mom_marg(obj, ind_lag, d)
            %MOM_MARG return the marginals of the joint occupation measure
            %for use in data consistency constraints
            
            mom_out = obj.meas_occ.mom_monom_marg(ind_lag, d);
            
        end
        
        function obj = dual_process(v, phi)
            %DUAL_PROCESS store dual functions and compute nonnegative
            %functions for this subsystem
        end
    end
end

