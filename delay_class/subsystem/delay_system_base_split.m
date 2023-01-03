classdef delay_system_base_split < delay_system_interface
    %DELAY_SYSTEM_SPLIT A subsystem for time-delay system containing the
    %split-joint occupation measure in the absence of uncertainty or input
    
    properties
        varnames = {'t', 'x', 'x_lag'};        
    end    
    
    methods
        %% setup
        function obj = delay_system_base_split(delay_supp, f)
            %DELAY_SYSTEM_BASE Construct an instance of this class
            %   Detailed explanation goes here

            
            obj@delay_system_interface(delay_supp, f, [], @meas_joint_split);
            
            %correctly define the occupation measure
            obj.meas_occ = obj.meas_def_split(delay_supp);
        end
        
        function MO = meas_def(obj, vars, suffix, supp_new)
            %a technical kludge to define the split occupation measure
            MO = [];
        end
        
        function MO = meas_def_split(obj, delay_supp)
            %the actual definition of the occupation measure that matters
            MO = meas_joint_split(delay_supp);
        end
        
        %% get moments
        function Ay = cons_liou(obj,d)
            %CONS_LIOU Liouville Equation contribution of Lie derivative
            
            Ay = obj.meas_occ.mom_lie(d, obj.get_vars(), obj.f);
        end
        
        function mom_out = mom_marg(obj, ind_lag, d)
            %MOM_MARG return the marginals of the joint occupation measure
            %for use in data consistency constraints
            
            mom_out = obj.meas_occ.mom_monom_marg(ind_lag, d);
            
        end
        
        
        %% bugged duals for future work
        function obj = dual_process(v, phi)
            %DUAL_PROCESS store dual functions and compute nonnegative
            %functions for this subsystem
        end
    end
end

