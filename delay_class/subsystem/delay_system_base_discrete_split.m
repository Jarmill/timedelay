classdef delay_system_base_discrete_split < delay_system_base_split
    %DELAY_SYSTEM_BASE_DISCRETE_SPLIT A subsystem for a discrete-time 
    %time-delay system containing the split-joint occupation measure in the
    %absence of uncertainty or input
        
    properties
        dt;
    end
    
    methods
        function obj = delay_system_base_discrete_split(delay_supp, f)
            %UNTITLED5 Construct an instance of this class
            %   Detailed explanation goes here
            obj@delay_system_base_split(delay_supp, f);
            
            obj.dt = delay_supp.dt;
        end
        
        function Ay = cons_liou(obj,d)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            %
            
            %use the internals of meas_occ.mom_push            
            y_push = obj.meas_occ.mom_push(d, obj.get_vars(), obj.f, obj.dt);
            y_orig = obj.meas_occ.mom_monom_marg(0, d);
            
            Ay = y_push - y_orig;
%             v_occ = obj.meas_occ.monom_proj(d);
%             push_occ = obj.meas_occ.var_sub(obj.get_vars(), push_f);
%             Rv = subs(v_occ, , ...
%                 push_occ);
%             mom_out = mom(Rv);
        end
    end
end

