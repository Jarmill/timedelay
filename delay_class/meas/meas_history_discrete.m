classdef meas_history_discrete < meas_history
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        dt;
    end
    
    methods
        function obj = meas_history_discrete(delay_supp)
            %UNTITLED Construct an instance of this class
            %   Detailed explanation goes here
            obj@meas_history(delay_supp);
            obj.dt = delay_supp.dt;
        end               
        
        function mom_out = shaping_mom_const(obj, d)
            %the histories are constant in time within their allowed set
            v = mmon([obj.vars.t; obj.vars.x], d);
            Lv = 0;
            if (isfield(obj.vars, 't') && ~isempty(obj.vars.t))
                Lv = subs(v, obj.vars.t, obj.vars.t + obj.dt) - v;
            end
            
            mom_out=0;
            for i = 1:length(obj.meas)
                Lv_curr = obj.meas{i}.var_sub(obj.get_vars(), Lv);
                
                mom_out = mom_out + mom(Lv_curr);
            end
            
        end
        
        function t_mom = hist_time_mom(obj, d, span_curr, x0)
            %HIST_TIME_MOM find the moments of the time distribution
            t_mom = TrainMom(d, span_curr, obj.dt, x0);
        end
        
        function mom_history = hist_integrate(obj, d, span_curr, X_history)
            %HIST_INTEGRATE find the moments of supplied (single) initial
            %history
            t_sample = span_curr(1):obj.dt:(span_curr(2)-obj.dt);
            x_sample = X_history(t_sample);

            dv = genPowGlopti(length(obj.vars.x)+1, d);
            
            mom_history = monom_rect_int(t_sample, x_sample, dv);  
        end
        
    end
end

