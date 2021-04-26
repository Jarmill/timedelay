classdef meas_joint_base < meas_interface
    %MEAS_JOINT_BASE Joint occupation measure in a time-delay system
    %   Only includes time and (delayed) states
    
%     properties
%         Property1
%     end
    
    methods
        function obj = meas_joint_base(vars, supp)
            %MEAS_JOINT_BASE Construct an instance of this class
            %   Detailed explanation goes here
%             obj.Property1 = inputArg1 + inputArg2;
            obj@meas_interface(vars, supp);
        end
        
        %% getters        
        function vars_out = get_vars(obj)
            %GET_VARS get variables in measure
            varnames = fields(obj.vars);
            vars_out = [];
            for i = 1:length(varnames)
                curr_var = varnames{i};
                vars_out = [vars_out; reshape(obj.vars.(curr_var), [], 1)];
            end
        end
        
        function mom_out = mom_monom_marg(obj, ind_lag, dmin, dmax)
            %MOM_MONOM_MARG get moments of the marginals of the joint 
            %occupation measure. used for data consistency constraints
            
            %(t, x(t+ lag)) marginal
            
            if nargin < 4
                dmax = dmin;
                dmin = 0;
            end                        
            
            t_curr = obj.vars.t;
            if ind_lag == 0                
                x_curr = obj.vars.x;
            else
                x_curr = obj.vars.x_lag(:, ind_lag);
            end
            
            mom_out = mom(mmon([t_curr; x_curr], dmin, dmax));
        end
        
        %% liouville moments 
        function mom_out = mom_lie(obj, d, vars_old, f_old)
            %MOM_LIE lie derivative moments
            v = mmon([obj.vars.t; obj.vars.x], d);
            f_curr = obj.var_sub(obj.stack_vars(vars_old), f_old);
            mom_out = mom(diff(v, obj.vars.x)*f_curr);
            
            if (isfield(obj.vars, 't') && ~isempty(obj.vars.t))
                mom_out = mom(diff(v, obj.vars.t)) + mom_out;
            end
        end
        
        function mom_out = mom_push(obj, d, vars_old, f_old)
            %MOM_PUSH pushforward moments v(f(x)) - v(x)
            v = mmon([obj.vars.t; obj.vars.x], d);
            f_curr = obj.var_sub(obj.stack_vars(vars_old), f_old);
            
            %TODO: let dt change as T scales
            %TODO: Also find out if digital systems work properly
            dt = 1;
            Rv = subs(v, [obj.vars.t + dt; obj.vars.x], ...
                [obj.vars.t; f_curr]);
            mom_out = mom(Rv);
        end
        
        
        function vars_out = stack_vars(obj, vars)
            %generate a stack of variables for mom_lie
            %inefficient code, fix later
            varnames = fields(vars);
            vars_out = [];
            for i = 1:length(varnames)
                curr_var = varnames{i};
                vars_out = [vars_out; reshape(vars.(curr_var), [], 1)];
            end
        end
        
        %% moment recovery (not needed)
        function [optimal, mom_out, corner] = recover(obj, tol) 
            
            %this class will always be an occupation measure, so recovery
            %is unnecesary
            if nargin < 2
                tol = 5e-4;
            end
            optimal = 0;
            mom_out = [];
            corner = obj.mmat_corner();
        end
    end
end

