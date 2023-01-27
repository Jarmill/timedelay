classdef meas_joint_split < meas_collection
    %MEAS_SPLIT Joint occupation measures corresponding to the time-delay
    %system. Measures are of the form (t, x0, x1) for times [0, T-tau] and
    %[T-tau, T] (in the case of a s ingle delay
    %   Detailed explanation goes here
    
    properties
        Tmax;       %terminal time
%         vars;
        lags;       %time delays in system
        lag_span;   %spans of component measures
        
    end
    
    methods
        %% setup
        function obj = meas_joint_split(delay_supp)
            %MEAS_COMPONENT Construct an instance of this class
            
            %fill in attributes
            obj@meas_collection(delay_supp, {'t', 'x', 'x_lag'}); 
            
            obj.Tmax = delay_supp.Tmax;
            obj.lags = delay_supp.lags;
%             obj.vars = struct('t', delay_supp.vars.t, 'x', delay_supp.vars.x);                                    
            
            %create measures
            obj.lag_span = obj.lag_intervals();

            obj = obj.generate_measures(delay_supp);
        end
        
        function lag_span = lag_intervals(obj, lags, Tmax)
            %find the time support of each component measure
            %assume the lags are in increasing order tau1 < tau2 < tau3...
            
            if nargin < 2
                lags = obj.lags;
            end
            
            if nargin < 3
                Tmax = obj.Tmax;
            end
            
            lag_rev = reshape(lags(end:-1:1), 1, []);
            
            %all time lag intervals together
            lag_all = [0, Tmax - lag_rev, Tmax];
            
            %stack up the spans together
            lag_span = [lag_all(1:end-1); lag_all(2:end)];
        end
        
        function obj = generate_measures(obj,delay_supp)
            %GENERATE_MEASURES Define and fill in Joint Occupation measures
            %based on the 
            
%             X = delay_supp.get_X();
%             X_history = delay_supp.get_X_history_supp();
            
            X_joint = delay_supp.get_X_joint();
            Nlag = length(obj.lags);
           
            
            obj.meas = cell(Nlag+1, 1);
            
            for i =1:Nlag+1
                %figure out the suffix of the component measure variables
%                 curr_ind = ind(i);
%                 si = num2str(curr_ind);
                span = obj.lag_span(:, i); 
                t_supp_curr = (obj.vars.t - span(1))*(span(2) - obj.vars.t) >= 0;   
                
%                 if curr_ind == 0 
%                     suffix = ['_z'];
%                     curr_supp = [t_supp_curr;X];
%                 else
                    suffix = ['_d', num2str(i-1)];
                    curr_supp = [t_supp_curr;X_joint];
%                 end
                
                %compute the time support and define the component measure

                                                
                obj.meas{i} = obj.meas_def({'t', 'x', 'x_lag'}, suffix, curr_supp);                                                                     
            end                                    
        end
        
        %% getters        

        
        function mom_out = mom_monom_marg(obj, ind_lag, dmin, dmax)
            %MOM_MONOM_MARG get moments of the marginals of the joint 
            %occupation measure. used for data consistency constraints
            
            %(t, x(t+ lag)) marginal
            
            if nargin < 4
                dmax = dmin;
                dmin = 0;
            end                        
            
            %perform the substitution
            t_curr = obj.vars.t;
            if ind_lag == 0                
                x_curr = obj.vars.x;
            else
                x_curr = obj.vars.x_lag(:, ind_lag);
            end
            v = mmon([t_curr; x_curr], dmin, dmax);
            vars_curr = obj.get_vars();
            
            %sum over the measures
            mom_out = 0;
            for i = 1:length(obj.meas)
                v_curr = obj.meas{i}.var_sub(vars_curr, v);
                
                mom_out = mom_out + mom(v_curr);
            end
        end
        
        function mom_out = mom_monom_shift(obj, ind_lag, dmin, dmax)
            %MOM_MONOM_SHIFT moments of monomials where the time is shifted
            %from 't' to 't-time_shift'
            
            if nargin < 5
                dmax = dmin;
                dmin = 0;
            end
            
            if ind_lag == 0
                %no delay shifting occurs
                mom_out = obj.mom_monom_marg(0, dmin, dmax);
            else
                %there is a nontrivial delay shift
                
                %TODO: check this. is it (t+tau) or (t-tau)?
                var_shift = [obj.vars.t+obj.lags(ind_lag); obj.vars.x];
                nvar = length(var_shift);
                
                %monomial generation copied over from @mpol/mmon
                vpow = [];
                for k = dmin:dmax
                    vpow = [vpow;genpow(nvar,k)];
                end
                %powers of the shifted times
                monom_shift = prod(var_shift'.^vpow, 2);              

                %return the output
                mom_out = 0;
                for i = 1:(length(obj.lags) - ind_lag+1)
                    v_curr = obj.meas{i}.var_sub(obj.get_vars(), monom_shift);
                
                    mom_out = mom_out + mom(v_curr);
                end
            end
        end
        
        
        %% generator expressions
        
        function mom_out = mom_lie(obj, d, vars_old, f_old)
            %Lie derivative v -> v_t + f' * grad_x v
            mom_out = 0;
            v = mmon([obj.vars.t; obj.vars.x], d);
            Lv = diff(v, obj.vars.x)*f_old;
            if (isfield(obj.vars, 't') && ~isempty(obj.vars.t))
                Lv = Lv + diff(v, obj.vars.t);
            end
            
            for i = 1:length(obj.meas)
                Lv_curr = obj.meas{i}.var_sub(vars_old, Lv);
                
                mom_out = mom_out + mom(Lv_curr);
            end
        end
        
        function mom_out = mom_hess(obj, d, vars_old, g_old)
            %hessian v -> g'*hess_x*g
            mom_out = 0;
            v = mmon([obj.vars.t; obj.vars.x], d);
%             v_hess = v;
            for k = 1:length(v)
                v_partial = diff(v(k), obj.vars.x);
                hess_curr = diff(v_partial', obj.vars.x);
                v_hess(k) = g_old'* hess_curr *g_old;
            end
            
            for i = 1:length(obj.meas)
                hess_curr = obj.meas{i}.var_sub(vars_old, v_hess);
                
                mom_out = mom_out + mom(hess_curr);
            end
        end
        
        function mom_out = mom_push(obj, d, vars_old, f_old, tshift)
            
            if nargin < 5
                tshift = 0;
            end
            mom_out = 0;
            v = mmon([obj.vars.t; obj.vars.x], d);
%             v_hess = v;
            push = subs(v, [obj.vars.t; obj.vars.x], [obj.vars.t+tshift; f_old]);
%             for k = 1:length(v)
%                 push
%             end
            
            for i = 1:length(obj.meas)
%                 hess_curr = obj.meas{i}.var_sub(vars_old, v_hess);
                
                mom_out = mom_out + mom(push);
            end
        end
        
        
        
    end
end

