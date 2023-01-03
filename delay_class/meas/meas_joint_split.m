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
        
        function mom_out = mom_lie(obj, d, vars_old, f_old)
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
        
        
        
    end
end

