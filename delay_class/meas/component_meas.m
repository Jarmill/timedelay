classdef component_meas
    %COMPONENT_MEAS Collection of component measures for time delay
    %systems. Enforces data consistency by shifting component measures
    %   Detailed explanation goes here
    
    properties
        Tmax;
        vars;
        lags;
        lag_span;
        meas;
        
        %TODO: hybrid systems with time delays? Does it make sense to have
        %more than one location?
        
        %TODO: Treat the shaping constraints
    end
    
    methods
        function obj = component_meas(delay_supp)
            %COMPONENT_MEAS Construct an instance of this class
            
            %fill in attributes
            obj.Tmax = delay_supp.Tmax;
            obj.lags = delay_supp.lags;
            obj.vars = struct('t', delay_supp.vars.t, 'x', delay_supp.vars.x);                                    
            
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
            lag_all = [-lag_rev, 0, Tmax - lag_rev, Tmax];
            
            %stack up the spans together
            lag_span = [lag_all(1:end-1); lag_all(2:end)];
        end
        
        function obj = generate_measures(obj,delay_supp)
            %GENERATE_MEASURES Define and fill in component measures based
            %on support information
            
            X_history = delay_supp.get_X_history();
            
            Nlag = length(obj.lags);
            ind = -Nlag:Nlag;
            
            obj.meas = cell(length(ind), 1);
            
            for i =1:length(ind)
                %figure out the suffix of the component measure variables
                curr_ind = ind(i);
%                 si = num2str(curr_ind);
                if curr_ind < 0
                    suffix = ['_n', num2str(-curr_ind)];
                elseif curr_ind == 0 
                    suffix = ['_z'];
                else
                    suffix = ['_p', num2str(curr_ind)];
                end
                
                %compute the time support and define the component measure
                span = obj.lag_span(:, i); 
                
                t_supp_curr = (obj.vars.t - span(1))*(span(2) - obj.vars.t) >= 0;
                
                obj.meas{i} = obj.meas_def(suffix, [t_supp_curr; X_history]);                                                                     
            end                                    
        end

        function meas_new = meas_def(obj, suffix, supp_ref)           
            %declare a variable for each measure (index ind in the union)
            vars_new = struct('t', [], 'x', []);           
            varnames = fields(vars_new);
            for i = 1:length(varnames)
                curr_name = varnames{i};
                curr_var = obj.vars.(curr_name);
                
                if ~isempty(curr_var)
                    %declare a new variable
                    new_name = [curr_name, suffix];
                    mpol(new_name, length(curr_var), 1);
                    %load the new variable into vars_new
                    vars_new.(curr_name) = eval(new_name);
                end
            end
            
           supp_new = subs_vars(supp_ref, [obj.vars.t; obj.vars.x], ...
                                [vars_new.t; vars_new.x]);
           
            
            %define the measure
            %TODO: extend to uncertain systems
            %or rather, decide if there is dead time and the control/input
            %w should be stored in the component measures
            meas_new = meas_base(vars_new, supp_new);
        end
        
        function cons = shaping_cons(obj, d)
            %TODO: Treat the shaping constraints
            cons = [];
        end
        
        function shift_out = mom_shift(obj, ind_lag, dmin, dmax)
            %MOM_SHIFT moments of shifted copies of component measures for
            %the data at time lag 'ind_lag'

            if nargin < 4
                dmax = dmin;
                dmin = 0;
            end
            
            
            
            %identify the current lag
            Nlag = length(obj.lags);
            if ind_lag == 0 
                lag_curr = 0;
            else
                lag_curr = obj.lags(ind_lag);
            end
            
            %iterate through component measures     
            shift_out = 0;
            for i = 0:Nlag
                ind_curr = Nlag + i - ind_lag + 1;
                
                mom_shift_curr = obj.meas{ind_curr}.mom_monom_shift(...
                        -lag_curr, dmin, dmax);
                    
                shift_out = shift_out + mom_shift_curr;
            end            
        end
    end
end

