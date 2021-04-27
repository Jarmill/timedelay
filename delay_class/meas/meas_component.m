classdef meas_component < meas_collection
    %MEAS_COMPONENT Collection of component measures for time delay
    %systems. Enforces data consistency by shifting component measures
    %   Detailed explanation goes here
    
    properties
        Tmax;       %terminal time
%         vars;
        lags;       %time delays in system
        lag_span;   %spans of component measures
%         meas;
        
        %TODO: hybrid systems with time delays? Does it make sense to have
        %more than one location?
        
        %TODO: Treat the history shaping constraints
    end
    
    methods
        function obj = meas_component(delay_supp)
            %MEAS_COMPONENT Construct an instance of this class
            
            %fill in attributes
            obj@meas_collection(delay_supp, {'t', 'x'}); 
            
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
            lag_all = [-lag_rev, 0, Tmax - lag_rev, Tmax];
            
            %stack up the spans together
            lag_span = [lag_all(1:end-1); lag_all(2:end)];
        end
        
        function obj = generate_measures(obj,delay_supp)
            %GENERATE_MEASURES Define and fill in component measures based
            %on support information
            
            X = delay_supp.get_X();
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
                    curr_supp = X_history;
                elseif curr_ind == 0 
                    suffix = ['_z'];
                    curr_supp = X;
                else
                    suffix = ['_p', num2str(curr_ind)];
                    curr_supp = X;
                end
                
                %compute the time support and define the component measure
                span = obj.lag_span(:, i); 
                
                t_supp_curr = (obj.vars.t - span(1))*(span(2) - obj.vars.t) >= 0;
                
                obj.meas{i} = obj.meas_def({'t', 'x'}, suffix, [t_supp_curr; curr_supp]);                                                                     
            end                                    
        end
        
        function mmmon_out = mom_monom_shift(obj, meas_in, time_shift, dmin, dmax)
            %MOM_MONOM_SHIFT moments of monomials where the time is shifted
            %from 't' to 't-time_shift'
            if nargin < 5
                dmax = dmin;
                dmin = 0;
            end
            
%             monom = obj.monom(dmin, dmax);
%             monom_shift = subs(monom, obj.vars.t, obj.vars.t - time_shift);
%           @mpol/subs takes a long time, mostly due to the 'locate' call
            
            if time_shift == 0
                mmmon_out = meas_in.mom_monom(dmin,dmax);
            else
                var_shift = [meas_in.vars.t-time_shift; meas_in.vars.x];
                nvar = length(var_shift);
                
                %monomial generation copied over from @mpol/mmon
                vpow = [];
                for k = dmin:dmax
                    vpow = [vpow;genpow(nvar,k)];
                end

                monom_shift = prod(var_shift'.^vpow, 2);              

                mmmon_out = mom(monom_shift);
            end
        end
        
        function cons = shaping_cons(obj, d)
            %TODO: Treat the history shaping constraints
            cons = [];
        end
        
        function shift_out = mom_shift(obj, ind_lag, dmin, dmax)
            %MOM_SHIFT moments of shifted copies of component measures for
            %the data at time lag 'ind_lag'. This is used to constrain the
            %marginals of the joint occupation measure

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
                
                mom_shift_curr = obj.mom_monom_shift(obj.meas{ind_curr},...
                        -lag_curr, dmin, dmax);
                    
                shift_out = shift_out + mom_shift_curr;
            end            
        end
    end
end

