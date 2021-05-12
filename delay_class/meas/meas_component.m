classdef meas_component < meas_collection
    %MEAS_COMPONENT Collection of component measures for time delay
    %systems. Enforces data consistency by shifting component measures
    %   Detailed explanation goes here
    
    properties
        Tmax;       %terminal time
%         vars;
        lags;       %time delays in system
        lag_span;   %spans of component measures
        
        Npts_history = 1000;
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
            X_history = delay_supp.get_X_history_supp();
            
            Nlag = length(obj.lags);
            ind = -Nlag:Nlag;
            
            obj.meas = cell(length(ind), 1);
            
            for i =1:length(ind)
                %figure out the suffix of the component measure variables
                curr_ind = ind(i);
%                 si = num2str(curr_ind);
                span = obj.lag_span(:, i); 
                t_supp_curr = (obj.vars.t - span(1))*(span(2) - obj.vars.t) >= 0;   
                
                if curr_ind < 0
                    suffix = ['_n', num2str(-curr_ind)];
                    curr_supp = [t_supp_curr]; %assume single trajectory with moment substitution
                    if isa(X_history, 'supcon')
                        curr_supp = [t_supp_curr; X_history]; %multiple histories
                    end
                elseif curr_ind == 0 
                    suffix = ['_z'];
                    curr_supp = [t_supp_curr;X];
                else
                    suffix = ['_p', num2str(curr_ind)];
                    curr_supp = [t_supp_curr;X];
                end
                
                %compute the time support and define the component measure

                                                
                obj.meas{i} = obj.meas_def({'t', 'x'}, suffix, curr_supp);                                                                     
            end                                    
        end
                        
        
        function mom_con = history_traj_con(obj, d, X_history)
            %MOMS: pin down moments of history measures according to single
            %trajectory in X_history            
             
            mom_con = [];
            if ~isa(X_history, 'supcon')
                %history is either a trajectory (function handle) or a
                %constant trajectory (double)                                                        
                for i = 1:length(obj.lags)
                    
                    span_curr = obj.lag_span(:, i);
                    if isa(X_history, 'function_handle')
                        %numerically integrate trajectories
                        t_sample = linspace(span_curr(1), span_curr(2), obj.Npts_history);
                        x_sample = X_history(t_sample);
                                                
                        
                        dv = genPowGlopti(length(obj.vars.x)+1, d);
                        mom_history = monom_int(t_sample, x_sample, dv);                        
                    else
                        %X_history is a double
                        %constant history trajectory
                        mom_history= ConstMom(d, obj.lag_span(:, i), X_history);
                    end
                    
                    %time marginals of history component measure
                    moms = obj.meas{i}.mom_monom(d);

                    mom_con = [mom_con; moms == mom_history];
                end
            end
            
        end
        
        function mom_con = history_free_con(obj, d, X_history)
            %MOM_HISTORY pin down t-marginal moments of history measures to
            %the lebesgue measure in time on the appropriate interval            
             
            mom_con = [];
            if strcmp(class(X_history), 'supcon')
                %only perform this bounding histories if there are multiple
                %histories allowed. This occurs when X_history is a support
                %set.
                for i = 1:length(obj.lags)

                    %moments of lebesgue distribution on time support
                    t_mom = ConstMom(d, obj.lag_span(:, i), []);

                    %time marginals of history component measure
                    moms = mom(obj.meas{i}.vars.t.^genPowGlopti(1, d));

                    mom_con = [mom_con; moms == t_mom];
                end
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
