classdef meas_history < meas_collection
    %MEAS_HISTORY A collection of measures that define the history of
    %trajectories in the estimation and control programs
    
    properties
        Tmax;       %terminal time

        lags;       %time delays in system
        lag_span;   %spans of component measures
        
        Npts_history = 1000;
    end
    
    methods
        function obj = meas_history(delay_supp)
            %MEAS_HISTORY Construct an instance of this class
            %   Detailed explanation goes here
            obj@meas_collection(delay_supp, {'t', 'x'}); 
            
            obj.Tmax = delay_supp.Tmax;
            obj.lags = delay_supp.lags;
%             obj.vars = struct('t', delay_supp.vars.t, 'x', delay_supp.vars.x);                                    
            
            %create measures
            obj.lag_span = obj.lag_intervals();

            obj = obj.generate_measures(delay_supp);
        end
        
        function lag_span = lag_intervals(obj, lags)
            %find the time support of each component measure
            %assume the lags are in increasing order tau1 < tau2 < tau3...
            
            if nargin < 2
                lags = obj.lags;
            end           
            
            lag_rev = reshape(lags(end:-1:1), 1, []);
            
            %all time lag intervals together
            %focus on the histories
            lag_all = [-lag_rev, 0];
            
            %stack up the spans together
            lag_span = [lag_all(1:end-1); lag_all(2:end)];
        end
    
    
        function obj = generate_measures(obj,delay_supp)
            %GENERATE_MEASURES Define and fill in history measures based
            %on support information
            %
            %index ordering: {[-tau1, 0], [-tau2, -tau1], [-tau3, -tau2]}

%             X = delay_supp.get_X();
            X_history = delay_supp.get_X_history_supp();

            Nlag = length(obj.lags);
%             ind = -Nlag:Nlag;

            obj.meas = cell(Nlag, 1);

            for i =1:Nlag
                %figure out the suffix of the history measure variables
                span = obj.lag_span(:, i); 
                t_supp_curr = (obj.vars.t - span(1))*(span(2) - obj.vars.t) >= 0;   

                suffix = ['_h', num2str(Nlag-i+1)];
                curr_supp = t_supp_curr; %assume single trajectory with moment substitution
                if isa(X_history, 'supcon')
                    curr_supp = [t_supp_curr; X_history]; %multiple histories
                end

                %compute the time support and define the component measure
                obj.meas{Nlag-i+1} = obj.meas_def({'t', 'x'}, suffix, curr_supp);                                                                     
            end                                    
        end
        
        %% moment computation
        function mom_out = mom_shift(obj, ind_lag, dmin, dmax)
            %MOM_SHIFT moments of shifted copies of component measures for
            %the data at time lag 'ind_lag'. This is used to constrain the
            %marginals of the joint occupation measure

            if nargin < 4
                dmax = dmin;
                dmin = 0;
            end
                   
            mom_out = 0;
            if ind_lag >= 1
                %identify the current lag
                Nlag = length(obj.lags);
                ind_lag = min(ind_lag, Nlag);

                %perform the shift
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

                %iterate through component history measures     
                for i = 1:(ind_lag)
                    v_curr = obj.meas{i}.var_sub(obj.stack_vars(obj.vars), monom_shift);

                    mom_out = mom_out + mom(v_curr);
                end      
            end
        end
        
        %% pin down moment costraints based on the history
        function mom_con = history_traj_con(obj, d, X_history)
            %MOMS: pin down moments of history measures according to 
            %*single* trajectory in X_history            
            %
            %
             
            mom_con = [];
            if ~isa(X_history, 'supcon')
                %history is either a trajectory (function handle) or a
                %constant trajectory (double)                                                        
                for i = 1:length(obj.lags)
                    
                    span_curr = obj.lag_span(:, length(obj.lags)-i+1);
                    if isa(X_history, 'function_handle')
                        %numerically integrate trajectories
                        mom_history = obj.hist_integrate(d, span_curr, X_history);
                        
%                         t_sample = linspace(span_curr(1), span_curr(2), obj.Npts_history);
%                         x_sample = X_history(t_sample);
%                                                 
%                         
%                         dv = genPowGlopti(length(obj.vars.x)+1, d);
%                         mom_history = monom_int(t_sample, x_sample, dv);                        
                    else
                        %X_history is a double
                        %constant history trajectory
                        mom_history= obj.hist_time_mom(d, span_curr, X_history);
                    end
                    
                    %time marginals of history component measure
                    moms = obj.meas{i}.mom_monom(d);

                    mom_con = [mom_con; moms == mom_history];
                end
            end
            
        end
        
        function cons = shaping_cons(obj, d)
            %TODO: Treat the history shaping constraints
            cons = [];
            %likely restrict to constants only
        end
        
        function mom_out = shaping_mom_const(obj, d)
            %the histories are constant in time within their allowed set
            v = mmon([obj.vars.t; obj.vars.x], d);
            Lv = 0;
            if (isfield(obj.vars, 't') && ~isempty(obj.vars.t))
                Lv = Lv + diff(v, obj.vars.t);
            end
            
            mom_out=0;
            for i = 1:length(obj.meas)
                Lv_curr = obj.meas{i}.var_sub(obj.get_vars(), Lv);
                
                mom_out = mom_out + mom(Lv_curr);
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
                    span_curr = obj.lag_span(:, length(obj.lags)-i+1);
                    %moments of lebesgue distribution on time support
%                     t_mom = ConstMom(d, span_curr, []);
                    t_mom = obj.hist_time_mom(d, span_curr, []);
                    %time marginals of history component measure
                    moms = mom(obj.meas{i}.vars.t.^genPowGlopti(1, d));

                    mom_con = [mom_con; moms == t_mom];
                end
            end
            
        end
        
        function t_mom = hist_time_mom(obj, d, span_curr, x0)
            %HIST_TIME_MOM find the moments of the time distribution
            t_mom = ConstMom(d, span_curr, x0);
        end
        
    
        function mom_history = hist_integrate(obj, d, span_curr, X_history)
            %HIST_INTEGRATE find the moments of supplied (single) initial
            %history
            t_sample = linspace(span_curr(1), span_curr(2), obj.Npts_history);
            x_sample = X_history(t_sample);


            dv = genPowGlopti(length(obj.vars.x)+1, d);
            mom_history = monom_int(t_sample, x_sample, dv);  
        end
        
    end
end

