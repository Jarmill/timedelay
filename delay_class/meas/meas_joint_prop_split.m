classdef meas_joint_prop_split < meas_joint_split
    %MEAS_JOINT_PROP_SPLIT Joint occupation measures corresponding to the proportional time-delay
    %system. Measures are of the form (t, x0, x1) for times [0, s T] and
    %[sT, T] (in the case of a single delay)
    %   Detailed explanation goes here
    
    %
    %storage: lags in increasing order
    %lags: [s1, s2]: coords (0, x(t), x(s1 t), x(s2 t))
        
    methods
        %% setup
        function obj = meas_joint_prop_split(delay_supp)
            %MEAS_COMPONENT Construct an instance of this class
            
            %fill in attributes
            obj@meas_joint_split(delay_supp);             
        end
        
        function lag_span = lag_intervals(obj, lags, Tmax)
            %find the time support of each component measure
            %assume the lags are in increasing order tau1 < tau2 < tau3...
            
            %proportional-time-delay
            %[0, T tau1, T tau2, T tau3, T]
            
            if nargin < 2
                lags = obj.lags;
            end
            
            if nargin < 3
                Tmax = obj.Tmax;
            end
            
            %all time lag intervals together
            lag_all = [0, Tmax*lags, Tmax];
            
            %stack up the spans together
            lag_span = [lag_all(1:end-1); lag_all(2:end)];
        end
               
        %% getters                             
        function mom_out = mom_monom_scale(obj, ind_lag, dmin, dmax)
            %MOM_MONOM_SCALE moments of monomials where the time is scaled
            %from 't' to 't*scale' where scale<1
            
            if nargin < 5
                dmax = dmin;
                dmin = 0;
            end
            
            if ind_lag == 0
                %no delay shifting occurs
                mom_out = obj.mom_monom_marg(0, dmin, dmax);
            else
                %there is a nontrivial delay shift
                
                %TODO: check this. is it t+tau) or (t-tau)?
                lag_curr = obj.lags(ind_lag);
                var_shift = [obj.vars.t/lag_curr; obj.vars.x];
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
                for i = 1:ind_lag
                    v_curr = obj.meas{i}.var_sub(obj.get_vars(), monom_shift);
                
                    mom_out = mom_out + mom(v_curr/lag_curr);
                end
            end
        end

    end
end

