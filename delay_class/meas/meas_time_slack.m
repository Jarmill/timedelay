classdef meas_time_slack < meas_collection
    %MEAS_TIME_SLACK Collection of slack measures in case the time-delay
    %system has a free terminal time. The marginals of the joint occupation 
    %measure are absolutely continuous with respect to shifted copies of
    %the component measures. The nonnegative slack of each absolute
    %continuity constraint are stored in this structure.
    
    
    properties
        Tmax;
%         vars;
        lags;
%         lag_span;
%         meas;
    end
    
    methods
        function obj = meas_time_slack(delay_supp)
            %MEAS_TIME_SLACK Construct an instance of this class
            %   Detailed explanation goes here
            %fill in attributes
            obj@meas_collection(delay_supp, {'t', 'x'}); 
            
            obj.Tmax = delay_supp.Tmax;
            obj.lags = delay_supp.lags;
%             obj.vars = struct('t', delay_supp.vars.t, 'x', delay_supp.vars.x);                                    
            
            %create measures
%             obj.lag_span = obj.lag_intervals();

            obj = obj.generate_measures(delay_supp);
        end
                
        function obj = generate_measures(obj,delay_supp)
            %GENERATE_MEASURES Define and fill in component measures based
            %on support information
            
            X = delay_supp.get_X();
            
            Nlag = length(obj.lags);
            ind = 1:Nlag;
            
            obj.meas = cell(length(ind), 1);
            
            for i =1:length(ind)
                %figure out the suffix of the component measure variables
                curr_ind = ind(i);
%                 si = num2str(curr_ind);

                suffix = ['_slack', num2str(curr_ind)];
                
                %compute the time support and define the component measure                
                
                t_supp_curr = obj.vars.t*(obj.Tmax-obj.vars.t) >= 0;
                
                obj.meas{i} = obj.meas_def({'t', 'x'}, suffix, [t_supp_curr; X]);                                                                     
            end                                    
        end        
        
        function mom_out = mom_index(obj, ind_lag, d)
            %get the monomial sequence of this constraint
            mom_out = obj.meas{ind_lag}.mom_monom(d);
        end
    end
end

