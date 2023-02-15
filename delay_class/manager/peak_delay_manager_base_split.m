classdef peak_delay_manager_base_split < manager_interface
    %PEAK_DELAY_MANAGER perform peak estimation single time-delay trajectory 
    %through the use of an LMI relaxation 
    %
    %this involves the joint-split occupation measure formulation
    

    
    methods
        function obj = peak_delay_manager_base_split(loc_supp,f, objective)
            %WEAK_MANAGER Construct an instance of this class
            %   Detailed explanation goes here
%             obj.Property1 = inputArg1 + inputArg2;
            
            %weak solution: fixed terminal time
            loc_supp.FREE_TERM = 1;
%             loc_supp.SCALE_TIME = 0; %change this later
            if loc_supp.PROPORTIONAL
                loc_curr = delay_location_prop_split(loc_supp, f, objective);            
            else
                loc_curr = delay_location_base_split(loc_supp, f, objective);            
            end
            
            obj@manager_interface(loc_curr);
        end
        
        function [objective, mom_con, supp_con, len_dual] = cons(obj,d, Tmax)
            %PEAK_CONS formulate support and measure constraints for peak
            %program at degree d
            %Input:
            %   d:      Monomials involved in relaxation (2*order)
            %   Tmax:   Maximum time (only when time-independent)
            %
            %Output:
            %   objective:  target to maximize  (@mom)
            %   mom_con:    moment constraints  (@momcon)
            %   supp_con:   support constraints (@supcon)
          

            supp_con = obj.loc.supp_con();       %support constraint                 
                        
            %gather all constraints in each location
            %(loop over locations)
            [objective, cons_eq, cons_ineq, len_dual] = ...
                obj.loc.all_cons(d);
            %finalize moment constraints
            
            %mass of initial measure sums to one
            mass_init_con = (obj.loc.mass_init() - 1 == 0);            

            mom_con = [cons_eq; cons_ineq; mass_init_con];            
        end          
        
        function obj = dual_process(obj, d, dual_rec, len_dual)
            %DUAL_PROCESS dispatch the dual variables from solution to
            %locations and measures, turn the variables into nonnegative
            %functions along trajectories
            
            %v: coefficients of liouville
            %beta: coefficients of cost (if able)
            %alpha: dual of zeno gaps
            
            %TODO: Fill this in
            
        end
        
        function sol = run(obj, order)
            %RUN the main call, the full peak program at the target order
            
            %Tmax is set by the constraints and the support
            %so it should be the default in obj.loc.supp.Tmax
            sol = run@manager_interface(obj, order, obj.loc.supp.Tmax);
        end
    end
end

