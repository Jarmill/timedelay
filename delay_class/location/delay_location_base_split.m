classdef delay_location_base_split < location_interface
    %DELAY_LOCATION A location (space) of a dynamical system
    %   includes descriptions of the space as well as measures
    %   used for continuous or discrete time-delay systems
    
    %right now this is a 'base' version without external inputs
    
    properties
        varnames = {'t', 'x', 'x_lag'};
        
        history = [];
        time_slack = [];
    end
    
    methods
        function obj = delay_location_base_split(delay_supp,f, objective)
            %DELAY_LOCATION Construct an instance of this class
            %   Detailed explanation goes here
            
            Tmax = delay_supp.Tmax;
            
            if nargin < 3             
                %by default, no objective
                objective = [];    
            end
            obj@location_interface(delay_supp, f, objective, []);

            %scale down the lags
            if obj.supp.SCALE_TIME
                obj.supp.lags = obj.supp.lags/Tmax;    
                delay_supp.lags = delay_supp.lags/Tmax;
                delay_supp.Tmax = 1;
                if obj.supp.DISCRETE_TIME
                    delay_supp.dt = 1/Tmax;
                end
            end
            
            Nsys = length(obj.f);
            obj.sys = cell(Nsys, 1);
            %subsystems
            
            %history measure
            if obj.supp.DISCRETE_TIME
                obj.history = meas_history_discrete(delay_supp);
            else
                obj.history = meas_history(delay_supp);
            end
            
            
            for i = 1:Nsys                
                %TODO: implement digital system
                if obj.supp.DISCRETE_TIME
                    obj.sys{i} = delay_system_base_discrete_split(obj.supp, obj.f{i});
                else
                    obj.sys{i} = delay_system_base_split(delay_supp, obj.f{i});
                end
            end  

            
%             obj.history = meas_history(delay_supp);

            %slack measure for free terminal time
            if obj.supp.FREE_TERM
                obj.time_slack = meas_time_slack(obj.supp);
            end
        end

        
        %% getters 
        
        %% support
        function supp_con_out = supp_con(obj)
            %SUPP_CON get support constraints of measures
            
            %support of initial, terminal, occupation measures
            supp_con_out = supp_con@location_interface(obj);
            
            %consistency measure
            supp_con_out = [supp_con_out; obj.history.supp()];
            
            %time slack
            if obj.supp.FREE_TERM
                supp_con_out = [supp_con_out; obj.time_slack.supp()];
            end
            
        end
        
        %% constraints        
        function [cons, len_consistency] = consistency_con(obj, d)
            %CONSISTENCY_CON Data consistency constraints between joint
            %occupation measure and marginals of the component measure
            Nlag = length(obj.supp.lags);
            
            cons = [];
            len_consistency = zeros(Nlag, 1);
            for i = 1:Nlag
%                 cons_curr = 0;
%                 for j = 1:length(obj.sys)
%                     cons_curr = 
                %marginal of the joint occupation measure
                term_marg =  obj.sys{1}.meas_occ.mom_monom_marg(i, d);
                
                %shifted marginals of the component measures
                term_shift= obj.sys{1}.meas_occ.mom_monom_shift(i, d);
                
                %shifted marginals of the history measure
                term_history = obj.history.mom_shift(i, d);
                
                %slack in case there is free terminal time
                if obj.supp.FREE_TERM
                    term_slack = obj.time_slack.mom_index(i, d);
                else
                    term_slack = 0;
                end
                
                term_lhs = (term_marg + term_slack);
                term_rhs = (term_shift + term_history);
                
                cons = [cons;  (term_lhs - term_rhs) == 0];
                
                %should all be the same length
                len_consistency(i) = length(term_lhs);
            end                        
        end
        
        function [history_con] = history_con(obj, d)
            %HISTORY_CON Pin down moments of the history component
            %measures. 
            %
            %If X_history is given as a vector or a function
            %handle, set moments of the history component measures to
            %moments of this distribution.
            %
            %Else if there is a free terminal time, ensure that the
            %t-marginal of the history component measures are distributed
            %according to the Lebesgue distribution
            %
            %Else do nothing (fixed terminal time and allowable measures in
            %a set)
            
            %TODO: This is some bad code. untangle the X_history type checks 
            %it should be mom_con_traj or mom_con_free but not both.
                       
            %TODO: Also figure out dual recovery
            X_history = obj.supp.get_X_history_supp();
            history_con = [];
            
            %fixed supplied histories (single given trajectory)
            %either a function handle or a constant point
            if ~isa(X_history, 'supcon')
                history_con = obj.history.history_traj_con(d, X_history);                  
            end

            %free-time multiple histories
            if obj.supp.FREE_TERM && isa(X_history, 'supcon')
                history_con  = obj.history.history_free_con(d, X_history);               
            end
            
            if obj.supp.CONSTANT_HIST
                %shaping constraint to ensure that the histories in
                %X_history are constant in time
                shape_con = obj.shape_constant_con(d);                    
                history_con = [history_con; shape_con];
            end
            
            %TODO: add shaping here
            
            
        end
        
        function [objective, cons_eq, cons_ineq, len_dual] = all_cons(obj, d)
            %ALL_CONS all constraints involving solely location
            %does not include sum of mass of initial measures
            %Output:
            %   cons_eq: equality constraints
            %   cons_ineq: inequality constraints (objective)
            
            %gather all constraints together
            %liouville
            liou = obj.liou_con(d);
            len_liou = length(liou);
            
            %box absolute continuity
            [abscont_box, len_abscont] = obj.abscont_box_con(d);
            
            %objective
            [objective, cons_ineq] = obj.objective_con();
            
            %data consistency
            [consistency, len_consistency] = obj.consistency_con(d);
            

            [cons_history] = obj.history_con(d);
            
            %package up the output
            len_dual = struct;
            len_dual.v = len_liou;
            len_dual.zeta = len_abscont;
            len_dual.phi  = len_consistency;
            len_dual.beta = length(cons_ineq);
            
            %ensure this is the correct sign
%             cons_eq = [-liou; abscont_box; consistency]==0;                        
            cons_eq = [-[liou; abscont_box]==0; consistency];                        
            
            %history constraints perform moment substitution 
            %this is permissible, because we know exactly what moments are
            %being substituted.
            cons_eq = [cons_eq; cons_history];
        end     
        
        function [abscont_box, len_abscont] = abscont_box_con(obj, d)
            %absolute continuity constraints for the box
            %not used here
            abscont_box = [];
            len_abscont = 0;
        end
        
        function [shape_con, len_shape_con] = shape_constant_con(obj, d)
            %shaping constraint to ensure that the histories are constant
            %in time. Adding more complicated shaping constraints will be
            %the subject of future work.
            shape_con = [];
            
            %Lie derivative of history
            Lhist = obj.history.shaping_mom_const(d);
            
            %monomials for initial measure
            monom = mmon([obj.vars.t; obj.vars.x], d);
            
            monom_tau = subs(monom, obj.vars.t, -max(obj.supp.lags));
            monom_0 = subs(monom, obj.vars.t, 0);
            
            vars_reduced = struct('t', obj.vars.t, 'x', obj.vars.x);
            mom_tau = obj.init.var_sub_mom(vars_reduced, monom_tau);
            mom_0 = obj.init.var_sub_mom(vars_reduced, monom_0);

            
            shape_con = (mom_tau + Lhist - mom_0 == 0);

            len_shape_con = length(shape_con);
            
            
        end
        
        function [len_out] = len_eq_cons(obj)
            %LEN_EQ_CONS Number of equality constraints strictly in this
            %location 
            len_out = obj.len_dual.v + sum(obj.len_dual.zeta)+ sum(obj.len_dual.phi);
        end
        
        function [obj_max, obj_con] = objective_con(obj, objective)
            %OBJECTIVE_CON deal with the objective, which may be maximin
                                    
            %TODO: This should maybe go in the manager
            %The current implementation is only for peak estimation
            
            %TODO: include support for putting objectives on initial and
            %occupation measures as well as the terminal measure
            if nargin == 1
                objective = obj.objective;
            end
                                    
            obj_con = [];
            
            %TODO: add 'th' to this 
            var_end = obj.var_index(obj.vars, {'t', 'x'});
            if isempty(objective)
                obj_max = 0;
            elseif length(objective) == 1    
                obj_subs = obj.term.var_sub_mom(var_end, objective);
                obj_max = (obj_subs);                            
            else
                obj_subs = obj.term.var_sub_mom(var_end, objective);
                q_name = ['q_', num2str(obj.id)];
                mpol(q_name, 1, 1);
                q = eval(q_name);
                muq = meas(q);
                obj.cost_q = q;
                
                obj_max = mom(q);
                obj_con = [mass(q) == 1; (mom(q) <= obj_subs);];
            end            
        end
        
        %% Dual variables
        function obj = dual_process(obj, d, rec_eq, rec_ineq, gamma)
            %TODO: fill this in
            %pass right now
        end
    end
end

