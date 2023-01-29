classdef delay_location_prop_split < location_interface
    %DELAY_LOCATION_PROP_SPLIT  A location (space) of a dynamical system
    %   includes descriptions of the space as well as measures
    %   used for continuous or discrete time-delay systems
    %
    % proportional delay: x'(t) = f(t, x(t), x(s t)) for s in [0, 1].
    %
    %
    % there is no history measure involved here. Negative time is never
    % needed in the proportional delay model
    properties
        varnames = {'t', 'x', 'x_lag'};
        
        time_slack = [];
    end
    
    methods
        function obj = delay_location_prop_split(delay_supp,f, objective)
            %DELAY_LOCATION_PROP_SPLIT Construct an instance of this class
            %   Detailed explanation goes here
            
            if nargin < 3             
                %by default, no objective
                objective = [];    
            end
            obj@location_interface(delay_supp, f, objective, []);

            %scale down the lags
            if obj.supp.SCALE_TIME                
                delay_supp.Tmax = 1;
            end
            
            Nsys = length(obj.f);
            obj.sys = cell(Nsys, 1);
            
            %subsystems
            for i = 1:Nsys                
                obj.sys{i} = delay_system_base_split(delay_supp, obj.f{i}, @meas_joint_prop_split);                
            end  

            %slack measure for free terminal time
            if obj.supp.FREE_TERM
                obj.time_slack = meas_time_slack(obj.supp);
            end
            
        end
        
        %% support
        function supp_con_out = supp_con(obj)
            %SUPP_CON get support constraints of measures
            
            %support of initial, terminal, occupation measures
            supp_con_out = supp_con@location_interface(obj);

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
                %marginal of the joint occupation measure
                term_marg =  obj.sys{1}.meas_occ.mom_monom_marg(i, d);
                
                %shifted marginals of the component measures
                term_shift= obj.sys{1}.meas_occ.mom_monom_scale(i, d);
                
                
                %slack in case there is free terminal time
                if obj.supp.FREE_TERM
                    term_slack = obj.time_slack.mom_index(i, d);
                else
                    term_slack = 0;
                end
                
                term_lhs = (term_marg + term_slack);
                term_rhs = term_shift;
                
                cons = [cons;  (term_lhs - term_rhs) == 0];
                
                %should all be the same length
                len_consistency(i) = length(term_lhs);
            end                        
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

            %package up the output
            len_dual = struct;
            len_dual.v = len_liou;
            len_dual.zeta = len_abscont;
            len_dual.phi  = len_consistency;
            len_dual.beta = length(cons_ineq);
            
            %ensure this is the correct sign
%             cons_eq = [-liou; abscont_box; consistency]==0;                        
            cons_eq = [-[liou; abscont_box]==0; consistency];                        
            
        end     
        
        function [abscont_box, len_abscont] = abscont_box_con(obj, d)
            %absolute continuity constraints for the box
            %not used here
            abscont_box = [];
            len_abscont = 0;
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

