classdef delay_support < loc_support
    %DELAY_SUPPORT Support of a location of a system with bounded
    %discrete time delays
    
    
    properties
        lags;       %An array of lags for bounded discrete time delays (double)
        X_history;  %support of history component occupation measures 
                    %Type options:
                    %   @supcon:    support for multiple histories
                    %   double:     Constant in time, supported at a single
                    %               state point
                    %   function:   A handle evaluating the single history 
                    %               as  X_history(t,x)
                    %TODO: allow for multiple histories (?)
        dt = 0;
        DISCRETE_TIME=0;
        CONSTANT_HIST=0; %shaping constraint to ensure that the histories 
                         %inside X_history are constant in time
        %vars will include 'x_lag' and possibly 'w_lag'
        
    end
    
    methods
        function obj = delay_support(vars,loc_ref)
            %DELAY_SUPPORT Construct an instance of this class
            %   Detailed explanation goes here
            if nargin == 1
                loc_ref = [];
            end
            
            obj@loc_support(vars, loc_ref);
            
            if nargin > 1                
                obj.lags = loc_ref.lags;
                obj.X_history = subs_vars(loc_ref.X_history, loc_ref.vars.x, obj.vars.x);                    
            end
        end
        
        %% support getters        
        function X_joint = get_X_joint(obj, X_in)
            %Get state support of joint occupation measure
            
            if nargin == 1
                X_in = obj.X;
            else
                X_in = obj.X;
            end
            
            X_joint = X_in;
            
            for i = 1:length(obj.lags)
                X_curr = subs_vars(X_in, obj.vars.x, obj.vars.x_lag(:, i));
                X_joint = [X_joint; X_curr];
            end            
        end
        
        function X_sys = get_X_sys_single(obj, X_sys_in)
            if isempty(X_sys_in)
                X_sys = obj.get_X_joint();
            else
                X_sys = obj.get_X_joint(X_sys_in);
            end
        end
        
        function X_sys = get_X_sys_ind(obj, ind)
            %get support of subsystem's joint occupation measure
            if isempty(ind) || isempty(obj.X_sys) || isempty(obj.X_sys{ind})
                X_sys = obj.get_X_joint();
            else
                if iscell(obj.X_sys)
                    X_sys = obj.get_X_joint(obj.X_sys{ind});
                else
                    X_sys = obj.get_X_joint(obj.X_sys);
                end
            end
        end

        function X_hist = get_X_history_supp(obj)
            %constraint on support of history
            if isempty(obj.X_history)
                %by default use the support of the initial measure 
                %at time t=0
                X_hist = obj.get_X_init();
            else
                X_hist = obj.X_history;
            end
            
            %return support of history if X_sys is a supcon
            %otherwise ignore support. substitute in moments of history 
            %measures later            
            
        end

        %TODO: shaping constraints on histories
        %such as weak derivative constraints
        
        function W_joint = get_W_joint(obj)
            %Get input support of joint occupation measure
            %useful for dead-time control
            %TODO
            W_joint = 0;
        end
    end
end

