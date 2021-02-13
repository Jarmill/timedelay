%An SIR model of an epidemic for varying time delay
%
%Only SI, R can be recovered later by integrating I
%Jared Miller, 2/12/2021


%% Set up model parameters
p = struct;
p.beta = 0.4;
p.gamma = 0.1;

N_lags = 9;
lag_list = 1:N_lags; %delay time


Tmax = 30; %max time of simulation


%initial infection rate
I0= 0.2;
x0 = [1-I0; I0];    %initial condition at t=0
% xh0 = [1; 0];       %constant history
xh0 = x0;

%% Solve Differential Equations

%ode
Trange = linspace(0, Tmax, 200);
sir_curr =@(t,x) Tmax * sir(t,x,p);
[T, X] = ode45(sir_curr, Trange/Tmax, x0);

I_max = zeros(N_lags+1, 1);
I_max(1) = max(X(:, 2));

%dde
options = ddeset('AbsTol', 1e-11, 'RelTol', 1e-9, 'Jumps', 0, 'MaxStep', 1/Tmax);
sir_delay_curr = @(t,x,Z) Tmax * sir_delay(t,x,Z,p);
sir_history_curr = @(t) x0*(t ==0) + xh0*(t ~= 0);

sol_lags = cell(N_lags, 1);

for i = 1:N_lags
    sol_lags{i} = dde23(sir_delay_curr, lag_list(i)/Tmax, sir_history_curr, Trange/Tmax, options);
    I_max(i+1) = max(sol_lags{i}.y(2, :));
end

%% Plot

figure(1)
clf
hold on



plot(Tmax * T, X(:, 2), 'DisplayName', 'Delay=0')

for i = 1:N_lags
    plot(Tmax * sol_lags{i}.x, sol_lags{i}.y(2, :), 'DisplayName', ['Delay=', num2str(lag_list(i))])
end
title('Infection Rate of Epidemic')
xlabel('time (days)')
ylabel('infection rate')

legend('location', 'northeast')

%% Functions
function dydt = sir(t, y, p)
    dydt = [-p.beta*y(1)*y(2);
            p.beta*y(1)*y(2) - p.gamma*(y(2))];
end

function dydt = sir_delay(t, y, Z, p)
    %Z = y(t - lag)
    dydt = [-p.beta*y(1)*y(2);
            p.beta*Z(1)*Z(2) - p.gamma*(y(2))];
end

function s = sir_history(t, x0)
    s = x0;
end
