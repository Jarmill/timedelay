%An SIR model of an epidemic for varying time delay
%
%
%Jared Miller, 2/12/2021

%% Set up model parameters
p = struct;
p.beta = 0.4;
p.gamma = 0.08;


N_lags = 5;
lag_list = 1:N_lags; %delay time


Tmax = 40; %max time of simulation


%initial infection rate
I0= 0.2;
x0 = [1-I0; I0; 0];

%% Solve Differential Equations

%ode
Trange = linspace(0, Tmax, 200);
sir_curr =@(t,x) sir(t,x,p);
[T, X] = ode45(sir_curr, Trange, x0);

I_max = zeros(N_lags+1, 1);
I_max(1) = max(X(:, 2));

%dde
dde_opts = ddeset('MaxStep', 0.5);
sir_delay_curr = @(t,x,Z) sir_delay(t,x,Z,p);
sir_history_curr = @(t) x0;

sol_lags = cell(N_lags, 1);

for i = 1:N_lags
    sol_lags{i} = dde23(sir_delay_curr, lag_list(i), sir_history_curr, Trange, dde_opts);
    I_max(i+1) = max(sol_lags{i}.y(2, :));
end

%% Plot

figure(1)
clf
hold on



plot(T, X(:, 2), 'DisplayName', 'Delay=0')

for i = 1:N_lags
    plot(sol_lags{i}.x, sol_lags{i}.y(2, :), 'DisplayName', ['Delay=', num2str(lag_list(i))])
end
title('Infection Rate of Epidemic')
xlabel('time (days)')
ylabel('infection rate')

legend('location', 'northeast')

%% Functions
function dydt = sir(t, y, p)
    dydt = [-p.beta*y(1)*y(2);
            p.beta*y(1)*y(2) - p.gamma*(y(2));
            p.gamma*y(2)];
end

function dydt = sir_delay(t, y, Z, p)
    %Z = y(t - lag)
    dydt = [-p.beta*y(1)*y(2);
            p.beta*Z(1)*Z(2) - p.gamma*(y(2));
            p.gamma*y(2)];
end

function s = sir_history(t, x0)
    s = x0;
end
