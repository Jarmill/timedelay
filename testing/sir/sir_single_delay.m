%An SIR model of an epidemic for varying time delay
%
%Only SI, R can be recovered later by integrating I
%Jared Miller, 2/12/2021

%% model parameter
beta = 0.4;
gamma = 0.1;

Tmax = 30;  %max time of simulation (days)
tau = 9;    %incubation period (days)

Trange = linspace(0, Tmax, 200);


I0 = 0.2;   %initial infection rate, constant history
xh = [1-I0; I0];

options = ddeset('AbsTol', 1e-9, 'RelTol', 1e-7, 'Jumps', [0], 'MaxStep', 0.1);

sir_delay = @(t,y,Z) Tmax * [-beta*y(1)*y(2);
            beta*Z(1)*Z(2) - gamma*(y(2))];

sir_history = @(t) xh;

sol = dde23(sir_delay, tau/Tmax, sir_history, Trange/Tmax, options);

plot([-tau, Tmax * sol.x], [xh(2), sol.y(2, :)], 'DisplayName', ['Delay=', num2str(tau)])

title('Infection Rate of Epidemic')
xlabel('time (days)')
ylabel('infection rate')