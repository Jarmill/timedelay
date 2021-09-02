xh0 = 0.5;
tau = 14;
T = 300;
xstar = 0.75;

r = 0.15;
K = 1;

options = ddeset('AbsTol', 1e-11, 'RelTol', 1e-9, 'Jumps', 0, 'MaxStep', 1/T);

% sol_open = dde23(@(t,y,z) r*y*(xstar-z/K), [tau],@(t) xh0,[0,T]);

sol_scaled = dde23(@(t,y,z) T*r*y*(xstar-z/K), [tau/T],@(t) xh0,[0,1], options);

figure(1)
clf

% plot(sol_open.x, sol_open.y)
plot(sol_scaled.x, sol_scaled.y)

