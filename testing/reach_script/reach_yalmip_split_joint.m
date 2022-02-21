%% problem properties
tau = 0.3;
K0 = 0.5;
K1 = 4;
T = 1;
boxlim = 1;
order = 5;
d = 2*order;

% xh0_max = -0.6;
xh0_max = -0.9;
xh0_min = -1;

%% sdpvariables
t = sdpvar(1,1);
x = sdpvar(1, 1);
x0 = sdpvar(1,1);
x1 = sdpvar(1,1);

f = -K0*x0-K1*x1;

[v, cv, mv] = polynomial([t; x], d);
[w, cw, mw] = polynomial([x], d);
[phi1, cphi1, mphi1] = polynomial([t; x], d);

v0 = replace(v, t, 0);
vT = replace(v, t, T);

cons = [];
coeff = [cv; cw; cphi0; cphi1];


%% support sets
Xbase = [boxlim^2-x^2];
XT = struct('ineq', Xbase, 'eq', []);
Xhbase = (xh0_max-x)*(x-xh0_min);
X0 = struct('ineq', [Xhbase], 'eq', []);

%history [-tau, 0]
Omeganeg1 = struct('ineq', [(0-t)*(t-(-tau)); Xhbase], 'eq', []);
%present [0, T-tau]
% Omega0 = struct('ineq', [t*((T-tau) - t); Xbase], 'eq', []);
%present [T-tau, T]
% Omega1 = struct('ineq', [(t-(T-tau))*(T-t); Xbase], 'eq', []);

Xall0 = struct('ineq', [t*(T-tau-t); boxlim^2-[x0; x1].^2], 'eq', []);
Xall1 = struct('ineq', [(t-(T-tau))*(T-t); boxlim^2-[x0; x1].^2], 'eq', []);
% Xall = struct('ineq', [t*(T-t); boxlim^2-[x0; x1].^2], 'eq', []);

%% cost
leb_mom = LebesgueBoxMom( d, [-boxlim; boxlim], 1);
cost = cw'*leb_mom;

%% constraints

%initial constraint
[p0, cons0, Gram0] = psatz(-v0, X0, order, x);

%terminal constraint
[pT, consT, GramT] = psatz(vT + w -1, XT, order, x);

%lie constraint
v_lie = replace(v, x, x0);
Lie = jacobian(v_lie, t) + f * jacobian(v_lie, x0);
Lie_term = -Lie + replace(phi0, x, x0) + replace(phi1, x, x1);
% [plie, conslie, Gramlie] = psatz(-Lie_term, Xall, order, [t; x0; x1]);

%Omega 0 constraint [0, T-tau]
phi1_shift = replace(phi1, t, t+tau);
[pOm0, consOm0, GramOm0] = psatz(Lie_term + phi1_shift, Omega0, order, [t; x]);

%Omega 1 constraint [T-tau, 1]

[pOm1, consOm1, GramOm1] = psatz(phi0, Omega1, order, [t; x]);

%Omega -1 constraint [-tau, 1]
% phi0_shift = replace(phi0, t, t+tau);
[pOmneg, consOmneg, GramOmneg] = psatz(phi1_shift, Omeganeg1, order, [t; x]);

%w positivity constraint
[pw, consw, Gramw] = psatz(w, XT, order, x);

%% wrap up and solve the program


cons = [cons0:'Time 0'; consT:'Time T'; conslie:'Lie Derivative'; ...
    consOm1:'Time [T-tau, T]'; consOm0:'Time[0, T-tau]'; consOmneg:'Time [-tau, 0]'; ...
    consw:'w nonneg'];

opts = sdpsettings('solver', 'mosek');

sol = optimize(cons, cost, opts);


%% recover solutions
value(cost)
cw_rec = value(cw)'
w_rec = value(cw)'*mw;
v_rec = value(cv)'*mv;
phi0_rec = value(cphi0)'*mphi0;
phi1_rec = value(cphi1)'*mphi1;

w_func = polyval_func(w_rec, x);

%% plot the estimated ROA
figure(1)
clf
Npts = 101;
x_sample = linspace(-boxlim, boxlim, Npts);
w_sample = w_func(x_sample);
hold on
plot(x_sample, w_sample, 'LineWidth', 3);
plot([-boxlim, boxlim], [1,1], 'k', 'LineWidth', 3)
xlim([-boxlim, boxlim])
xlabel('x')
ylabel('w(x)')
title('Estimated T-reachable set above black line', 'FontSize', 14)
hold off
