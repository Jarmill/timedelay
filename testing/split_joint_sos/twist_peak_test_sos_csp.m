
%dual formulation for peak estimation of dde
%the Twist system with time delay
%
%CSP: Apply Correlative Sparsity to the Lie Derivative constraint
%% setup

%variables
t = sdpvar(1,1);
x0 = sdpvar(3,1);
x1 = sdpvar(3,1);
beta = sdpvar(1, 1);

%% system properties
Tmax = 8;
% T = 5;
% tau = 0.5;
tau = 0.5;
ts = tau/Tmax;

C0 = [-0.5; 0; 0];
R0 = 0.2;

A_true = [-1 1 1; -1 0 -1; 0 1 -2];
B_true = [-1 0 -1;
          0 1 1;
          1 1 0]/2;

f = A_true*x0 - B_true*(4*(x0).*(x1.^2) - 3*x0);
      
      
%set the sos tightening order      
order = 2;



d = 2*order;


%objective
p = x0(3);

%functions
[v, cv] = polynomial([t; x0], d);
[phi, cphi] = polynomial([t; x0], d);
[xi, cxi] = polynomial(t, d); 
gamma = sdpvar(1,1);
beta = sdpvar(1,1);

leb_mom = LebesgueBoxMom( d, [-ts; 0], 1);

Lv = jacobian(v, t) + Tmax * jacobian(v,x0)*f;

v0 = replace(v, t, 0);


%% support sets
% X0 = 
% Xu = 
box = [-1.1, 0.75; -1.5, 2; -1, 1.2];
X0 = [(x0 - box(:, 1)).*(box(:, 2) - x0)];
X1 = replace(X0, x0, x1);
tp = t*(1-t);
tn = (t + ts)*(-t);

t0 = t*((1-ts) - t);
t1 = (t - (1-ts))* (1-t);

H0 = R0^2 - sum((x0-C0).^2);
H = [H0; tn];

phi1 = replace(phi, x0, x1);
phishift = replace(phi, t, t+ts);


%% psatz constraints

[p0, cons0, coeff0] = constraint_psatz(-v0 + gamma, H0, [x0], d);
% 
[pc, consc, coeffc] = constraint_psatz(-p + v, [tp; X0], [t; x0], d);

% [plie0, conslie0, coefflie0] = constraint_psatz(Lv , [t0; X0], [t; x0], d+2);
% [plie0, conslie0, coefflie0] = constraint_psatz(Lv , [t0; X0; X1], [t; x0; x1], d+2);
% 
[ph, consh, coeffh] = constraint_psatz(-phishift + xi, H, [t; x0], d);
% 
[plie0, conslie0, coefflie0] = constraint_psatz(-Lv - phishift + phi1, [t0; X0; X1], [t; x0; x1], d+2);
% 
[plie1, conslie1, coefflie1] = constraint_psatz(-Lv + 0 + phi1, [t1; X0; X1], [t; x0; x1], d+2);
% 
[pfree, consfree, coefffree] = constraint_psatz(phi, [tp; X0], [t; x0], d);

consbeta = [sum(beta)==1; beta>=0];


%% accumulate constraints
objective = (gamma + leb_mom'*cxi);

% coeff = [cv; cphi; cxi; gamma;
%         coeff0; coeffc; coeffh; coefflie0;coefflie1; coefffree];
%     

% cons = [cons0; consc; consh; conslie0; conslie1; consfree];

%progressively adding constraints
coeff = [cv; gamma; cxi; cphi; beta;
    coeff0; coeffc; coeffh; coefflie0; coefflie1; coefffree];
cons = [cons0; consc; consh; conslie0; conslie1; consfree; consbeta];

%% solve
opts = sdpsettings('solver','mosek');


[sol, monom, Gram, residual] = solvesos(cons, objective, opts, [coeff]);

value(objective)
