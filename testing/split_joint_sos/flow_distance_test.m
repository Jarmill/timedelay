
%dual formulation for peak estimation of dde

%% setup

%variables
t = sdpvar(1,1);
x0 = sdpvar(2,1);
x1 = sdpvar(2,1);
y = sdpvar(2, 1);

%parameters and dynamics 
% Tmax = 5; %safe at order 3, p^* = -0.0343 < 0
          %dist_rec = 0.3046, order 3
          
%the safety margin fails, but the distance is successful
Tmax = 8; 


%Tmax = 8 
% order = 1; %1.1897e-04
% order = 2; %4.0420e-04
%order=3;% 0.1572, 
%order=4; % 0.1820
order=5; % 0.1820
tau = 0.5;
ts = tau/Tmax;

f = [x0(2); -x1(1) + (1/3).*x0(1).^3- x0(2)];
% f = [x0(2); -x0(1) - x0(2)];

% f = [x0(2); -x0(1) - x0(2) - x1(1)];

%initial set
C0 = [1.5; 0];
R0 = 0.4;

%unsafe set
Cu = [-0.5; -1];
Ru = 0.5;
theta_c = 5*pi/4;       %p* = -0.1417, beta = [0, 1], safe

w_c = [cos(theta_c); sin(theta_c)];
c1f = Ru^2 - (x0(1) - Cu(1)).^2 - (x0(2) - Cu(2)).^2;
% w_c = [cos(theta_c); sin(theta_c)];
c2f = w_c(1)*(x0(1) - Cu(1)) + w_c(2) * (x0(2) - Cu(2)); 

cf = [c1f; c2f];
cfy = replace(cf, x0, y);

c = sum((x0-y).^2);


d=2*order;

%objective
p = x0(2);

%functions
[v, cv] = polynomial([t; x0], d);
[w, cw] = polynomial([x0], d);
[phi, cphi] = polynomial([t; x0], d);
[xi, cxi] = polynomial(t, d); 
gamma = sdpvar(1,1);


leb_mom = LebesgueBoxMom( d, [-ts; 0], 1);

Lv = jacobian(v, t) + Tmax * jacobian(v,x0)*f;

v0 = replace(v, t, 0);


%% support sets
% X0 = 
% Xu = 
% box = [-1.25, 2.5; -1.25, 1.5];
% box = [-1.25, 2.5; -1.5, 1.1];
box = [-1, 2.5; -1.5, 1.1];
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

[p0, cons0, coeff0] = constraint_psatz(v0 - gamma, H0, [x0], d);
% 
[pc, consc, coeffc] = constraint_psatz(w- v, [tp; X0], [t; x0], d);

[pdist, consdist, coeffdist] = constraint_psatz(c - w, [X0; cfy], [x0; y], d);

% [plie0, conslie0, coefflie0] = constraint_psatz(Lv , [t0; X0], [t; x0], d+2);
% [plie0, conslie0, coefflie0] = constraint_psatz(Lv , [t0; X0; X1], [t; x0; x1], d+2);
% 
[ph, consh, coeffh] = constraint_psatz(phishift - xi, H0, [t; x0], d);
% 
[plie0, conslie0, coefflie0] = constraint_psatz(Lv + phishift - phi1, [t0; X0; X1], [t; x0; x1], d+2);
% 
[plie1, conslie1, coefflie1] = constraint_psatz(Lv + 0 - phi1, [t1; X0; X1], [t; x0; x1], d+2);
% 
[pfree, consfree, coefffree] = constraint_psatz(-phi, [tp; X0], [t; x0], d);


%% accumulate constraints
objective = -(gamma + leb_mom'*cxi);

% coeff = [cv; cphi; cxi; gamma;
%         coeff0; coeffc; coeffh; coefflie0;coefflie1; coefffree];
%     

% cons = [cons0; consc; consh; conslie0; conslie1; consfree];

%progressively adding constraints
coeff = [cv; gamma; cxi; cphi; cw; 
    coeff0; coeffc; coeffh; coefflie0; coefflie1; coefffree; coeffdist];
cons = [cons0; consc; consh; conslie0; conslie1; consfree; consdist];

%adding the cxi term causes infeasibility

% %no delay
% coeff = [cv; gamma;
%     coeff0; coeffc; coefflie0];
% cons = [cons0; consc; conslie0; cxi==0; cphi==0];
%% solve
opts = sdpsettings('solver','mosek');

% opts.sos.model = 2;

[sol, monom, Gram, residual] = solvesos(cons, objective, opts, [coeff]);


dist_rec = sqrt(-value(objective))
% value(objective)
