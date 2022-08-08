
%dual formulation
%come back to this later

%variables
t = sdpvar(1,1);
x0 = sdpvar(2,1);
x1 = sdpvar(2,1);

%parameters and dynamics 
T = 5;
tau = 0.5;
ts = tau/T;

order = 2;
d=2*order;

%functions
v = polynomial([t; x0], d);
phi = polynomial([t; x0], d);
xi = polynomial(t, d); 
 

leb_mom = -(ts).^(1:(d+1)) ./ (1:(d+1));

%support sets
% X0 = 
% Xu = 
X = struct('ineq', 