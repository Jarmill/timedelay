load('marg_confuse.mat')

%% create polynomials
mset clear
mpol('x0', 1, 1)
mpol('x1', 1, 1)
mpol('t', 1, 1)

order = 2;

d = 2*order;
dv = genPowGlopti(2,d);

v = mmon([t; x0], 0, d);



K0 = 2;
K1 = 3;
f = T*(-K0*x0 - K1 * x1);


liou_term = diff(v, t) + f * diff(v, x0);


%% evaluate polynomials
t_curr = tsample;
x0_curr = x_traj(1, :);
x1_curr = x_delay(1, :);

mom_init = eval(v, [t; x0], [0; x0_curr(1)]);
mom_term = eval(v, [t; x0], [1; x0_curr(end)]);


mom_d_liou = eval(liou_term, [t; x0; x1], [tsample/T; x0_curr; x1_curr]);

mom_d_mix = eval(liou_term, [t; x0; x1],  [tsample/T; x0_curr; x_delay(2, :)]);

mom_liou = simps(tsample/T, mom_d_liou, 2);

mom_mix = simps(tsample/T, mom_d_mix, 2);


%this term should be zero for every monomial
liou_eq = mom_term - mom_init - mom_liou;