mpol('t', 1, 1);
mpol('x', n, 1);
lags = [1; 2];

Tmax = 5;
supp = [t*(Tmax-t) >= 0; sum(x) <= 1; x >= 0];

vars = struct('t', t, 'x', x);
cm = meas_base(vars, supp);

cm.supp

cm.mom_monom_shift(-2, 3)