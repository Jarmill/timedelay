function em = monom_int(t, x, dv)
    %numerically integrate t^alpha x(t)^beta in span [0, T]
    %em: empirical moments
    %there is probably a more efficient way to implement this
    n_monom = size(dv, 1);
    em = zeros(n_monom, 1);
    for i = 1:n_monom        
        v_curr = (t.^dv(i, 1)) .* prod(x.^(dv(i, 2:end)'), 1);
%         em(i) = trapz(t, v_curr);
        em(i) = simps(t, v_curr);
    end
end