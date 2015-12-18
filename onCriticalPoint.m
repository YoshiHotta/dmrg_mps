N=100;
Ds = 2:20;
isGS = 1;
J = 1;
h = 1; 
errors = [];
for D = Ds
    E_GS = minimizeE(J, h, N, D, isGS);
    E_GS_exact = exactTransverseIsing(J, h ,N);
    errors = [errors, abs(E_GS - E_GS_exact)];
    disp(abs(E_GS - E_GS_exact));
end

plot(Ds, errors);