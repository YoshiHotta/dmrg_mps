N=100;
D = 8;
J = 1;
h_range = 0.2:0.2:2;

gaps = [];
for h = h_range
    E_GS = minimizeE(J, h, N, D, 1);
    E_excited = minimizeE(J,h,N,D,0);
    gaps = [gaps, E_excited - E_GS ];
    disp(E_excited - E_GS );
end

plot(h_range, gaps);