function E_GS = exactTransverseIsing(J, h, N)
E_GS = J * ( 1 - csc(pi / 2 / (2 * N + 1 )) );
end