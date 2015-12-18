function hamiltonians = HamiltonianXY( N, h, r )
% Hamiltonian of XY model 
%%
% <<HamiltonianXY.jpg>>
%%
% =
%%
% <<def_Hamiltonian.jpg>>
M = 3 * N - 2;
hamiltonians = cell(M, N);

sx=[0,1;1,0];  sy=[0,-1i;1i,0];
sz=[1,0;0,-1]; id=eye(2);
for m = 1:M
    for n = 1:N
        hamiltonians{m, n } = id;
    end
end

for j = 1:N-1
    m = j;
    hamiltonians{m, j}     = - sx * (1 + r)/2 ;
    hamiltonians{m, j + 1} = - sx * (1 + r)/2 ;  
    m = j + N - 1;
    hamiltonians{m, j}     = - sy * (1 - r)/2 ;
    hamiltonians{m, j + 1} = - sy * (1 - r)/2 ;
    m = j + 2 * (N - 1);
    hamiltonians{m, j}     = - h * sz ; 
end
hamiltonians{M, N}         = - h * sz;
end