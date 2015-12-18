function HamiltonianMPO = transverseIsingMPO(N, J, h)
% 横磁場イジングのハミルトニアンをMatrix product operator形式で返す
%%
% <<transverseIsingMPO.jpg>>
sx=[0,1;1,0] ; sy=[0,-1i;1i,0] ; sz=[1,0;0,-1] ; % Spin operators.

W = zeros(3,3,2,2);
W(1,1,:,:) = eye(2);
W(2,1,:,:) = - J * sx;
W(3,1,:,:) = - h * sz;
W(3,2,:,:) = sx;
W(3,3,:,:) = eye(2);

W1 = zeros(1,3,2,2);
W1(1,1,:,:) = - h * sz;
W1(1,2,:,:) = sx;
W1(1,3,:,:) = eye(2);

WN = zeros(3,1,2,2);
WN(1,1,:,:) = eye(2);
WN(2,1,:,:) = -J * sx;
WN(3,1,:,:) = -h * sz;

HamiltonianMPO = cell(1,N);
HamiltonianMPO{1} = W1;
HamiltonianMPO{N} = WN;
for i = 2:N-1
    HamiltonianMPO{i} = W;
end

end