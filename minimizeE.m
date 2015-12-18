% Hamiltonian��MPO�Œ�`����BMPO�Ƃ͎��Œ�`�����悤�Ȃ��̂̂��Ƃł���F
% $$H(\sigma', \sigma) = \prod_{l} W_{b_{l-1},b_l}^{\sigma_{l}',\sigma_l}$$
% ����������
% $$H = \sum_{m=1}^{M} h^{m}_1\otimes \cdots \otimes h^{m}_N$$
% mps = { M(1), ..., M(N) }. M{i}��i�Ԗڂ̗��q�ɑΉ�����D *�@D *�@d�̃e���\���B
% �e���\���̓Y������M{i}(p, q, r)�Ƃ���ƁAp, q, r�͂��ꂼ�ꉺ�}�̂P�ԁA�Q�ԁA�R�ԂɑΉ�����F


function E_GS = minimizeE(J, h, N, D, isGS)  
% ������C�W���O���f���̊��G�l���M�[�����߂�֐��B
% ����
%       J, h : ���f���̃p�����[�^
%       N : �T�C�g��
%       D : MPS�̃{���h����
% �Ԃ�l
%       E_GS : isGS=1�Ȃ���G�l���M�[�BisGS=0�Ȃ����N���ʁB

num_sweep = 10;
% N = 20;
% D = 2;
% J = 1;
% h = 1;
HamiltonianMPO = transverseIsingMPO(N, J, h );

d = size(HamiltonianMPO{1}, 3);

% ������Ԃ̓����_���ɂ���B
initialRandomMPS = randomMPS(D, d, N);

% Start from some initial guess for |\phi>, which is right-normalized, i.e. consists of B-matrices only.
initialRandomMPSRightNormalized = rightNormalize( initialRandomMPS ); % B-matrix�ɕϊ�����
M = initialRandomMPSRightNormalized;

% Calculate the R-expressions iteratively for all site positions N - 1 through 1 iteratively.
R = cell(1, N-1); % R{l}��a_{l-1}, b_{-1}, a'_{l-1}���J���Ă���e���\��
% �C�ӂ̃n�~���g�j�A���i�����j��������悤�ɂ��悤�Ƃ�����������Ȃ�������߂��B
%% 
% <<def_Hamiltonian.jpg>>
%R_N_temp = contractTensors( B{N}, 3, 
% for l = N-1:-1:1
%     R_l_temp = contractTensors( R{l+1}, 3, 1, B{l+1}, 3, 2 );
%     
%     R_l_temp_temp = [];
%     for m = 1:M
%         W = hamiltonians{m, l+1};
%         R_l_temp_temp{m} = contractTensors( R_l_temp, 4, [1 4], W, 4, [2 3] )...
%                            + R_l_temp_temp;
%     end
%     R_l_temp = R_l_temp_temp;
%     R{l} = contractTensors( R_l_temp, 4, [1 4], conj(B{l+1}), [2 3] ); 
% end

% MPO�\����Hamiltonian����������B
R_N_temp = contractTensors( M{N}, 3, 3, HamiltonianMPO{N}, 4, 4 );
% R�̓Y�����̒�`�͉��}�̒ʂ�
%%
% <<defOfRIndex.jpg>>

% R{N-1}�����͓��ʂȏ��������āA�Y�����̒�`����}�ɍ��v����悤�ɂ��Ȃ��Ă͂����Ȃ��B
% �Ō�̓Y�����̎������P�Ȃ̂Ń����N��������āArank-6����Ȃ���rank-5�ɂȂ�B
R_N_temp = contractTensors( R_N_temp, 5, 5, conj( M{N} ), 3, 3 );
% �Ō�̓�̓Y�����̎�����1�ɂ��āA�����N��2���Ƃ��B����ɂ����R_N_temp��rank-3�̃e���\���ɂȂ�B
R_N_temp = permute( R_N_temp, [1 3 5 2 4] );
R{N - 1} = R_N_temp; 
for l = N-2:-1:1
    R_l_temp = contractTensors( R{l+1}, 3, 1, M{l+1}, 3, 2 );
    R_l_temp = contractTensors( R_l_temp, 4, [1 4], HamiltonianMPO{l+1}, 4, [2 3] ) ;
    R_1_temp = contractTensors( R_l_temp, 4, [1 4], conj(M{l+1}), 3, [2 3] ); 
    R{l} = R_1_temp;
end

L = cell(1,N -1);
E_GSs = [];
for i_sweep = 1:num_sweep % �����I�Ɏ����𔻕ʂ��ă��[�v�𔲂���悤�ɂ���B

    % Right sweep
    for l=1:N-1
        if l == 1
            Heff = contractTensors( HamiltonianMPO{1}, 4, 2, R{1}, 3, 2);
            Heff = reshape(Heff, [1, size(Heff)]);
        else
            Heff = contractTensors( L{l-1}, 3,2, HamiltonianMPO{l}, 4, 1);
            Heff = contractTensors( Heff, 5, 3, R{l}, 3, 2 );
        end
        HeffMatrixForm = reshape( permute( Heff, [3 2 6 4 1 5] ), ...
            [size(HamiltonianMPO{l}, 3) * size(Heff, 2) * size(R{l}, 3), size(HamiltonianMPO{l}, 4) * size(Heff, 1 ) * size(R{l}, 1) ] );
        MlVectorForm = reshape( permute(M{l}, [3 1 2]), [numel(M{l}), 1] );
        options.v0 = MlVectorForm ;
        
        % �{����P' * Heff *P�Ƃ����������m�ɋ��܂�炵���B
        if isGS
            [MlVectorForm, E_GS] = eigs(HeffMatrixForm, 1,  'sr', options ) ;
        else
            [MlVectorForm2, E_GS2] = eigs(HeffMatrixForm, 2,  'sr', options ) ;
            MlVectorForm = MlVectorForm2(:,2);
            E_GS = E_GS2(2,2);
        end
        
        E_GSs = [E_GSs E_GS];
        M{l} = permute( reshape( MlVectorForm, [size(M{l},3), size(M{l},1), size(M{l},2)] ), [2 3 1 ] ) ;
        
        MlMatrixForm = reshape( permute(M{l}, [3 1 2]), [size(M{l},3) * size(M{l},1), size(M{l},2)] );
        [AlMatrixForm, S, V] = svd(MlMatrixForm, 'econ') ;
        dimS = size( S, 1);
        M{l} = permute( reshape(AlMatrixForm, [size(M{l},3), size(M{l},1),  dimS]), [2 3 1] );
        VdagM = contractTensors( V', 2, 2, M{l+1}, 3, 1 );
        %    M{l+1} = [];
        %     for s = 1:dimS
        %         M{l+1}(s, :, :) = S(s,s) * VdagM(s, :, :);
        %     end
        M{l+1} = contractTensors(S, 2, 2, VdagM, 3, 1 );
        
        if l == 1
            L{1} = contractTensors( M{1}, 3, 3, HamiltonianMPO{1}, 4, 4 );
            L{1} = contractTensors( L{1}, 5, 5, conj( M{1} ), 3, 3 );
            L{l} = permute( L{1}, [2 4 6 1 3 5] ); % rank-3�̃e���\���ɂ���B
        else
            L{l} = contractTensors( L{l-1}, 3, 1, M{l}, 3, 1 );
            L{l} = contractTensors( L{l}, 4, [1 4], HamiltonianMPO{l}, 4, [1 4] );
            L{l} = contractTensors( L{l}, 4, [1 4], conj( M{l} ), 3, [1 3] );
        end
        
    end
    
    
    % Left sweep
    for l = N:-1:2
        if l == N
            Heff = contractTensors( L{N-1}, 3, 2, HamiltonianMPO{N}, 4, 1 );
            Heff = permute(Heff, [1 2 4 5 3]) ;
            HeffMatrixForm = reshape( permute( Heff, [3 2 4 1] ), ...
                [size(HamiltonianMPO{N}, 3) * size(Heff, 2), size(HamiltonianMPO{l}, 4) * size(Heff, 1)]) ;
        else
            Heff = contractTensors( L{l-1}, 3,2, HamiltonianMPO{l}, 4, 1);
            Heff = contractTensors( Heff, 5, 3, R{l}, 3, 2 );
            HeffMatrixForm = reshape( permute( Heff, [3 2 6 4 1 5] ), ...
                [size(HamiltonianMPO{l}, 3) * size(Heff, 2) * size(R{l}, 3), size(HamiltonianMPO{l}, 4) * size(Heff, 1 ) * size(R{l}, 1) ] );
        end
        MlVectorForm = reshape( permute(M{l}, [3 1 2]), [numel(M{l}), 1] );
        options.v0 = MlVectorForm ;
        
        if isGS
            [MlVectorForm, E_GS] = eigs(HeffMatrixForm, 1,  'sr', options ) ;
        else
            [MlVectorForm2, E_GS2] = eigs(HeffMatrixForm, 2,  'sr', options ) ;
            MlVectorForm = MlVectorForm2(:,2);
            E_GS = E_GS2(2,2);
        end
        
        E_GSs = [E_GSs E_GS];
        M{l} = permute( reshape( MlVectorForm, [size(M{l},3), size(M{l},1), size(M{l},2)] ), [2 3 1 ] ) ;
        
        % ��������rightNormalize()�̃R�s�y
        Dleft = size( M{l}, 1 );
        Dright = size( M{l},2);
        d = size( M{l}, 3 );
        MlMatrixForm = reshape( permute( M{l}, [1,3,2]), [Dleft, Dright*d] ) ;
        [U,S,V] = svd(MlMatrixForm, 'econ');
        BlMatrixForm = V';
        dimS = size(S, 1); % S�͑Ίp�s��
        M{l} = permute( reshape(BlMatrixForm, [dimS, d, Dright]), [1,3,2]) ;
        MU = contractTensors(M{l-1}, 3, 2, U, 2, 1) ;
        MU = permute(MU, [1 3 2]);
        M{l-1} = zeros(size(MU,1), dimS, size(MU,3));
        %     for s = 1:dimS
        %         M{l-1}(:,s,:) = MU(:,s,:) * S(s,s);
        %     end
        M{l-1} = permute( contractTensors( MU, 3, 2, S, 2, 1 ), [1 3 2] );
        
        
        % rightNormalize()�̃R�s�y�I��
        
        if l == N
            R_N_temp = contractTensors( M{N}, 3, 3, HamiltonianMPO{N}, 4, 4 );
            R_N_temp = contractTensors( R_N_temp, 5, 5, conj( M{N} ), 3, 3 );
            % �Ō�̓�̓Y�����̎�����1�ɂ��āA�����N��2���Ƃ��B����ɂ����R_N_temp��rank-3�̃e���\���ɂȂ�B
            R_N_temp = permute( R_N_temp, [1 3 5 2 4] );
            R{N - 1} = R_N_temp;
        else
            R_lm1_temp = contractTensors( R{l}, 3, 1, M{l}, 3, 2 );
            R_lm1_temp = contractTensors( R_lm1_temp, 4, [1 4], HamiltonianMPO{l}, 4, [2 3] ) ;
            R_1m1_temp = contractTensors( R_lm1_temp, 4, [1 4], conj(M{l}), 3, [2 3] );
            R{l-1} = R_1m1_temp;
        end     
    end
    
end

% i_sweep = 1:length(E_GSs) ;
% plot(i_sweep, E_GSs);
% disp('calculated:');
% disp(E_GSs(end));
% disp('Exact:')
% disp(exactTransverseIsing(J,h,N));
end

function mps = randomMPS(D, d, N)
mps = cell(1, N);
mps{1} = randn( 1, D, d );
mps{N} = randn( D, 1, d );
for i = 2:N-1
    mps{i} = randn(D,D,d);
end
end
