function Bs = rightNormalize( mps )
% �C�ӂ�MPS��normalize���ꂽright-canonical form�ɕϊ�����B
% ���� : mps = {M_1, M_2, ..., M_N}.
% �o�� : Bs  = {B_1, B_2, ..., B_N}.
% ������N�͗��q���ł���AM_i��rank-3�̃e���\���ł���B
% B_i��Schollwoeck��B-matrix�ł���B
%
% �e�X�g��test_rightNormalized�ōs�����B
N = length( mps );
Bs = cell(1, N);
for i = N:-1:2
    Dleft = size( mps{i}, 1 );
    Dright = size(mps{i},2);
    d = size( mps{i}, 3 );
    MMatrixForm = reshape( permute( mps{i}, [1,3,2]), [Dleft, Dright*d] ) ;
    [U,S,V] = svd(MMatrixForm, 'econ');
    BMatrixForm = V';
    dimS = size(S, 1); % S�͑Ίp�s��
    Bs{i} = permute( reshape(BMatrixForm, [dimS, d, Dright]), [1,3,2]) ;
    MU = contractTensors(mps{i-1}, 3, 2, U, 2, 1) ;
    MU = permute(MU, [1 3 2]);
    mps{i-1} = zeros(size(MU,1), dimS, size(MU,3));
%     for s = 1:dimS
%         mps{i-1}(:,s,:) = MU(:,s,:) * S(s,s);
%     end
    mps{i-1} = permute( contractTensors( MU, 3, 2, S, 2, 1 ), [1 3 2] );

end
B1_unnormalized = mps{1};

% �f�o�b�O�̎��͐U����ς��Ȃ��B�{�Ԃł͋K�i������B
% Bs{1} = B1_unnormalized;
% return;
 
% �k�񂷂�Ƃ��ɍŌ�̓Y�����̎�����1�ɂ��ă����N���������B
B1_unnormalizedMatrixForm = permute(B1_unnormalized, [2 3 1]);
normOfMPS = sqrt( contractTensors(B1_unnormalizedMatrixForm, 2, [1 2], conj(B1_unnormalizedMatrixForm), 2, [1 2]) ) ;
B1_normalized = B1_unnormalized / normOfMPS;
Bs{1} = B1_normalized;
end