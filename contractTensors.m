function Z = contractTensors(X, rankX, indXContr, Y, rankY, indYContr)
% -------------------------------------------
%  [X,numindX]=contractTensors(X,numindX,indX,Y,numindY,indY)
%
%  �e���\��X�̓Y����indXContr�ƃe���\��Y�̓Y����indYContr���k�񂷂�B
%
%  Parameters
% X = A tensor
%    indXContr     = The index of X to contract
%    rankX         = X��rank
% Y = A tensor
%    indYContr     = The index of Y to contract
%    rankY         = Y��rank
% Returns
%    Z        = The contracted tensor of X and Y
%
% Z�̓Y�����̏��Ԃ�X�̏k�񂵂Ă��Ȃ��Y��������ɗ��āAY�̏k�񂵂Ă��Ȃ��Y��������ɗ���B
% ��)
% X(i_1, ..., i_p, j_1, ..., j_q)
% Y(k_1, ..., k_r, j_1, ..., j_q)
% Z = sum_j X Y = Z(i_1, ..., i_p, k_1, ..., k_r)   
%
% test_testTensors�œ�����e�X�g�����B
%
% ������
% MATLAB��size(T)=[n m 1]�̃e���\���������n*m�̍s��ɕς��Ă��܂��̂ŁA���̊֐��͐���ɓ����Ȃ����Ƃ�����B
% -------------------------------------------
numIndContr = length( indXContr ); % �����length( indYContr )�Ɠ����ł���K�v������B
indXNotContr = 1:rankX;  indXNotContr( indXContr )= [];
indYNotContr = 1:rankY;  indYNotContr( indYContr )= [];

% X��Y�̓Y������move����B���Ƃ�XPermuted��XVectorForm��XMatrixForm�ɂ��ēY�������u���˂�v�B
if rankX ~= 1
    XPermuted = permute(X, [indXNotContr, indXContr]);
else
    XPermuted = X;
end
if rankY ~= 1
    YPermuted = permute(Y, [indYContr, indYNotContr]);
else
    YPermuted = Y;
end
sizeXPermuted = size( XPermuted );
sizeYPermuted = size( YPermuted );    
if rankX == numIndContr && rankY == numIndContr
    XVectorForm = reshape(XPermuted, [1, numel(X)]);
    YVectorForm = reshape(YPermuted, [numel(Y), 1]);
    Z = XVectorForm * YVectorForm;
elseif rankX == numIndContr && rankY ~= numIndContr
    XVectorForm = reshape(XPermuted, [1, numel(X)]);
    YMatrixFormSize =  [sizeYPermuted(1), prod(sizeYPermuted(2 : end) )];
    YMatrixForm = reshape(YPermuted, YMatrixFormSize );
    ZMatrixForm = XVectorForm * YMatrixForm ;
    if length(indYNotContr) == 1 
        % reshape�̑������͒����Q�ȏ�̃x�N�g���łȂ��Ƃ����Ȃ�����l�ߕ�������B
        Z = reshape(ZMatrixForm, [1, sizeYPermuted(2)] ); 
    else
        Z = reshape(ZMatrixForm, sizeYPermuted(2:end) );
    end
elseif rankX ~= numIndContr && rankY == numIndContr
    XMatrixFormSize = [prod(sizeXPermuted(1:rankX - 1)), sizeXPermuted(rankX) ];
    XMatrixForm = reshape(XPermuted, XMatrixFormSize);
    YVectorForm = reshape(YPermuted, [numel(Y), 1]);
    ZMatrixForm = XMatrixForm * YVectorForm;
    if length(indXNotContr) == 1
        % reshape�̑������͒����Q�ȏ�̃x�N�g���łȂ��Ƃ����Ȃ�����l�ߕ�������B
        Z = reshape(ZMatrixForm, [1, sizeXPermuted(1)] );
    else
        Z = reshape(ZMatrixForm, sizeXPermuted(1:length(indXNotContr)));
    end
else % rankX ~= ndimContr && rankY ~= ndimContr
    XMatrixFormSize = [prod(sizeXPermuted(1:rankX - numIndContr)), prod(sizeXPermuted(rankX - numIndContr + 1:end)) ];
    YMatrixFormSize = [prod(sizeYPermuted(1:numIndContr)), prod(sizeYPermuted( numIndContr + 1 : end))  ];
    XMatrixForm = reshape(XPermuted, XMatrixFormSize);
    YMatrixForm = reshape(YPermuted, YMatrixFormSize);
    ZMatrixForm = XMatrixForm * YMatrixForm;
    Z = reshape(ZMatrixForm, [sizeXPermuted(1:length(indXNotContr)), sizeYPermuted(rankY - length(indYNotContr) + 1 : end )] );
end
end
