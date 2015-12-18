function Z = contractTensors(X, rankX, indXContr, Y, rankY, indYContr)
% -------------------------------------------
%  [X,numindX]=contractTensors(X,numindX,indX,Y,numindY,indY)
%
%  テンソルXの添え字indXContrとテンソルYの添え字indYContrを縮約する。
%
%  Parameters
% X = A tensor
%    indXContr     = The index of X to contract
%    rankX         = Xのrank
% Y = A tensor
%    indYContr     = The index of Y to contract
%    rankY         = Yのrank
% Returns
%    Z        = The contracted tensor of X and Y
%
% Zの添え字の順番はXの縮約していない添え字が先に来て、Yの縮約していない添え字が後に来る。
% 例)
% X(i_1, ..., i_p, j_1, ..., j_q)
% Y(k_1, ..., k_r, j_1, ..., j_q)
% Z = sum_j X Y = Z(i_1, ..., i_p, k_1, ..., k_r)   
%
% test_testTensorsで動作をテストした。
%
% ＊注意
% MATLABはsize(T)=[n m 1]のテンソルを勝手にn*mの行列に変えてしまうので、この関数は正常に動かないことがある。
% -------------------------------------------
numIndContr = length( indXContr ); % これはlength( indYContr )と同じである必要がある。
indXNotContr = 1:rankX;  indXNotContr( indXContr )= [];
indYNotContr = 1:rankY;  indYNotContr( indYContr )= [];

% XとYの添え字をmoveする。あとでXPermutedをXVectorFormやXMatrixFormにして添え字を「束ねる」。
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
        % reshapeの第二引数は長さ２以上のベクトルでないといけないから詰め物をする。
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
        % reshapeの第二引数は長さ２以上のベクトルでないといけないから詰め物をする。
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
