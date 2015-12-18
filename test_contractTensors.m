function test_contractTensors
% test case 1 
X = randi(100, 1,5);
Y = randi(100, 1,5);
if contractTensors(X, 1, 1, Y, 1, 1) == X*Y'
    disp('test case 1 : succeed');
else
    disp('test case 1 : fail');
end

% test case 2 
X = randi(100, 1,5);
Y = randi(100,5,10);
if contractTensors(X, 1, 1, Y, 2, 1) == X*Y
    disp('test case 2 : succeed');
else
    disp('test case 2 : fail');
end

% test case 3 
X = randi(100,1,5);
Y = randi(100,10, 5);
if contractTensors(X, 1, 1, Y, 2, 2) == X*Y'
    disp('test case 3 : succeed');
else
    disp('test case 3 : fail');
end
    
% test case 4 
X = randi(100,10, 5);
Y = randi(100,5, 1);
if contractTensors(X, 2, 2, Y, 1, 1) == (X*Y)'
    disp('test case 4 : succeed');
else
    disp('test case 4 : fail');
end

% test case 5 
X = randi(100,1,5);
Y = randi(100,10, 5);
if contractTensors(X, 1, 1, Y, 2, 2) == X*Y'
    disp('test case 5 : succeed');
else
    disp('test case 5 : fail');
end

% test case 6
X = randi(100, 1,5);
Y = randi(100, 1,5);
if contractTensors(X, 2, 2, Y, 2, 2) == X*Y'
    disp('test case 6 : succeed');
else
    disp('test case 6 : fail');
end

% test case 7
X = randi(100, 1);
Y = randi(100, 1);
if contractTensors(X, 1,1,Y,1,1)
    disp('test case 7 : succeed');
else
    disp('test case 7 : fail');
end

% test case 8
X = randi(100, 1);
Y = randi(100, 1);
if contractTensors(X, 2,1,Y,2,1)
    disp('test case 8 : succeed');
else
    disp('test case 8 : fail');
end

% test case 9
X = randi(100, 4,5,6);
Y = randi(100, 6,5,4);
if contractTensors(X, 3, 1, Y, 3, 3) == ansTensor(X, 3, 1, Y, 3, 3)
    disp('test case 9 : succeed');
else
    disp('test case 9 : fail');
end

% test case 10
X = randi(100, 4,5,6);
Y = randi(100, 6,5,4);
if contractTensors(X, 3, [1,3], Y, 3, [3,1]) == ansTensor(X, 3, [1,3], Y, 3, [3,1]) 
    disp('test case 10 : succeed');
else
    disp('test case 10 : fail');
end

% test case 11
X = randi(100, 4,5,6);
Y = randi(100, 6,5,4);
if contractTensors(X, 3, [1,3,2], Y, 3, [3,1,2]) == ansTensor(X, 3, [1,3,2], Y, 3, [3,1,2]) 
    disp('test case 11 : succeed');
else
    disp('test case 11 : fail');
end

% test case 12
X = randi(100, 4,5,6);
Y = randi(100, 6,5,4);
if contractTensors(X, 3, [2,3,1], Y, 3, [2,1,3]) == ansTensor(X, 3, [2,3,1], Y, 3, [2,1,3]) 
    disp('test case 12 : succeed');
else
    disp('test case 12 : fail');
end

% test case 13
X = randi(100, 4,5,6);
Y = randi(100, 6,5,4);
if contractTensors(X, 3, 1, Y, 3, 3) == ansTensor(X, 3, 1, Y, 3, 3)  
    disp('test case 13 : succeed');
else
    disp('test case 13 : fail');
end

% test case 14
X = randi(100, 4,5,6);
Y = randi(100, 6,5,4);
if contractTensors(X, 3, 2, Y, 3, 2) == ansTensor(X, 3, 2, Y, 3, 2)  
    disp('test case 14 : succeed');
else
    disp('test case 14 : fail');
end

% test case 15
X = randi(100, 4,5,6);
Y = randi(100, 6,5,4);
if contractTensors(X, 3, [1,2], Y, 3, [2,3]) == ansTensor(X, 3, [1,2], Y, 3, [2,3]) 
    disp('test case 15 : succeed');
else
    disp('test case 15 : fail');
end

% test case 16
X = randi(100, 4,5);
Y = randi(100, 5,4);
if contractTensors(X, 2, [1,2], Y, 2, [2,1]) == ansTensor(X, 2, [1,2], Y, 2, [2,1]) 
    disp('test case 16 : succeed');
else
    disp('test case 16 : fail');
end


% test case 17
X = randi(100, 1, 2, 3);
Y = randi(100, 1, 2, 3);

if contractTensors(X, 3, [2,3], Y, 3, [2,3]) == ansTensor(X, 3, [2,3], Y, 3, [2,3])
    disp('test case 17 : succeed');
else
    disp('test case 17 : fail');
end

end



function [X,numindX]=ansTensor(X,numindX,indX,Y,numindY,indY)
% Contraction of index indX of tensor X with index indY of tensor Y (X
% and Y have a number of indices corresponding to numindX and numindY
% respectively)
% indXとindYがスカラーの場合だけでなく、配列の場合も想定している。
Xsize=ones(1,numindX); Xsize(1:length(size(X)))=size(X);
Ysize=ones(1,numindY); Ysize(1:length(size(Y)))=size(Y);
indXl=1:numindX; indXl(indX)=[]; % indXlはindX以外の添え字
indYr=1:numindY; indYr(indY)=[]; 
sizeXl=Xsize(indXl);
sizeX=Xsize(indX);
sizeYr=Ysize(indYr);
sizeY=Ysize(indY);
if prod(sizeX)~=prod(sizeY) %prodは要らないんじゃないか？
    error('indX and indY are not of same dimension.');
end
if isempty(indYr)
    if isempty(indXl)
        X=permute(X,[indX]);
        X=reshape(X,[1,prod(sizeX)]);
        Y=permute(Y,[indY]);
        Y=reshape(Y,[prod(sizeY),1]);
        X=X*Y;
        Xsize=1;
        return;
    else
        X=permute(X,[indXl,indX]);
        X=reshape(X,[prod(sizeXl),prod(sizeX)]);
        Y=permute(Y,[indY]);
        Y=reshape(Y,[prod(sizeY),1]);
        X=X*Y;
        Xsize=Xsize(indXl);
        X=reshape(X,[Xsize,1]);
        return
    end
end
X=permute(X,[indXl,indX]);
X=reshape(X,[prod(sizeXl),prod(sizeX)]);
Y=permute(Y,[indY,indYr]);
Y=reshape(Y,[prod(sizeY),prod(sizeYr)]);
X=X*Y;
Xsize=[Xsize(indXl),Ysize(indYr)];
numindX=length(Xsize);
X=reshape(X,[Xsize,1]);
end

