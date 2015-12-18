function test_rightNormalize
% デバッグの時はrightNormalized()の一番最後で規格化しない方がよい。
N=2; % N=2で固定
D=5;
d=4;
mps = cell(1, N);
mps{1} = randn(1,D,d);
mps{N} = randn(D,1,d);
for i = 2:length(mps)-1
    mps{i} = randn(D,D,d);
end
Bs = rightNormalize( mps );

for i = 1:N
    temp = Bs{i}(:,:,1) * Bs{i}(:,:,1)';
    for sigma = 2:d
        temp = Bs{i}(:,:,sigma) * Bs{i}(:,:,sigma)' + temp;
    end
    disp(temp);
end

for dd = 1:d
    amplitude_original = 1;
    amplitude_canonical = 1;
    for n =1:N
        amplitude_original = amplitude_original * mps{n}(:,:,dd);
        amplitude_canonical = amplitude_canonical * Bs{n}(:,:,dd);
    end
    disp(amplitude_original - amplitude_canonical);
end

amplitude_canonical = 0;
for d1=1:d
    for d2=1:d
        amplitude_canonical = ( Bs{1}(:,:,d1) * Bs{2}(:,:,d2) )^2 + amplitude_canonical;
    end
end
disp(amplitude_canonical)

end


