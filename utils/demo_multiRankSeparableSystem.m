clear 
rng(0)

addpath(genpath('utils/'));

global m n T maxrank
vec = @(z) z(:);

n = 32; m = 64;
maxrank = 2; % rank of the measurement model
T = 16*(2*n*maxrank/m); % number of random images sampled
matType = 'calib';
SNR = inf;

sensorSize = [m,m]; % sensor size
screenPatchSize = n; % image size

switch matType
    case 'random'
        A_cube = double(rand(m,n,maxrank)>.5);
        B_cube = double(rand(m,n,maxrank)>.5);
        weights = [1 0.5 0.25 0.1]>0;        
        for r = 1:maxrank
            A_cube(:,:,r) = A_cube(:,:,r)*weights(r);
            B_cube(:,:,r) = B_cube(:,:,r)*weights(r);
        end
    case 'calib'
        %
        calibrationSequence = 'hadamardSeparable'; % hadamardSeparable
        maskType = '127M060';
        downSamplingFactor = 16;
        downSamplingMode = 'average';
        calibDir = sprintf('Data/mask%s_%s%d/',maskType,calibrationSequence,screenPatchSize);
        calib_mat = matfile(sprintf('%sCalibCompute_downSampled%02d_%s.mat',calibDir,downSamplingFactor,downSamplingMode));
        A_cube = calib_mat.Phi_rec_left(:,:,1:maxrank);
        B_cube = calib_mat.Phi_rec_right(:,:,1:maxrank);
end

f_h = @(z) op_separableSystem(z, A_cube, B_cube,n,sensorSize);
ft_h = @(z) adj_separableSystem(z, A_cube, B_cube,n,sensorSize);

screenImages = zeros(n^2,T);
sensorImages = zeros(prod(sensorSize),T);
for t = 1:T
    X = sign(randn(n));
    screenImages(:,t) = X(:);
    sensorImages(:,t) = f_h(X(:));
end
sigma = norm(mean(sensorImages,2))/(m*10^(SNR/10));
noise = sigma*randn(prod(sensorSize),T); 
sensorImages_noisy = sensorImages+noise; 
%%

LInit = randn(m*n,maxrank); 
RInit = randn(m*n,maxrank);

[Lr,Rr] = recover_multiRankSeparableSystem(sensorImages(:),screenImages(:),LInit,RInit);

%% post processing 
Ar_cube = reshape(Lr,m,n,maxrank);
Br_cube = reshape(Rr,m,n,maxrank);

tst = 8;
testImages = randn(n,n,tst);
% testImages = reshape(screenImages,n,n,tst);
measuredImages = zeros(m^2,tst);
simulatedImages = zeros(m^2,tst);
for t = 1:tst
    x = testImages(:,:,t);
    measuredImages(:,t) = op_separableSystem(x(:), A_cube, B_cube,n,sensorSize);
    simulatedImages(:,t) = op_separableSystem(x(:), Ar_cube, Br_cube,n,sensorSize);
end

figure(1); imagesc([reshape(measuredImages,m,[]); reshape(simulatedImages,m,[]); reshape(measuredImages-simulatedImages,m,[])*1]); title('original and simultated measurements and error')
% figure(2); imagesc([reshape(A_cube,m,[]); reshape(Ar_cube,m,[]); reshape(A_cube-Ar_cube,m,[])*1]); title('original and reconstructed left matrices and errors');
% figure(3); imagesc([reshape(B_cube,m,[]); reshape(Br_cube,m,[]); reshape(B_cube-Br_cube,m,[])*1]); title('original and reconstructed left matrices and errors');

