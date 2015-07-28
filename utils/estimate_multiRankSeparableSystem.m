clear 

addpath(genpath('utils/'));
addpath 'misc/'

global m n T maxrank
vec = @(z) z(:);

% calibration sequence parameters
calibrationSequence = 'hadamardSeparable'; % hadamardSeparable
maskType = '127M060'; 
op_makeSeparable = @(z) subtract_mean(z);

downSamplingFactor = 16;
downSamplingMode = 'average';

n = 32; m = 64;
maxrank = 1; % rank of the measurement model
T = 20*(2*n*maxrank/m); % number of random images sampled

sensorSize = [m,m]; % sensor size
screenPatchSize = n; % image size

% convert 0/1 images to -1/1
screenImages_all = double(imread(sprintf('misc/random%d.bmp',n)));
screenImages_all(screenImages_all==0) = -1; 
screenImages_all = reshape(screenImages_all,n,n,[]);
screenImages = screenImages_all(:,:,1:T);

% read observed images for the random screen patterns
testDir = sprintf('Data/mask%s_%dx%d_testImages/',maskType,screenPatchSize,screenPatchSize);
test_mat = matfile(sprintf('%stestImages_downSampled%02d_%s.mat',testDir,downSamplingFactor,downSamplingMode));
%
% subtract positive and negative images
sensorImages_all = test_mat.sensorTest_random;
sensorImages_all = reshape(sensorImages_all,sensorSize(1),2*sensorSize(2),[]);
sensorImages_all = sensorImages_all(:,1:sensorSize(2),:)-sensorImages_all(:,sensorSize(2)+1:end,:); 
for t = 1:size(sensorImages_all,3)
    sensorImages_all(:,:,t) = op_makeSeparable(sensorImages_all(:,:,t));
end 
sensorImages = sensorImages_all(:,:,1:T); 

% Ground truth using hadamardSeparable calibration
calibDir = sprintf('Data/mask%s_%s%d/',maskType,calibrationSequence,screenPatchSize);
calib_mat = matfile(sprintf('%sCalibCompute_downSampled%02d_%s.mat',calibDir,downSamplingFactor,downSamplingMode));
A_cube = calib_mat.Phi_rec_left(:,:,1:maxrank);
B_cube = calib_mat.Phi_rec_right(:,:,1:maxrank);
% verify
simulatedImages1 = 0*sensorImages;
for t = 1:T
    simulatedImages1(:,:,t) = reshape(op_separableSystem(screenImages(:,:,t), A_cube, B_cube,n,sensorSize),m,m);
    figure(1); imagesc([sensorImages(:,:,t) simulatedImages1(:,:,t) sensorImages(:,:,t)-simulatedImages1(:,:,t)]); shg;
    pause(1/60); 
end

%%

LInit = randn(m*n,maxrank); 
RInit = randn(m*n,maxrank);

[Lr,Rr] = recover_multiRankSeparableSystem(sensorImages(:),screenImages(:),LInit,RInit);

%% post processing Ar_cube = reshape(Lr,m,n,maxrank);
Ar_cube = reshape(Lr,m,n,maxrank)*(sqrt(norm(Lr(:))*norm(Rr(:)))/norm(Lr(:)));
Br_cube = reshape(Rr,m,n,maxrank)*(sqrt(norm(Lr(:))*norm(Rr(:)))/norm(Rr(:)));

tst = 10; 
tst_range = randi([size(sensorImages,3)+1 size(sensorImages_all,3)],1,tst);
% tst_range = [T T+1]
testImages = screenImages_all(:,:,tst_range);
measuredImages = sensorImages_all(:,:,tst_range);

simulatedImages1 = zeros(m,m,tst);
simulatedImages2 = zeros(m,m,tst);
for t = 1:tst
    simulatedImages1(:,:,t) = reshape(op_separableSystem(testImages(:,:,t), A_cube, B_cube,n,sensorSize),m,m);
    simulatedImages2(:,:,t) = reshape(op_separableSystem(testImages(:,:,t), Ar_cube, Br_cube,n,sensorSize),m,m);
end

figure(101); imagesc([reshape(measuredImages,m,[]); reshape(simulatedImages1,m,[]); reshape(measuredImages-simulatedImages1,m,[])*1]); title('original and simultated-Calibrate measurements and error')
figure(102); imagesc([reshape(measuredImages,m,[]); reshape(simulatedImages2,m,[]); reshape(measuredImages-simulatedImages2,m,[])*1]); title('original and simultated-Random measurements and error')
figure(2); imagesc([reshape(A_cube,m,[]); reshape(Ar_cube,m,[]); reshape(A_cube-Ar_cube,m,[])*1]); title('original and reconstructed left matrices and errors');
figure(3); imagesc([reshape(B_cube,m,[]); reshape(Br_cube,m,[]); reshape(B_cube-Br_cube,m,[])*1]); title('original and reconstructed left matrices and errors');

for t=1:tst; figure(103); plot([vec(abs(measuredImages(:,:,t)-simulatedImages1(:,:,t))) vec(abs(measuredImages(:,:,t)-simulatedImages2(:,:,t)))]); shg; pause(1/60);  end

% reconstruct test images using the calibrated/estimated matrices... 
damp = 2; 
atol = 1e-6;
btol = 1e-6;
conlim = 1e6;
itnlim = 1000;
show = 0;

fprintf('estimating test images ... ');
% calib matrices
f_h = @(z) op_separableSystem(z, A_cube, B_cube,screenPatchSize,sensorSize);
ft_h = @(z) adj_separableSystem(z, A_cube, B_cube,screenPatchSize,sensorSize);
fhandle = @(z,mode) fhandle_mode(z,mode, f_h, ft_h);
iPhi_rec1 = @(z) lsqrSOL(numel(z), screenPatchSize^2, fhandle, z(:), damp, atol, btol, conlim, itnlim, show );
rec1 = zeros(n,n,tst);
for t = 1:tst 
    rec1(:,:,t) = reshape(iPhi_rec1(measuredImages(:,:,t)),n,n); 
end

% estimated matrices
f_h = @(z) op_separableSystem(z, Ar_cube, Br_cube,screenPatchSize,sensorSize);
ft_h = @(z) adj_separableSystem(z, Ar_cube, Br_cube,screenPatchSize,sensorSize);
fhandle = @(z,mode) fhandle_mode(z,mode, f_h, ft_h);
iPhi_rec2 = @(z) lsqrSOL(numel(z), screenPatchSize^2, fhandle, z(:), damp, atol, btol, conlim, itnlim, show );
rec2 = zeros(n,n,tst);
for t = 1:tst 
    rec2(:,:,t) = reshape(iPhi_rec2(measuredImages(:,:,t)),n,n); 
end

fprintf('done! \n');

figure(201); imagesc([reshape(testImages,n,[]); reshape(rec1,n,[]); reshape(testImages(:)-rec1(:),n,[])*1]); title('original images and their reconstruction using calibrated matrix and error')
figure(202); imagesc([reshape(testImages,n,[]); reshape(rec2,n,[]); reshape(testImages(:)-rec2(:),n,[])*1]); title('original images and their reconstruction using estimated matrix and error')
for t=1:tst; err1 = testImages(:,:,t)-rec1(:,:,t); err2 = testImages(:,:,t)-rec2(:,:,t); figure(203); plot([err1(:) err2(:)]); title(sprintf('error1 = %3.4g    error2 = %3.4g',norm(err1(:)), norm(err2(:)))); shg; pause;  end

