% test separable multi-rank least-squares... 
clear

addpath 'lsqrSOL'

screenPatchSize = 256; 
sensorSize = [256 320];
rank = 4; 

X = imresize(double(imread('cameraman.tif')),[screenPatchSize screenPatchSize]);
A_cube = randn(sensorSize(1),screenPatchSize,rank);
B_cube = randn(sensorSize(2),screenPatchSize,rank);
weights = [1 0.1 0.09 0.05];
for r = 1:rank
    A_cube(:,:,r) = A_cube(:,:,r)*weights(r);
    B_cube(:,:,r) = B_cube(:,:,r)*weights(r);
end   
y = op_separableSystem(X(:), A_cube, B_cube,screenPatchSize,sensorSize);

%% matlab lsqr
lambda = 0; 
f_h = @(z) op_separableSystem(z, A_cube, B_cube,screenPatchSize,sensorSize,lambda);
ft_h = @(z) adj_separableSystem(z, A_cube, B_cube,screenPatchSize,sensorSize,lambda);
fhandle = @(z,mode) fhandle_mode(z,mode, f_h, ft_h);
tol = 1e-6;
maxit = 1000;
rec = lsqr(fhandle,[y;zeros(screenPatchSize^2*(lambda>0),1)],tol,maxit);

%% Saunders lsqr
damp = 0; 
atol = 1e-6;
btol = 1e-6;
conlim = 1e6;
itnlim = 1000;
show = 0;
f_h = @(z) op_separableSystem(z, A_cube(:,:,1), B_cube(:,:,1),screenPatchSize,sensorSize);
ft_h = @(z) adj_separableSystem(z, A_cube(:,:,1), B_cube(:,:,1),screenPatchSize,sensorSize);
fhandle = @(z,mode) fhandle_mode(z,mode, f_h, ft_h);
rec2 = lsqrSOL( numel(y), numel(X), fhandle, y, damp, atol, btol, conlim, itnlim, show );


%%
Xh = reshape(rec,screenPatchSize,screenPatchSize);
Xh2 = reshape(rec2,screenPatchSize,screenPatchSize);

figure(1); imagesc([X Xh Xh2])
title(sprintf('rank = %d, rel. rec. error1 = %3.4g, error2 = %3.4g',rank, norm(X(:)-Xh(:))/norm(X(:)), norm(X(:)-Xh2(:))/norm(X(:))));