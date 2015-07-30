% Post processing on averaged images

clc
close all
clear

addpath '../misc/'
addpath '../utils/'
addpath '../utils/lsqrSOL'
addpath '../utils/export_fig/'
addpath 'unlocbox/solver'

saveFigs = false;

% calibration sequence parameters
calibrationSequence = 'hadamardSeparable'; % hadamardSeparable
screenPatchSize = 128;
switch calibrationSequence
    case 'hadamard'
        calibPM = 1;
        N = screenPatchSize^2;
    case 'hadamardSeparable'
        calibPM = 1;
        N = screenPatchSize*2;
end
maskType = '255M030';
switch maskType
    case {'16Mseq','08Mseq','031M100','127M060','255M030'}
        op_makeSeparable = @(z) subtract_mean(z);
    case {'16Sseq'}
        op_makeSeparable = @(z) z;
end 
N = screenPatchSize^2;
%
downSamplingFactor = 1;
downSamplingMode = 'subsample';

msg1 = sprintf('computeCalib separable: calibration-%s, screenPatchSize=%d, maskType=%s, downSampling=%02d (%s)',...
    calibrationSequence, screenPatchSize,maskType,downSamplingFactor,downSamplingMode);
fprintf([msg1,'\n']);

calibDir = sprintf('../Data/mask%s_%s%d/',maskType,calibrationSequence,screenPatchSize);

%%
calibCompute_mat = sprintf('%sCalibCompute_downSampled%02d_%s.mat',calibDir,downSamplingFactor,downSamplingMode);
if ~exist(calibCompute_mat)
    %% Read stored images
    load(sprintf('%sCalibData_downSampled%02d_%s.mat',calibDir,downSamplingFactor,downSamplingMode));
    calibDir = sprintf('../Data/mask%s_%s%d/',maskType,calibrationSequence,screenPatchSize);
    
    fprintf('reading images ... ');
    % calibration images
    sensorImageP = double(sensorImageP);
    sensorImageN = double(sensorImageN);
    sensorZero = double(sensorZero);
    sensorImageP = (sensorImageP-repmat((mean(sensorZero,2)),1,size(sensorImageP,2)));
    sensorImageN = (sensorImageN-repmat((mean(sensorZero,2)),1,size(sensorImageN,2)));
    %
    % remove hot/dark pixels
    remove_pixels = [];
    % sensorImageP(remove_pixels,:) = [];
    % sensorImageN(remove_pixels,:) = [];
    % script_removeFaultyPixels
    %
    % scene_h = [screenImageP screenImageN];
    % sensor_h = [sensorImageP sensorImageN];
    %
    % rearrange hadamard sequence in the order or increasing frequency...
    % script_convertHadamardIndices;
    %
    % read screen image from the pattern file..
    calibMatrix = imread(sprintf('../misc/%s%d.bmp',calibrationSequence,screenPatchSize));
    screenImage = double(calibMatrix)/double(max(calibMatrix(:)));
    %
    screenImageP = screenImage;
    % screenImageP(screenImageP<0) = 0;
    screenImageN = -screenImage+1;
    % screenImageN(screenImageN<0) = 0;
    screenImage = [];
    
    fprintf('done! \n');
    
    %% compute inverse matrix ..
    invType = 'pm'; 
    
    fprintf('computing calibration matrix using %s ... \n ',invType);
    sensorImage_left = sensorImageP(:,1:end/2)-sensorImageN(:,1:end/2);
    sensorImage_right = sensorImageP(:,end/2+1:end)-sensorImageN(:,end/2+1:end);
    screenImage_left = screenImageP(:,1:end/2)-screenImageN(:,1:end/2); % e.g., hadamard pattern repeated in columns (Im = H(:,i)*ones(1,n))
    screenImage_right = screenImageP(:,end/2+1:end)-screenImageN(:,end/2+1:end); % e.g., hadamard pattern repeated in rows
    sensorImageP = []; sensorImageN = [];
    screenImageP = []; screenImageN = [];
    sensorZero = [];
    
    % different ways to calculate separable components
    % rank = 4; op_makeSeparable = @(z) z; % rank-2 without mean subtraction
    rank = 1; op_makeSeparable = @(z) subtract_mean(z); % rank-1 after mean substraction..
    
    % the "right" way to compute calibration matrix using svd...
    % fprintf('compute joint svd of all the left images for left matrix and all the right images for the right matrix\n');
    ItI_left = 0; ItI_right = 0;
    Itl_stack = zeros(sensorSize(1)*screenPatchSize,sensorSize(2)); 
    Itr_stack = zeros(sensorSize(2)*screenPatchSize,sensorSize(1)); 
    for ii =1:screenPatchSize;
        It = reshape(sensorImage_left(:,ii),sensorSize);
        Itl = op_makeSeparable(It);
        ItI_left = ItI_left + Itl'*Itl;
        Itl_stack((ii-1)*sensorSize(1)+1:ii*sensorSize(1),:) = Itl;
        
        It = reshape(sensorImage_right(:,ii),sensorSize)';
        Itr = op_makeSeparable(It);
        ItI_right = ItI_right + Itr'*Itr;
        Itr_stack((ii-1)*sensorSize(2)+1:ii*sensorSize(2),:) = Itr;
    end
    [Vl Sl Vl] = svd(ItI_left); Sl = diag(Sl);
    [Vr Sr Vr] = svd(ItI_right); Sr = diag(Sr);
    Itl_Vl = Itl_stack*Vl;
    Itr_Vr = Itr_stack*Vr;
    Itl_stack = []; Itr_stack = [];
    
    Phi_left = reshape(Itl_Vl(:,1:rank),sensorSize(1),[],rank);
    Phi_right = reshape(Itr_Vr(:,1:rank),sensorSize(2),[],rank);
    
    % test rank..
    %         for ii=1:screenPatchSize;
    %             r=1;
    %             figure(1); subplot(211);
    %             imagesc([Itl_stack((ii-1)*sensorSize(1)+1:ii*sensorSize(1),:)  Itl_Vl((ii-1)*sensorSize(1)+1:ii*sensorSize(1),1:r)*Vl(:,1:r)']);
    %             title(ii);
    %             subplot(212);
    %             imagesc([Itr_stack((ii-1)*sensorSize(2)+1:ii*sensorSize(2),:) Itr_Vr((ii-1)*sensorSize(2)+1:ii*sensorSize(2),1:r)*Vr(:,1:r)']);
    %             shg; drawnow; pause;
    %         end
    Itr_Vr = []; Itl_Vl = [];
    
    Phi_rec_left = zeros(sensorSize(1),screenPatchSize,rank);
    Phi_rec_right = zeros(sensorSize(2),screenPatchSize,rank);
    for r=1:rank
        Phi_rec_left(:,:,r) = Phi_left(:,:,r)*(screenImage_left(1:screenPatchSize,:)')/norm((screenImage_left(:,1)))^2;
        Phi_rec_right(:,:,r) = Phi_right(:,:,r)*(screenImage_left(1:screenPatchSize,:)')/norm((screenImage_left(:,1)))^2;
    end
    for r = 1:rank
        sum_Phi_rec_left = sum(Phi_rec_left(:,:,r),2);
        sum_Phi_rec_right = sum(Phi_rec_right(:,:,r),2);
        flip_sign =  sign(sum_Phi_rec_left'*Vr(:,r));
        figure(r);
        subplot(211); plot([sum_Phi_rec_left/norm(sum_Phi_rec_left) Vr(:,r)/norm(Vr(:,r))])
        subplot(212); plot([sum_Phi_rec_right/norm(sum_Phi_rec_right) Vl(:,r)/norm(Vl(:,r))])
        Phi_rec_right(:,:,r) = flip_sign*Phi_rec_right(:,:,r);
    end
    
    % normalize the matrices...
    tst = 1;%randint(1,1,[1 screenPatchSize])
    xTemp = reshape(screenImage_left(:,tst),screenPatchSize,[]);
    It = reshape(sensorImage_left(:,tst),sensorSize);
    It = op_makeSeparable(It);
    yTemp = 0;
    for r = 1:rank
        yTemp = yTemp+ Phi_rec_left(:,:,r)*xTemp*Phi_rec_right(:,:,r)';
    end
    normalizationFactor = sqrt(norm(It(:))/norm(yTemp(:)));
    Phi_rec_left = Phi_rec_left*normalizationFactor;
    Phi_rec_right = Phi_rec_right*normalizationFactor;
    
    sensorImage_left = []; sensorImage_right = [];
    screenImage_left = []; screenImage_right = [];
    
    save(sprintf('%sCalibCompute_downSampled%02d_%s.mat',calibDir,downSamplingFactor,downSamplingMode),'-v7.3');
    fprintf('done! \n');
end

%% test real images
vec = @(z) z(:); 
tmp_mat = matfile(sprintf('%sCalibCompute_downSampled%02d_%s.mat',calibDir,downSamplingFactor,downSamplingMode));
sensorSize = tmp_mat.sensorSize;
Phi_rec_left = tmp_mat.Phi_rec_left;
Phi_rec_right = tmp_mat.Phi_rec_right;

invMode1 = 'svd';
tvReg = true; 
fprintf('computing inverse matrix using %s ...  ',invMode1);
switch invMode1
    case 'svd'
        % svd for left- and right-side matrices
        [Ul Sl Vl] = svd(Phi_rec_left,'econ'); Sl = diag(Sl);
        [Ur Sr Vr] = svd(Phi_rec_right','econ'); Sr = diag(Sr);
        
        T = round(screenPatchSize/5); lambda = 2;
        iSl = speye(screenPatchSize); iSl = spdiags(1./(Sl+lambda),0,iSl);
        iSr = speye(screenPatchSize); iSr = spdiags(1./(Sr+lambda),0,iSr);
        iPhi_rec_left = @(z) Vl(:,1:T)*iSl(1:T,1:T)*Ul(:,1:T)'*reshape(z,sensorSize(1),[]);
        iPhi_rec_right = @(z) reshape(z,[],sensorSize(2))*Vr(:,1:T)*iSr(1:T,1:T)*Ur(:,1:T)';

        Phi_rec_left_SVD = Ul(:,1:T)*diag(Sl(1:T))*Vl(:,1:T)';
        Phi_rec_right_SVD = Ur(:,1:T)*diag(Sr(1:T))*Vr(:,1:T)';
        op_Phi_recSVD_TV = @(z) Phi_rec_left_SVD*z*Phi_rec_right_SVD;
        opt_Phi_recSVD_TV = @(z) Phi_rec_left_SVD'*z*Phi_rec_right_SVD';
        
        op_Phi_rec_TV = @(z) Phi_rec_left*z*Phi_rec_right';
        opt_Phi_rec_TV = @(z) Phi_rec_left'*z*Phi_rec_right;
        
        iPhi_rec = @(z) vec(iPhi_rec_right(iPhi_rec_left(z)));

        
        op_Phi_rec = @(z) vec(Phi_rec_left*(reshape(z,screenPatchSize,screenPatchSize)*Phi_rec_right'));
        
        
    case 'lsqr'
        % lambda = sqrt(1);
        % f_h = @(z) op_separableSystem(z, Phi_rec_left, Phi_rec_right, screenPatchSize,sensorSize,lambda);
        % ft_h = @(z) adj_separableSystem(z, Phi_rec_left, Phi_rec_right, screenPatchSize,sensorSize,lambda);
        %
        % fhandle = @(z,mode) fhandle_mode(z,mode, f_h, ft_h);
        % tol = 1e-6;
        % maxit = 1000;
        % iPhi_rec = @(z) lsqr(fhandle,[z;zeros(screenPatchSize^2*(lambda>0),1)],tol,maxit);
        
        op_Phi_rec = @(z) op_separableSystem(z, Phi_rec_left, Phi_rec_right, screenPatchSize,sensorSize);
        
        damp = 2;
        atol = 1e-6;
        btol = 1e-6;
        conlim = 1e6;
        itnlim = 1000;
        show = 0;
        f_h = @(z) op_separableSystem(z, Phi_rec_left, Phi_rec_right, screenPatchSize,sensorSize);
        ft_h = @(z) adj_separableSystem(z, Phi_rec_left, Phi_rec_right, screenPatchSize,sensorSize);
        fhandle = @(z,mode) fhandle_mode(z,mode, f_h, ft_h);
        iPhi_rec = @(z) lsqrSOL(numel(z), screenPatchSize^2, fhandle, z, damp, atol, btol, conlim, itnlim, show );
end 

% test images
testDir = sprintf('../Data/mask%s_realImages/',maskType);
tmp_mat = matfile(sprintf('%srealImages_downSampled%02d_%s.mat',testDir,downSamplingFactor,downSamplingMode));
sensorTest = double(tmp_mat.sensorTest);

fprintf('estimating real images ... ');

% test images
if saveFigs
    test_indices = find(sum(abs(sensorTest),1)>0);
else
    test_indices = [1];
%     test_indices = find(sum(abs(sensorTest),1)>0);
end

sensor_test = sensorTest(:,test_indices);
%
% remove faulty pixels
remove_pixels = [];
sensor_test(remove_pixels,:) = [];

% simulated measurements and reconstruction for test images
scene_test_rec = zeros(screenPatchSize^2,length(test_indices));
scene_test_recT = zeros(screenPatchSize^2,length(test_indices));
for tst = 1:length(test_indices)
    It = reshape(sensor_test(:,tst),sensorSize);
    if strcmpi(calibrationSequence,'hadamardSeparable')
        It = op_makeSeparable(It);
        sensor_test(:,tst) = It(:);
    end
    scene_test_rec(:,tst) = iPhi_rec(sensor_test(:,tst));
    if(tvReg)
        epsilon = 3000;
        param.maxit = 6;
        scene_test_recTV(:,tst) = vec(solve_tvdn(reshape(sensor_test(:,tst),sensorSize), epsilon, op_Phi_rec_TV, opt_Phi_rec_TV,param));
    end
end
% scene_test_rec = bsxfun(@rdivide, scene_test_rec, max(abs(scene_test_rec),[],1));

f881=figure(881);
if(tvReg)
    subplot(2,1,2);
    imagesc(reshape(scene_test_recTV,screenPatchSize,[]),[0 1]);
    title('+TV Regularization');
    colormap gray; axis image; axis off;
    subplot(2,1,1);
end
imagesc(reshape(scene_test_rec,screenPatchSize,[]),[0 1]);
title(sprintf([msg1, '\n real test images']));
colormap gray; axis image; axis off;
%
% maximize(f881);
set(gcf,'Color','w');
if saveFigs
    figName = sprintf('recovery_calib%s_%dx%d_%s_DS%02d%s_realImages', ...
        calibrationSequence, screenPatchSize, screenPatchSize,maskType,downSamplingFactor,downSamplingMode);
    export_fig(f881,[testDir,figName],'-png','-pdf','-native');
    saveas(f881,[testDir,figName],'fig');
end
fprintf('done! \n');
