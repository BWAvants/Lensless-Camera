clc
% close all
% clear
clearvars -except f101 

addpath '../misc/'
addpath '../utils/'
addpath '../utils/lsqrSOL'
addpath '../utils/export_fig/'
addpath './fastacodes'
addpath(genpath('unlocbox'))
addpath(genpath('TVAL3_beta2.4'))
addpath(genpath('WAVELAB850'))

global Phi_rec_left Phi_rec_right sensorSize screenPatchSize waveletLevel QMF

targetImages = 'testVideos';
saveFigs = true;  %test_indices = [] ;%[1 3 5] [7 149 159]

screenPatchSize = 256;
numOfChannels = 1;
downSamplingFactor = 4;
downSamplingMode = 'average';
%%
invMode = {'svd','lsqr'};

%Wavelet
waveletType = 'Daubechies';
waveletPar = 8;
waveletLevel = 5;
QMF = MakeONFilter(waveletType,waveletPar);

%FastaLSQR
optsFastaLSQR = [];
optsFastaLSQR.tol = 1e-5;  % Use super strict tolerance
optsFastaLSQR.recordObjective = true; %  Record the objective function so we can plot it
optsFastaLSQR.verbose = true;
optsFastaLSQR.maxIters = 20;
muFastaLSQR = 3;

%FastaLasso
optsFastaLasso = [];
optsFastaLasso.tol = 1e-5;  % Use super strict tolerance
optsFastaLasso.recordObjective = true; %  Record the objective function so we can plot it
optsFastaLasso.verbose = true;
optsFastaLasso.maxIters = 20;
muFastaLasso = 3400;

%Wave+TVAL3
optsTVAL3Wave.TVL2 = false;
% optsTVAL3Wave.nonneg = true;
optsTVAL3Wave.mu = 1000;
optsTVAL3Wave.beta = 100000;
optsTVAL3Wave.disp = true;
optsTVAL3Wave.maxit = 100;
optsTVAL3Wave.TVnorm = 1;

%TV
epsilonTV = 50000;
paramTV.maxit = 5;
% paramTV.tight_b2 = 0;
% paramTV.tol_b2 = 0.5;
% paramTV.nu_b2 = 3; 

%%TVAL3
optsTVAL3.TVL2 = true;
% optsTVAL3.nonneg = true;
optsTVAL3.mu = 1000;
optsTVAL3.beta = 100000;
optsTVAL3.disp = true;
optsTVAL3.maxit = 100;
optsTVAL3.isreal = true;
optsTVAL3.TVnorm = 1;
% optsTVAL3.scale_A = false;
% optsTVAL3.scale_b = optsTVAL3.scale_A;
%% calibration sequence parameters
vec = @(z) z(:); 
calibrationSequence = 'hadamardSeparable'; % hadamardSeparable
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

msg1 = sprintf('computeCalib separable: calibration-%s, screenPatchSize=%d, maskType=%s, downSampling=%02d (%s)\n',...
    calibrationSequence, screenPatchSize,maskType,downSamplingFactor,downSamplingMode);
fprintf([msg1,'\n']);

calibDir = sprintf('../Data/mask%s_%s%d/',maskType,calibrationSequence,screenPatchSize);

%% load hadamard
tmp_mat = matfile(sprintf('%sCalibCompute_downSampled%02d_%s.mat',calibDir,downSamplingFactor,downSamplingMode));
sensorSize = tmp_mat.sensorSize;
Phi_rec_left = tmp_mat.Phi_rec_left;
Phi_rec_right = tmp_mat.Phi_rec_right;

for i = 1:length(invMode)
    fprintf('computing inverse matrix using %s ...  \n',invMode{i});
    switch invMode{i}
        case {'svd', 'svdTV','TV-unlocbox','FTVd'}
            % svd for left- and right-side matrices
            [Ul,Sl,Vl] = svd(Phi_rec_left,'econ'); Sl = diag(Sl);
            [Ur,Sr,Vr] = svd(Phi_rec_right','econ'); Sr = diag(Sr);
            switch invMode{i}
                case 'svd'    
                    T = round(screenPatchSize/2); lambda = 2;
                    iSl = speye(screenPatchSize); iSl = spdiags(1./(Sl+lambda),0,iSl);
                    iSr = speye(screenPatchSize); iSr = spdiags(1./(Sr+lambda),0,iSr);
                    iPhi_rec_left = @(z) Vl(:,1:T)*iSl(1:T,1:T)*Ul(:,1:T)'*reshape(z,sensorSize(1),[]);
                    iPhi_rec_right = @(z) reshape(z,[],sensorSize(2))*Vr(:,1:T)*iSr(1:T,1:T)*Ur(:,1:T)';
                    iPhi_rec_svd = @(z) vec(iPhi_rec_right(iPhi_rec_left(z)));        
                case 'svdTV'
                    Phi_rec_left_SVD = Ul(:,1:T)*diag(Sl(1:T))*Vl(:,1:T)';
                    Phi_rec_right_SVD = Ur(:,1:T)*diag(Sr(1:T))*Vr(:,1:T)';
                    op_Phi_rec_svdTV = @(z) Phi_rec_left_SVD*z*Phi_rec_right_SVD;
                    opt_Phi_rec_svdTV = @(z) Phi_rec_left_SVD'*z*Phi_rec_right_SVD';
                case {'TV-unlocbox', 'FTVd'}
                    T = round(screenPatchSize/1); lambda = 2;
                    iSl_TV = speye(screenPatchSize); iSl_TV = spdiags(1./(Sl+lambda),0,iSl_TV);
                    iSr_TV = speye(screenPatchSize); iSr_TV = spdiags(1./(Sr+lambda),0,iSr_TV);
                    iPhi_rec_left_TV = @(z) Vl(:,1:T)*iSl_TV(1:T,1:T)*Ul(:,1:T)'*z;
                    iPhi_rec_right_TV = @(z) z*Vr(:,1:T)*iSr_TV(1:T,1:T)*Ur(:,1:T)';

                    op_Phi_rec_TV = @(z) Phi_rec_left*z*Phi_rec_right';
%                     opt_Phi_rec_TV = @(z) Phi_rec_left'*z*Phi_rec_right;
                    opt_Phi_rec_TV = @(z) iPhi_rec_right_TV(iPhi_rec_left_TV(z));
                
            end
        case {'Wave+TVAL3','Wave+FastaLasso','Wave+FastaLSQR'}
            op_Phi_rec_WaveF = @(z) fun_Phi_rec_Wave(z,1);
            op_Phi_rec_WaveB = @(z) fun_Phi_rec_Wave(z,2);
            op_Phi_rec_Wave = @(z,mode) fun_Phi_rec_Wave(z, mode);
        case {'DCT+FastaLasso','DCT+FastaLSQR'}
            op_Phi_rec_DCTF = @(z) fun_Phi_rec_DCT(z,1);
            op_Phi_rec_DCTB = @(z) fun_Phi_rec_DCT(z,2);
            op_Phi_rec_DCT = @(z,mode) fun_Phi_rec_DCT(z, mode);
        case 'TVAL3'
            op_Phi_rec_TVAL3 = @(z,mode) fun_Phi_rec_TVAL3(z, mode);
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
            iPhi_rec_lsqr = @(z) lsqrSOL(numel(z), screenPatchSize^2, fhandle, z, damp, atol, btol, conlim, itnlim, show );
    end 
end
%% test images
switch targetImages
    case 'realImages'
        if numOfChannels > 1
            testDir = sprintf('../Data/mask%s_realImagesC/',maskType);
        else 
            testDir = sprintf('../Data/mask%s_realImages/',maskType);
        end
        tmp_mat = matfile(sprintf('%srealImages_downSampled%02d_%s.mat',testDir,downSamplingFactor,downSamplingMode));
        sensorTest = double(tmp_mat.sensorTest);
    case {'testImages', 'simulatedImages'}
        testDir = sprintf('../Data/mask%s_%dx%d_testImages/',maskType,screenPatchSize,screenPatchSize);
        tmp_mat = matfile(sprintf('%stestImages_downSampled%02d_%s.mat',testDir,downSamplingFactor,downSamplingMode));
        sceneTest = double(imread(sprintf('../misc/lowresTest%d.bmp',screenPatchSize)));
        sceneTest = reshape([sceneTest; max(sceneTest(:))-sceneTest],size(sceneTest,1),[]);
        sensorTest = double(tmp_mat.sensorTest_lowres);
    case 'testVideos'
        testDir = sprintf('../Data/mask%s_testVideos/',maskType);
        tmp_mat = matfile(sprintf('%stestVideos_downSampled%02d_%s.mat',testDir,downSamplingFactor,downSamplingMode));
%         sceneTest = double(imread(sprintf('../misc/lowresTest%d.bmp',screenPatchSize)));
%         sceneTest = reshape([sceneTest; max(sceneTest(:))-sceneTest],size(sceneTest,1),[]);
        sensorTest = double(cell2mat(tmp_mat.sensorTest'));  
        videoFrames = tmp_mat.videoFrames;
end


fprintf('estimating %s ... \n',targetImages);

%% test images
if ~exist('test_indices','var')
    test_indices = find(sum(abs(sensorTest),1)>0);
end

if strcmp(targetImages,'simulatedImages');
    sensor_test = zeros(prod(sensorSize),length(test_indices));
    for i = 1:length(test_indices)
        sensor_test(:,i) = vec(Phi_rec_left*reshape(sceneTest(:,test_indices(i)),screenPatchSize,[])*Phi_rec_right');
    end
else
    sensor_test = sensorTest(:,test_indices);
end

if exist('sceneTest','var')
    scene_test = sceneTest(:,test_indices);
end
%
% remove faulty pixels
remove_pixels = [];
sensor_test(remove_pixels,:) = [];

%% simulated measurements and reconstruction for test images
scene_test_rec = zeros(screenPatchSize^2,length(test_indices),length(invMode));
for tst = 1:length(test_indices)
    for i = 1:length(invMode)
        switch invMode{i}
            case 'svd'        
                scene_test_rec(:,tst,i) = iPhi_rec_svd(sensor_test(:,tst));
            case 'svdTV'
                scene_test_rec(:,tst,i) = vec(solve_tvdn(reshape(sensor_test(:,tst),sensorSize), epsilonTV, op_Phi_rec_svdTV, opt_Phi_rec_svdTV,paramTV));
            case 'TV-unlocbox'
                scene_test_rec(:,tst,i) = vec(solve_tvdn(reshape(sensor_test(:,tst),sensorSize), epsilonTV, op_Phi_rec_TV, opt_Phi_rec_TV,paramTV));
            case 'FTVd'
                scene_test_rec(:,tst,i) = FTVd_v4(op_Phi_rec_TV,reshape(sensor_test(:,tst),sensorSize),1,'L2');
            case 'TVAL3'
%                 optsTVAL3.init = scene_test_rec(:,tst,i-2);
                scene_test_rec(:,tst,i) = vec(TVAL3(op_Phi_rec_TVAL3,sensor_test(:,tst),screenPatchSize,screenPatchSize,optsTVAL3));
            case 'Wave+TVAL3'
                scene_test_rec(:,tst,i) = vec(IWT2_PO(TVAL3(op_Phi_rec_Wave,sensor_test(:,tst),screenPatchSize,screenPatchSize,optsTVAL3Wave),waveletLevel,QMF));
            case 'Wave+FastaLasso'
                scene_test_rec(:,tst,i) = vec(IWT2_PO(reshape(fasta_lasso(op_Phi_rec_WaveF,op_Phi_rec_WaveB,sensor_test(:,tst),muFastaLasso, zeros(screenPatchSize^2,1), optsFastaLasso),screenPatchSize,screenPatchSize),waveletLevel,QMF));
            case 'Wave+FastaLSQR'
                scene_test_rec(:,tst,i) = vec(IWT2_PO(reshape(fasta_sparseLeastSquares(op_Phi_rec_WaveF,op_Phi_rec_WaveB,sensor_test(:,tst),muFastaLSQR, zeros(screenPatchSize^2,1), optsFastaLSQR),screenPatchSize,screenPatchSize),waveletLevel,QMF));
            case 'DCT+FastaLasso'
                scene_test_rec(:,tst,i) = vec(idct2blk(reshape(fasta_lasso(op_Phi_rec_DCTF,op_Phi_rec_DCTB,sensor_test(:,tst),muFastaLasso, zeros(screenPatchSize^2,1), optsFastaLasso),screenPatchSize,screenPatchSize)));
            case 'DCT+FastaLSQR'
                scene_test_rec(:,tst,i) = vec(idct2blk(reshape(fasta_sparseLeastSquares(op_Phi_rec_DCTF,op_Phi_rec_DCTB,sensor_test(:,tst),muFastaLSQR, zeros(screenPatchSize^2,1), optsFastaLSQR),screenPatchSize,screenPatchSize)));
            case 'lsqr'
                scene_test_rec(:,tst,i) = iPhi_rec_lsqr(sensor_test(:,tst));
        end 
    end
end
%% 
switch targetImages
    case {'testImages', 'realImages', 'simulatedImages'}
        if ~saveFigs
            if exist('f101','var');
                figure(f101);
            else
                f101 = figure;
            end
            k = 1;
            for tst = 1:length(test_indices)/numOfChannels
                for i = 1:length(invMode)
                    subplot(length(test_indices)/numOfChannels,length(invMode),k);
                    thisImage = scene_test_rec(:,(tst-1)*numOfChannels+1:tst*numOfChannels,i);
                    thisImageNormalized = (thisImage-min(thisImage(:)));
                    thisImageNormalized = reshape(thisImageNormalized/max(thisImageNormalized(:)),screenPatchSize,screenPatchSize,[]);
                    imagesc(thisImageNormalized);title(invMode{i});
                    colormap gray; axis image; axis off;
                    k = k + 1;
                end
            end
        else
            k = 1;
            for tst = 1:length(test_indices)/numOfChannels
                for i = 1:length(invMode)
                    thisImage = scene_test_rec(:,(tst-1)*numOfChannels+1:tst*numOfChannels,i);
                    thisImageNormalized = (thisImage-min(thisImage(:)));
                    thisImageNormalized = reshape(thisImageNormalized/max(thisImageNormalized(:)),screenPatchSize,screenPatchSize,[]);
                    imagesc(thisImageNormalized);title(invMode{i});
                    colormap gray; axis image; axis off;
                    pause(1/10)
                    imwrite(thisImageNormalized,['Results\' targetImages '\result_' num2str(test_indices(tst)) '_' invMode{i} '_' num2str(screenPatchSize) '.bmp']);
                    k = k + 1;
                end
            end
        end
    case 'testVideos'
%         thisImage = scene_test_rec;
%         thisImageNormalized = (thisImage-min(thisImage(:)));
%         thisImageNormalized = thisImageNormalized/max(thisImageNormalized(:));
%         imwrite(255*reshape(scene_test_rec(:,1:2:end,2),screenPatchSize,screenPatchSize,1,[]),'Results/testVideos/LenslessCamera.gif','DelayTime',1/20,'LoopCount',3);
%         sceneTestCell = mat2cell(scene_test_rec);%,size(scene_test_rec,1),length(test_indices),2);
end