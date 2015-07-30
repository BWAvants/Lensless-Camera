clc
close all
clear
clearvars -except f101 

addpath '../misc/'
addpath '../utils/'
addpath '../utils/lsqrSOL'
addpath '../utils/export_fig/'
addpath './fastacodes'
addpath(genpath('unlocbox'))
addpath(genpath('TVAL3_beta2.4'))
addpath(genpath('WAVELAB850'))
addpath(genpath('BM3D'))

global Phi_rec channelIndex sensorSize screenPatchSize waveletLevel QMF

targetImages = 'testImages';
saveFigs = true; test_indices = [13];% 125 127] ;%[1 3 5] [7 149 159]

screenPatchSize = 64;
numOfChannels = 1;
downSamplingFactor = 1;
downSamplingMode = 'subsample';
%%
BM3DFlag =true;
iterationFlag = 1;
invMode = {'Wave+FastaLasso'};
%'TVAL3', 'Wave+TVAL3','Wave+FastaLasso','Wave+FastaLSQR'

titleBM3D = '';

%Wavelet
waveletType = 'Daubechies';
waveletPar = 8;
waveletLevel = 3;
QMF = MakeONFilter(waveletType,waveletPar);

%FastaLSQR
%min  mu|x|+.5||Ax-b||
optsFastaLSQR = [];
optsFastaLSQR.tol = 1e-5;  % Use super strict tolerance
optsFastaLSQR.recordObjective = true; %  Record the objective function so we can plot it
optsFastaLSQR.verbose = true;
optsFastaLSQR.maxIters = 20;
muFastaLSQR = 3;

%FastaLasso
%min  .5||Ax-b||^2 s.t. |x|<mu
optsFastaLasso = [];
optsFastaLasso.tol = 1e-5;  % Use super strict tolerance
optsFastaLasso.recordObjective = true; %  Record the objective function so we can plot it
optsFastaLasso.verbose = true;
optsFastaLasso.maxIters = 20;
muFastaLasso = screenPatchSize^2/10;

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
%%min sum ||D_i u|| + mu/2||Au-b||_2^2
optsTVAL3.TVL2 = false;
optsTVAL3.nonneg = true;
optsTVAL3.mu = 100000;
optsTVAL3.beta = 100000;
optsTVAL3.disp = true;
optsTVAL3.maxit = 100;
optsTVAL3.isreal = true;
optsTVAL3.TVnorm = 2;
% optsTVAL3.scale_A = false;
% optsTVAL3.scale_b = optsTVAL3.scale_A;
%% calibration sequence parameters
vec = @(z) z(:); 
calibrationSequence = 'hadamard'; % hadamardSeparable
switch calibrationSequence
    case 'hadamard'
        calibPM = 1;
        N = screenPatchSize^2;
    case 'hadamardSeparable'
        calibPM = 1;
        N = screenPatchSize*2;
end
maskType = '64SWIR';
switch maskType
    case {'64SWIR','16Mseq','08Mseq','031M100','127M060','255M030','255M030C'}
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
Phi_rec = tmp_mat.Phi_rec;

for i = 1:length(invMode)
    fprintf('computing inverse matrix using %s ...  \n',invMode{i});
    switch invMode{i}
        case {'svd', 'svdTV','TV-unlocbox','FTVd'}
            iPhi_rec_matrix = [];
            U = zeros(size(Phi_rec,1),N,numOfChannels);
            V = zeros(N,N,numOfChannels);
            S = zeros(N,numOfChannels);
            % svd for left- and right-side matrices
            for C = 1:numOfChannels
                [U(:,:,C),Stemp,V(:,:,C)] = svd(Phi_rec(:,:,C),'econ'); S(:,C) = diag(Stemp);
            end
            switch invMode{i}
                case 'svd'
                    T = round(N/5); lambda = 0;
                    iS = zeros(N,N,numOfChannels);
                    
                    for C = 1:numOfChannels
                        iS(:,:,C) = diag(1./(S(:,C)+lambda));
                        iPhi_rec_matrix = [iPhi_rec_matrix; V(:,1:T,C)*iS(1:T,1:T,C)*U(:,1:T,C)'];
                    end
                    trim = @(z) z(:,(1:numOfChannels).*(1:numOfChannels));
                    iPhi_rec = @(z) iPhi_rec_matrix *z;
                    iPhi_rec_svd = @(z) iPhi_rec(z);        
%                 case 'svdTV'
%                     Phi_rec_left_SVD = Ul(:,1:T)*diag(Sl(1:T))*Vl(:,1:T)';
%                     Phi_rec_right_SVD = Ur(:,1:T)*diag(Sr(1:T))*Vr(:,1:T)';
%                     op_Phi_rec_svdTV = @(z) Phi_rec_left_SVD*z*Phi_rec_right_SVD;
%                     opt_Phi_rec_svdTV = @(z) Phi_rec_left_SVD'*z*Phi_rec_right_SVD';
%                 case {'TV-unlocbox', 'FTVd'}
%                     T = round(screenPatchSize/1); lambda = 2;
%                     iSl_TV = speye(screenPatchSize); iSl_TV = spdiags(1./(Sl+lambda),0,iSl_TV);
%                     iSr_TV = speye(screenPatchSize); iSr_TV = spdiags(1./(Sr+lambda),0,iSr_TV);
%                     iPhi_rec_left_TV = @(z) Vl(:,1:T)*iSl_TV(1:T,1:T)*Ul(:,1:T)'*z;
%                     iPhi_rec_right_TV = @(z) z*Vr(:,1:T)*iSr_TV(1:T,1:T)*Ur(:,1:T)';
% 
%                     op_Phi_rec_TV = @(z) Phi_rec_left*z*Phi_rec_right';
% %                     opt_Phi_rec_TV = @(z) Phi_rec_left'*z*Phi_rec_right;
%                     opt_Phi_rec_TV = @(z) iPhi_rec_right_TV(iPhi_rec_left_TV(z)); 
            end
        case {'Wave+TVAL3','Wave+FastaLasso','Wave+FastaLSQR'}
            op_Phi_rec_WaveF = @(z) fun_Phi_rec_Wave_SWIR(z,1);
            op_Phi_rec_WaveB = @(z) fun_Phi_rec_Wave_SWIR(z,2);
            op_Phi_rec_Wave = @(z,mode) fun_Phi_rec_Wave(z, mode);
        case {'DCT+FastaLasso','DCT+FastaLSQR'}
            op_Phi_rec_DCTF = @(z) fun_Phi_rec_DCT_SWIR(z,1);
            op_Phi_rec_DCTB = @(z) fun_Phi_rec_DCT_SWIR(z,2);
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
        testDir = calibDir;
        tmp_mat = matfile(sprintf('%stestImages_downSampled%02d_%s.mat',testDir,downSamplingFactor,downSamplingMode));
        sensorTest = double(tmp_mat.sensorTest_lowres);
        sensorTestZero = double(tmp_mat.sensorTest_zero);
%         sensorTestOne = double(tmp_mat.sensorTest_one);
        sceneTest = double(tmp_mat.lowresImages);
        sensorTest = sensorTest-repmat(mean(sensorTestZero,2),1,size(sensorTest,2));
        sensorTestMean = mean(sensorTest,1);
        sensorTest = (sensorTest - ones(prod(sensorSize),1)*sensorTestMean);
%         sensorTestOne = (sensorTestOne - ones(prod(sensorSize),1)*mean(sensorTestOne,1));
%         sensorTest(:,1:2:24) = sensorTest(:,1:2:24)-sensorTest(:,2:2:24); 
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
while(iterationFlag<=2)
scene_test_rec = zeros(screenPatchSize^2*numOfChannels,length(test_indices),length(invMode));
scene_test_rec_neg = zeros(screenPatchSize^2*numOfChannels,length(test_indices),length(invMode));
for tst = 1:length(test_indices)
    for i = 1:length(invMode)
        switch invMode{i}
            case 'svd'        
                scene_test_rec(:,tst,i) = iPhi_rec_svd(sensor_test(:,tst));
        end
        for channelIndex = 1:numOfChannels
            sceneChannelRange = screenPatchSize^2*(channelIndex-1)+1:screenPatchSize^2*channelIndex;
            sensorChannelRange = prod(sensorSize(1:2))*(channelIndex-1)+1:prod(sensorSize(1:2))*channelIndex;
            switch invMode{i}
                case 'svdTV'
                    scene_test_rec(:,tst,i) = vec(solve_tvdn(reshape(sensor_test(:,tst),sensorSize), epsilonTV, op_Phi_rec_svdTV, opt_Phi_rec_svdTV,paramTV));
                case 'TV-unlocbox'
                    scene_test_rec(:,tst,i) = vec(solve_tvdn(reshape(sensor_test(:,tst),sensorSize), epsilonTV, op_Phi_rec_TV, opt_Phi_rec_TV,paramTV));
                case 'FTVd'
                    scene_test_rec(:,tst,i) = FTVd_v4(op_Phi_rec_TV,reshape(sensor_test(:,tst),sensorSize),1,'L2');
                case 'TVAL3'
                    scene_test_rec(sceneChannelRange,tst,i) = vec(TVAL3(op_Phi_rec_TVAL3,sensor_test(sensorChannelRange,tst),screenPatchSize,screenPatchSize,optsTVAL3));
                case 'Wave+TVAL3'
                    scene_test_rec(sceneChannelRange,tst,i) = vec(IWT2_PO(TVAL3(op_Phi_rec_Wave,sensor_test(sensorChannelRange,tst),screenPatchSize,screenPatchSize,optsTVAL3Wave),waveletLevel,QMF));
                case 'Wave+FastaLasso'
                    scene_test_rec(sceneChannelRange,tst,i) = vec(IWT2_PO(reshape(fasta_lasso(op_Phi_rec_WaveF,op_Phi_rec_WaveB,vec(sensor_test(sensorChannelRange,tst)),muFastaLasso, zeros(screenPatchSize^2,1), optsFastaLasso),screenPatchSize,screenPatchSize),waveletLevel,QMF));
                case 'Wave+FastaLSQR'
                    scene_test_rec(sceneChannelRange,tst,i) = vec(IWT2_PO(reshape(fasta_sparseLeastSquares(op_Phi_rec_WaveF,op_Phi_rec_WaveB,sensor_test(sensorChannelRange,tst),muFastaLSQR, zeros(screenPatchSize^2,1), optsFastaLSQR),screenPatchSize,screenPatchSize),waveletLevel,QMF));
                case 'DCT+FastaLasso'
                    scene_test_rec(sceneChannelRange,tst,i) = vec(idct2blk(reshape(fasta_lasso(op_Phi_rec_DCTF,op_Phi_rec_DCTB,sensor_test(sensorChannelRange,tst),muFastaLasso, zeros(screenPatchSize^2,1), optsFastaLasso),screenPatchSize,screenPatchSize)));
                case 'DCT+FastaLSQR'
                    scene_test_rec(sceneChannelRange,tst,i) = vec(idct2blk(reshape(fasta_sparseLeastSquares(op_Phi_rec_DCTF,op_Phi_rec_DCTB,sensor_test(sensorChannelRange,tst),muFastaLSQR, zeros(screenPatchSize^2,1), optsFastaLSQR),screenPatchSize,screenPatchSize)));
                case 'lsqr'
                    scene_test_rec(:,tst,i) = iPhi_rec_lsqr(sensor_test(:,tst));
            end
%             scene_test_rec(sceneChannelRange,tst,i) = sign(sum(scene_test_rec(sceneChannelRange,tst,i)))*scene_test_rec(sceneChannelRange,tst,i);
            if iterationFlag == 1
                rec1 = scene_test_rec(sceneChannelRange,tst,i);
                nrec1 = max((rec1-min(rec1)))-(rec1-min(rec1))+min(rec1);
                sensor_test(:,tst) = sensor_test(:,tst) - (Phi_rec*(nrec1)-mean(Phi_rec*(nrec1)));
            end
        end
    end
end
iterationFlag = iterationFlag+1;
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
            for tst = 1:length(test_indices)
                for i = 1:length(invMode)
                    subplot(length(test_indices),length(invMode),k);
                    thisImage = scene_test_rec(:,tst,i);
                    thisImageNormalized = (thisImage-min(thisImage(:)));
                    thisImageNormalized = reshape(thisImageNormalized/max(thisImageNormalized(:)),screenPatchSize,screenPatchSize,numOfChannels);
                    if BM3DFlag
                        if numOfChannels==1
                            [~, thisImageNormalized] = BM3D(1, thisImageNormalized);
                        else
                            [~, thisImageNormalized] = CBM3D(1, thisImageNormalized);
                        end
                        titleBM3D = '+BM3D';
                    end
                    imshow(thisImageNormalized);title([invMode{i} titleBM3D]);
                    k = k + 1;
                end
            end
        else
            k = 1;
            for tst = 1:length(test_indices)
                for i = 1:length(invMode)
                    subplot(length(test_indices),length(invMode),k);
                    thisImage = scene_test_rec(:,tst,i);
                    thisImageNormalized = (thisImage-min(thisImage(:)));
                    thisImageNormalized = reshape(thisImageNormalized/max(thisImageNormalized(:)),screenPatchSize,screenPatchSize,numOfChannels);
                    if BM3DFlag
                        if numOfChannels==1
                            [~, thisImageNormalized] = BM3D(1, thisImageNormalized);
                            thisImageNormalized = imsharpen(thisImageNormalized);
                        else
                            [~, thisImageNormalized] = CBM3D(1, thisImageNormalized);
                            thisImageNormalized = imsharpen(thisImageNormalized);
                        end
                        titleBM3D = '+BM3D';
                    end
                    imshow(thisImageNormalized);title([invMode{i} titleBM3D]);
                    saveDir = sprintf('Results/mask%s_%dx%d_%s_rec/',maskType,screenPatchSize,screenPatchSize,targetImages);
                    if ~exist(saveDir,'dir'); 
                        mkdir(saveDir);
                    end
                    imwrite(thisImageNormalized,[saveDir '\result_' num2str(test_indices(tst)) '_' invMode{i} titleBM3D '_' num2str(screenPatchSize) '.bmp']);
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