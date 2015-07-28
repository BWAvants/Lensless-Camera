% Post processing on averaged images

% clc
close all
clear

addpath 'misc/'
addpath 'utils/'
addpath 'utils/lsqrSOL/'
addpath 'utils/export_fig/';
% addpath 'misc\export_fig\';
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
maskType = '255M030C';
switch maskType
    case {'16Mseq','08Mseq','031M100','127M060','255M030','255M030C','127M060C','255M010C'}
        op_makeSeparable = @(z) subtract_mean(z);
    case {'16Sseq','2xPA060'}
        op_makeSeparable = @(z) z;
end
N = screenPatchSize^2;
%
downSamplingFactor = 4;
downSamplingMode = 'average'; 

msg1 = sprintf('computeCalib separable: calibration-%s, screenPatchSize=%d, maskType=%s, downSampling=%02d (%s)',...
    calibrationSequence, screenPatchSize,maskType,downSamplingFactor,downSamplingMode);
fprintf([msg1,'\n']);

% calibration images
calibDir = sprintf('Data/mask%s_%s%d/',maskType,calibrationSequence,screenPatchSize);

%% Read stored images
calibCompute_mat = sprintf('%sCalibCompute_downSampled%02d_%s.mat',calibDir,downSamplingFactor,downSamplingMode);

if ~exist(calibCompute_mat)
    tmp_mat = matfile(sprintf('%sCalibData_downSampled%02d_%s.mat',calibDir,downSamplingFactor,downSamplingMode));
    sensorSize = tmp_mat.sensorSize;
    numOfChannels = tmp_mat.numOfChannels;
    % load(sprintf('%sCalibData_downSampled%02d_%s.mat',calibDir,downSamplingFactor,downSamplingMode));
    calibDir = sprintf('Data/mask%s_%s%d/',maskType,calibrationSequence,screenPatchSize);
    
    
    %     if numOfChannels > 1;
    %         sensorSizeC = sensorSize;
    %         sensorSize = [sensorSizeC(1) sensorSizeC(2)*sensorSizeC(3)];
    %     end
    
    fprintf('reading images ... ');
    
    % calibration images
%     sensorImageP = double(sensorImageP);
%     sensorImageN = double(sensorImageN);
%     sensorZero = double(sensorZero);
%     sensorImageP = (sensorImageP-repmat((mean(sensorZero,2)),1,size(sensorImageP,2)));
%     sensorImageN = (sensorImageN-repmat((mean(sensorZero,2)),1,size(sensorImageN,2)));
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
     
    
    fprintf('done! \n');
     
    %% compute inverse matrix ..
    invType = 'pm';
    
    fprintf('computing calibration matrix using %s ... \n ',invType);
    switch calibrationSequence
        case 'hadamard'
            % Phi_pm=(sensorImageP-sensorImageN)/(screenImageP-screenImageN);
            Phi_rec=(sensorImageP-sensorImageN)*(((screenImageP-screenImageN)')/sum((screenImageP(:,1)-screenImageN(:,1)).^2));
        case 'hadamardSeparable'
%             sensorImageP = double(tmp_mat.sensorImageP);
%             sensorImageN = double(tmp_mat.sensorImageN);
%             sensorZero = double(tmp_mat.sensorZero);
            sensorImageP = (tmp_mat.sensorImageP-repmat((mean(tmp_mat.sensorZero,2)),1,size(tmp_mat.sensorImageP,2)));
            sensorImageN = (tmp_mat.sensorImageN-repmat((mean(tmp_mat.sensorZero,2)),1,size(tmp_mat.sensorImageN,2)));
            clear sensorZero
            sensorImage_leftC = sensorImageP(:,1:end/2);
            sensorImage_rightC = sensorImageP(:,end/2+1:end);
            clear sensorImageP
            sensorImage_leftC = sensorImage_leftC-sensorImageN(:,1:end/2);
            sensorImage_rightC = sensorImage_rightC - sensorImageN(:,end/2+1:end);
            clear sensorImageN
            
            % read screen image from the pattern file..
            calibMatrix = imread(sprintf('patterns/%s%d.bmp',calibrationSequence,screenPatchSize));
            screenImage = double(calibMatrix(1:screenPatchSize,:))/double(max(calibMatrix(:)));
            %
            clear calibMatrix
            screenImageP = screenImage;
            % screenImageP(screenImageP<0) = 0;
            screenImageN = -screenImage+1;
            % screenImageN(screenImageN<0) = 0;
            clear screenImage
            screenImage_left = screenImageP(:,1:end/2)-screenImageN(:,1:end/2); % e.g., hadamard pattern repeated in columns (Im = H(:,i)*ones(1,n))
            screenImage_right = screenImageP(:,end/2+1:end)-screenImageN(:,end/2+1:end); % e.g., hadamard pattern repeated in rows
            sensorImageP = []; sensorImageN = [];
            screenImageP = []; screenImageN = [];
            sensorZero = [];
            
            % different ways to calculate separable components
            % rank = 4; op_makeSeparable = @(z) z; % rank-2 without mean subtraction
            rank = 1; op_makeSeparable = @(z) subtract_mean(z); % rank-1 after mean substraction..
            Phi_rec_leftC = zeros(sensorSize(1),screenPatchSize,rank,numOfChannels);
            Phi_rec_rightC = zeros(sensorSize(2),screenPatchSize,rank,numOfChannels);
             
            for C =1:numOfChannels
                sensorImage_left = sensorImage_leftC((C-1)*prod(sensorSize(1:2))+1:C*prod(sensorSize(1:2)),:);
                sensorImage_right = sensorImage_rightC((C-1)*prod(sensorSize(1:2))+1:C*prod(sensorSize(1:2)),:);
                % the "right" way to compute calibration matrix using svd...
                % fprintf('compute joint svd of all the left images for left matrix and all the right images for the right matrix\n');
                ItI_left = 0; ItI_right = 0;
                Itl_stack = zeros(sensorSize(1)*screenPatchSize,sensorSize(2));
                Itr_stack = zeros(sensorSize(2)*screenPatchSize,sensorSize(1));
                for ii =1:screenPatchSize;
                    It = reshape(sensorImage_left(:,ii),sensorSize(1:2));
                    Itl = op_makeSeparable(It);
                    ItI_left = ItI_left + Itl'*Itl;
                    Itl_stack((ii-1)*sensorSize(1)+1:ii*sensorSize(1),:) = Itl;
                    
                    It = reshape(sensorImage_right(:,ii),sensorSize(1:2))';
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
                    Phi_rec_left(:,:,r) = Phi_left(:,:,r)*(screenImage_left(1:screenPatchSize,:)')/norm((screenImage_left(1:screenPatchSize,1)))^2;
                    Phi_rec_right(:,:,r) = Phi_right(:,:,r)*(screenImage_left(1:screenPatchSize,:)')/norm((screenImage_left(1:screenPatchSize,1)))^2;
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
                %% normalize the matrices...
                tst = 1;%randint(1,1,[1 screenPatchSize])
                xTemp = repmat(screenImage_left(1:screenPatchSize,tst),1,screenPatchSize);
                It = reshape(sensorImage_left(:,tst),sensorSize(1:2));
                It = op_makeSeparable(It);
                yTemp = 0;
                for r = 1:rank
                    yTemp = yTemp+ Phi_rec_left(:,:,r)*xTemp*Phi_rec_right(:,:,r)';
                end
                normalizationFactor = sqrt(norm(It(:))/norm(yTemp(:)));
                Phi_rec_left = Phi_rec_left*normalizationFactor;
                Phi_rec_right = Phi_rec_right*normalizationFactor;
                  
                Phi_rec_leftC(:,:,:,C) = Phi_rec_left;
                Phi_rec_rightC(:,:,:,C) = Phi_rec_right;
                
                
                
                %% test calibration matrix...
                tst = randint(1,1,[1 screenPatchSize])
                xTemp = repmat(screenImage_left(1:screenPatchSize,tst),1,screenPatchSize);
                It = reshape(sensorImage_left(:,tst),sensorSize(1:2));
                It = op_makeSeparable(It);
                yTemp = Phi_rec_left(:,:,1,1)*xTemp*Phi_rec_right(:,:,1,1)';
                figure(1); clf; imagesc([yTemp It It-yTemp]);
            end
            Phi_rec_left = []; Phi_rec_right = [];
%             sensorImage_left = []; sensorImage_right = [];
    end
    
    saveName = sprintf('%sCalibCompute_downSampled%02d_%s.mat',calibDir,downSamplingFactor,downSamplingMode);
    
    save(saveName,'-v7.3');
    fprintf('done! \n');
end
fprintf('done!\n');
% figure(555); imagesc([Phi_stack Phi_mux Phi_pm])

%% setup inverse problem
vec = @(z) z(:);
saveName = sprintf('%sCalibCompute_downSampled%02d_%s.mat',calibDir,downSamplingFactor,downSamplingMode);

tmp_mat = matfile(saveName);
sensorSize = tmp_mat.sensorSize;
Phi_rec_leftC = tmp_mat.Phi_rec_leftC;
Phi_rec_rightC = tmp_mat.Phi_rec_rightC;
screenImage_left = tmp_mat.screenImage_left;
screenImage_right = tmp_mat.screenImage_right;
numOfChannels = tmp_mat.numOfChannels;
% flip_sign_r = tmp_mat.flip_sign_r;
% flip_sign_l = tmp_mat.flip_sign_l;

invMode = 'svd';
fprintf('computing inverse matrix using %s ...  ',invMode);
switch invMode
    case 'svd'
        % svd for left- and right-side matrices
        iPhi_rec_left_matrix = [];
        iPhi_rec_right_matrix = [];
        Ul = zeros(size(Phi_rec_leftC,1),screenPatchSize,numOfChannels);
        Ur = zeros(size(Phi_rec_rightC,1),screenPatchSize,numOfChannels);
        Vl = zeros(screenPatchSize,size(Phi_rec_leftC,2),numOfChannels);
        Vr = zeros(screenPatchSize,size(Phi_rec_rightC,2),numOfChannels);
        iSl = zeros(screenPatchSize,screenPatchSize,numOfChannels);
        iSr = zeros(screenPatchSize,screenPatchSize,numOfChannels);
        Sl = zeros(screenPatchSize,numOfChannels);
        Sr = zeros(screenPatchSize,numOfChannels);
        for C = 1:numOfChannels
%             sum_Phi_rec_leftC = sum(Phi_rec_leftC(:,:,C),2);
%             sum_Phi_rec_rightC = sum(Phi_rec_rightC(:,:,C),2);
            [Ul(:,:,C) Sltemp Vl(:,:,C)] = svd(Phi_rec_leftC(:,:,C),'econ'); Sl(:,C) = diag(Sltemp);
            [Ur(:,:,C) Srtemp Vr(:,:,C)] = svd(Phi_rec_rightC(:,:,C),'econ'); Sr(:,C) = diag(Srtemp);
            
%             flip_sign_l(C) =  sign(sum_Phi_rec_leftC'*Vr(:,r));
%             flip_sign_r(C) =  sign(sum_Phi_rec_rightC'*Vl(:,r));
            
            T = round(screenPatchSize/4); lambda =1;
            %         iSl = speye(screenPatchSize);
            iSl(:,:,C) = diag(1./(Sl(:,C)+lambda));
            %         iSr = speye(screenPatchSize);
            iSr(:,:,C) = diag(1./(Sr(:,C)+lambda));
%             flipSignl = sign(sum(sum(Vl(:,1:T,C)*iSl(1:T,1:T,C)*Ul(:,1:T,C)')));
%             flipSignr = sign(sum(sum(Vr(:,1:T,C)*iSr(1:T,1:T,C)*Ur(:,1:T,C)')));
%             flipSign = min([flipSignl flipSignr]);
            iPhi_rec_left_matrix = [iPhi_rec_left_matrix; Vl(:,1:T,C)*iSl(1:T,1:T,C)*Ul(:,1:T,C)'];
            iPhi_rec_right_matrix = blkdiag(iPhi_rec_right_matrix, Vr(:,1:T,C)*iSr(1:T,1:T,C)*Ur(:,1:T,C)');
        end
        trim = @(z) z(:,(1:numOfChannels).*(1:numOfChannels));
        iPhi_rec_left = @(z) iPhi_rec_left_matrix *reshape(z,sensorSize(1),[]);
        iPhi_rec_right = @(z) reshape(z,sensorSize(1),[])*iPhi_rec_right_matrix';
        iPhi_rec = @(z) vec(trim(im2col(iPhi_rec_left(iPhi_rec_right(z)),[screenPatchSize screenPatchSize],'distinct')));
        
%         op_Phi_rec = @(z) vec(Phi_rec_left*(reshape(z,screenPatchSize,screenPatchSize)*Phi_rec_right'));
        
        
        %% Test resolution of reconstruction
        %         test_indices = [39 71 129 149 159 163 165]+0;
        %         T = round(screenPatchSize/1.0); lambda = 1;
        %         iSl = speye(screenPatchSize); iSl = spdiags(1./(Sl+lambda),0,iSl);
        %         iSr = speye(screenPatchSize); iSr = spdiags(1./(Sr+lambda),0,iSr);
        %         iPhi_rec_left = @(z) Vl(:,1:T)*(iSl(1:T,1:T)*(Ul(:,1:T)'*reshape(z,sensorSize(1),[])));
        %         iPhi_rec_right = @(z) reshape(z,[],sensorSize(2))*Vr(:,1:T)*iSr(1:T,1:T)*Ur(:,1:T)';
        %         iPhi_rec = @(z) vec(iPhi_rec_right(iPhi_rec_left(z)));
        %
        %         LSRec = zeros(screenPatchSize^2,length(test_indices));
        %         ProjRec = zeros(screenPatchSize^2,length(test_indices));
        %         tst=1;
        %         for ii=test_indices
        %         % for ii=1:length(test_indices);
        %             LSRec(:,tst) = iPhi_rec(sensor_test(:,ii));
        %             ProjRec(:,tst) = vec((Vl(:,1:T)*Vl(:,1:T)')*reshape(scene_test(:,ii),screenPatchSize,[])*(Ur(:,1:T)*Ur(:,1:T)'));
        %             tst = tst+1;
        %         end
        %         figure(1); imagesc([reshape(scene_test(:,test_indices),screenPatchSize,[]); reshape(LSRec,screenPatchSize,[]); reshape(ProjRec,screenPatchSize,[]);],[0 1]); colormap gray;
        %         % figure(1); imagesc([reshape(scene_test,screenPatchSize,[]); reshape(LSRec,screenPatchSize,[]); reshape(ProjRec,screenPatchSize,[]);],[0 1]); colormap gray;
        %         figure(2); imagesc([Vl Ur]);
        
    case 'lsqr'
        % lambda = sqrt(1);
        % f_h = @(z) op_separableSystem(z, Phi_rec_left, Phi_rec_right, screenPatchSize,sensorSize,lambda);
        % ft_h = @(z) adj_separableSystem(z, Phi_rec_left, Phi_rec_right, screenPatchSize,sensorSize,lambda);
        %
        % fhandle = @(z,mode) fhandle_mode(z,mode, f_h, ft_h);
        % tol = 1e-6;
        % maxit = 1000;
        % iPhi_rec = @(z) lsqr(fhandle,[z;zeros(screenPatchSize^2*(lambda>0),1)],tol,maxit);
        %
        op_Phi_rec = @(z) op_separableSystem(z, Phi_rec_left, Phi_rec_right, screenPatchSize,sensorSize);
        
        damp = 1;
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
fprintf('done! \n');

% test images
testDir = sprintf('Data/mask%s_%dx%d_testImages/',maskType,screenPatchSize,screenPatchSize);
%
% lowres test images
tmp_mat = matfile(sprintf('%stestImages_downSampled%02d_%s.mat',testDir,downSamplingFactor,downSamplingMode));
lowrestestImages = double(tmp_mat.lowresImages);
sensorTest_lowres = tmp_mat.sensorTest_lowres;
% highres test images
tmp_mat = matfile(sprintf('%stestImages_downSampled%02d_%s.mat',testDir,downSamplingFactor,downSamplingMode));
highrestestImages = double(tmp_mat.highresImages);
sensorTest_highres = tmp_mat.sensorTest_highres;


%% lowres
fprintf('estimating lowres test images ... ');
testDir = sprintf('Data/mask%s_%dx%d_testImages/',maskType,screenPatchSize,screenPatchSize);

% test images
% sensorTest_zero = sensorTest_lowres(:,end-1);
sensorTest = sensorTest_lowres(:,1:end);
screenTest = lowrestestImages(:,1:end);
if saveFigs
    test_indices = find(sum(abs(sensorTest),1)>0);
else
    test_indices = [21 125:2:142];
end

% sensorDC = sensorZero;
% sensorDC = sensorTest_zero;
sensorDC = zeros(size(sensorTest,1),1);
sensor_test = sensorTest(:,test_indices)-repmat(mean(sensorDC,2),1,length(test_indices));
%
% remove faulty pixels
remove_pixels = [];
sensor_test(remove_pixels,:) = [];
%
scene_test = screenTest(:,test_indices)/max(screenTest(:));

% simulated measurements and reconstruction for test images
sensor_test_simulated = 0*sensor_test;
scene_test_rec = zeros(screenPatchSize^2*numOfChannels,length(test_indices));
for tst = 1:length(test_indices)
    
    It = reshape(sensor_test(:,tst),sensorSize);
    if strcmpi(calibrationSequence,'hadamardSeparable')
        for C = 1:numOfChannels
            It(:,:,C) = op_makeSeparable(It(:,:,C));
        end
        sensor_test(:,tst) = It(:);
    end
%     sensor_test_simulated(:,tst) = op_Phi_rec(scene_test(:,tst));
    scene_test_rec(:,tst) = iPhi_rec(sensor_test(:,tst));
    for i = 1:numOfChannels
        scene_test_rec(screenPatchSize^2*(i-1)+1:screenPatchSize^2*i,tst) = sign(sum(scene_test_rec(screenPatchSize^2*(i-1)+1:screenPatchSize^2*i,tst)))*...
            scene_test_rec(screenPatchSize^2*(i-1)+1:screenPatchSize^2*i,tst);
    end
end
% scene_test_rec = bsxfun(@rdivide, scene_test_rec, max(abs(scene_test_rec),[],1));

% f771=figure(771);
% % subplot(411); images([sensor_test sensor_test_simulated sensor_test-sensor_test_simulated])
% subplot(411); plot([vec(sensor_test) vec(sensor_test_simulated) vec(sensor_test-sensor_test_simulated)]); axis tight
% title(sprintf([msg1, '\n lowres: comparison between observed and simulated measurements (original, simulated, error)']));
% subplot(412); plot([vec(scene_test) vec(scene_test_rec) vec(scene_test-scene_test_rec)]); axis tight;
% title(sprintf('normalized reconstruction error: %3.4g',norm(vec(scene_test-scene_test_rec))/norm(scene_test(:))));
% subplot(4,1,3:4); imagesc([reshape(scene_test,screenPatchSize,[]); reshape(scene_test_rec,screenPatchSize,[])])
% title('comparison between projected and reconstructed signal (original, simulated, error)'); colormap gray;
% %
% maximize(f771); set(gcf,'Color','w');

rearrangeToShow = @(z) reshape(permute(reshape(z,screenPatchSize,screenPatchSize,numOfChannels,length(test_indices)),[1 2 4 3]),screenPatchSize,[],numOfChannels);

f881=figure(881); imshow([rearrangeToShow(scene_test); rearrangeToShow(scene_test_rec)]);
title(sprintf([msg1, '\n lowres test images']));
if numOfChannels==1 
    colormap gray;
end
axis image; axis off;
%
maximize(f881); set(gcf,'Color','w');
if saveFigs
    %     figName = sprintf('recovery_calib%s_%dx%d_%s_misc_DS%02d%s', ...
    %         calibrationSequence, screenPatchSize, screenPatchSize,maskType,downSamplingFactor,downSamplingMode);
    %     export_fig(f771,[testDir,figName],'-png','-pdf');
    % saveas(f771,[testDir,figName],'fig');
    figName = sprintf('recovery_calib%s_%dx%d_%s_DS%02d%s', ...
        calibrationSequence, screenPatchSize, screenPatchSize,maskType,downSamplingFactor,downSamplingMode);
    export_fig(f881,[testDir,figName],'-png','-pdf');
    saveas(f881,[testDir,figName],'fig');
end
fprintf('done! \n');

%% highres
fprintf('estimating highres test images ... ');

% test images
% sensorTest_zero = sensorTest_highres(:,end-1);
sensorTest = sensorTest_highres(:,1:end);
screenTest = lowrestestImages(:,1:end);

% sensorDC = sensorZero;
% sensorDC = sensorTest_zero;
sensorDC = zeros(size(sensorTest,1),1);
sensor_test = sensorTest(:,test_indices)-repmat(mean(sensorDC,2),1,length(test_indices));
%
% remove faulty pixels
sensor_test(remove_pixels,:) = [];
%
scene_test = screenTest(:,test_indices)/max(screenTest(:));

% simulated measurements and reconstruction for test images
sensor_test_simulated = 0*sensor_test;
scene_test_rec = zeros(screenPatchSize^2*numOfChannels,length(test_indices));
for tst = 1:length(test_indices)
    It = reshape(sensor_test(:,tst),sensorSize);
    if strcmpi(calibrationSequence,'hadamardSeparable')
        for C = 1:numOfChannels
            It(:,:,C) = op_makeSeparable(It(:,:,C));
        end
        sensor_test(:,tst) = It(:);
    end
%     sensor_test_simulated(:,tst) = op_Phi_rec(scene_test(:,tst));
    scene_test_rec(:,tst) = iPhi_rec(sensor_test(:,tst));
    for i = 1:numOfChannels
        scene_test_rec(screenPatchSize^2*(i-1)+1:screenPatchSize^2*i,tst) = sign(sum(scene_test_rec(screenPatchSize^2*(i-1)+1:screenPatchSize^2*i,tst)))*...
            scene_test_rec(screenPatchSize^2*(i-1)+1:screenPatchSize^2*i,tst);
    end
end

% % f772=figure(772);
% % % subplot(411); imagesc([sensor_test sensor_test_simulated sensor_test-sensor_test_simulated])
% % subplot(411); plot([vec(sensor_test) vec(sensor_test_simulated) vec(sensor_test-sensor_test_simulated)]); axis tight
% % title(sprintf([msg1, '\n highres: comparison between observed and simulated measurements (original, simulated, error)']));
% % subplot(412); plot([vec(scene_test) vec(scene_test_rec) vec(scene_test-scene_test_rec)]); axis tight;
% % title(sprintf('normalized reconstruction error: %3.4g',norm(vec(scene_test-scene_test_rec))/norm(scene_test(:))));
% % subplot(4,1,3:4); imagesc([reshape(scene_test,screenPatchSize,[]); reshape(scene_test_rec,screenPatchSize,[])])
% % title('comparison between projected and reconstructed signal (original, simulated, error)'); colormap gray;
% % %
% % maximize(f772); set(gcf,'Color','w');
% 
% 
f882=figure(882); imshow([rearrangeToShow(scene_test); rearrangeToShow(scene_test_rec)],[0 1]);
title(sprintf([msg1, '\n highres test images']));
if numOfChannels==1 
    colormap gray;
end
axis image; axis off;
%
maximize(f882); set(gcf,'Color','w');

if saveFigs
    %     figName = sprintf('recovery_calib%s_%dx%d_%s_misc_DS%02d%s_highres', ...
    %         calibrationSequence, screenPatchSize, screenPatchSize,maskType,downSamplingFactor,downSamplingMode);
    %     export_fig(f772, [testDir,figName],'-png','-pdf');
    % saveas(f772,[testDir,figName],'fig');
    figName = sprintf('recovery_calib%s_%dx%d_%s_DS%02d%s_highres', ...
        calibrationSequence, screenPatchSize, screenPatchSize,maskType,downSamplingFactor,downSamplingMode);
    export_fig(f882, [testDir,figName],'-png','-pdf');
    saveas(f882,[testDir,figName],'fig');
end
fprintf('done! \n');
