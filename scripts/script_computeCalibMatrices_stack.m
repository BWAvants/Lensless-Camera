%% Read stored images
tmp_mat = matfile(sprintf('%sCalibData_downSampled%02d_%s.mat',calibDir,downSamplingFactor,downSamplingMode));
sensorSize = tmp_mat.sensorSize;
numOfChannels = tmp_mat.numOfChannels;
% load(sprintf('%sCalibData_downSampled%02d_%s.mat',calibDir,downSamplingFactor,downSamplingMode));
calibDir = sprintf('Data/mask%s_%s%d/',maskType,calibrationSequence,screenPatchSize);

skipH1 = 1; 

fprintf('reading images ... ');
sensorImageP = (tmp_mat.sensorImageP-repmat((mean(tmp_mat.sensorZero,2)),1,size(tmp_mat.sensorImageP,2)));
sensorImageN = (tmp_mat.sensorImageN-repmat((mean(tmp_mat.sensorZero,2)),1,size(tmp_mat.sensorImageN,2)));
clear sensorZero
sensorImage_leftC = sensorImageP(:,1:end/2);
sensorImage_rightC = sensorImageP(:,end/2+1:end);
clear sensorImageP
sensorImage_leftC = [sensorImage_leftC sensorImageN(:,1:end/2)];
sensorImage_rightC = [sensorImage_rightC sensorImageN(:,end/2+1:end)];
clear sensorImageN
%
% remove hot/dark pixels
% remove_pixels = [];
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
calibMatrix = imread(sprintf('patterns/%s%d.bmp',calibrationSequence,screenPatchSize));
screenImage = double(calibMatrix(1:screenPatchSize,:))/double(max(calibMatrix(:)));
%
clear calibMatrix
screenImageP = screenImage;
% screenImageP(screenImageP<0) = 0;
screenImageN = -screenImage+1;
% screenImageN(screenImageN<0) = 0;
clear screenImage
screenImage_left = [screenImageP(:,1:end/2) screenImageN(:,1:end/2)]; % e.g., hadamard pattern repeated in columns (Im = H(:,i)*ones(1,n))
screenImage_right = [screenImageP(:,end/2+1:end) screenImageN(:,end/2+1:end)]; % e.g., hadamard pattern repeated in rows
sensorImageP = []; sensorImageN = [];
screenImageP = []; screenImageN = [];
sensorZero = [];

fprintf('done! \n');


% different ways to calculate separable components
% rank = 4; op_makeSeparable = @(z) z; % rank-2 without mean subtraction
rank = 1; op_makeSeparable = @(z) subtract_mean(z); % rank-1 after mean substraction..

sensorImage_left = reshape(mean(reshape(sensorImage_leftC, sensorSize(1), sensorSize(2),sensorSize(3),[]),3),[],size(sensorImage_leftC,2));
sensorImage_right = reshape(mean(reshape(sensorImage_rightC, sensorSize(1), sensorSize(2),sensorSize(3),[]),3),[],size(sensorImage_rightC,2));
sensorImage_leftC = []; sensorImage_rightC = [];
% the "right" way to compute calibration matrix using svd...
% fprintf('compute joint svd of all the left images for left matrix and all the right images for the right matrix\n');
ItI_left = 0; ItI_right = 0;
Itl_stack = zeros(2*sensorSize(1)*screenPatchSize,sensorSize(2));
Itr_stack = zeros(2*sensorSize(2)*screenPatchSize,sensorSize(1));
for ii =1:size(sensorImage_left,2);
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

% % test rank..
% for ii=1:screenPatchSize;
%     r=1;
%     figure(1); subplot(211);
%     imagesc([Itl_stack((ii-1)*sensorSize(1)+1:ii*sensorSize(1),:)  Itl_Vl((ii-1)*sensorSize(1)+1:ii*sensorSize(1),1:r)*Vl(:,1:r)']);
%     title(ii);
%     subplot(212);
%     imagesc([Itr_stack((ii-1)*sensorSize(2)+1:ii*sensorSize(2),:) Itr_Vr((ii-1)*sensorSize(2)+1:ii*sensorSize(2),1:r)*Vr(:,1:r)']);
%     shg; drawnow; pause;
% end

Itr_Vr = []; Itl_Vl = [];

Phi_rec_left = zeros(sensorSize(1),screenPatchSize,rank);
Phi_rec_right = zeros(sensorSize(2),screenPatchSize,rank);
for r=1:rank
    % Phi_rec_left(:,:,r) = Phi_left(:,:,r)*(screenImage_left(1:screenPatchSize,:)')/norm((screenImage_left(1:screenPatchSize,1)))^2;
    % Phi_rec_right(:,:,r) = Phi_right(:,:,r)*(screenImage_left(1:screenPatchSize,:)')/norm((screenImage_left(1:screenPatchSize,1)))^2;
    %Phi_rec_left(:,:,r) = Phi_left(:,:,r)/screenImage_left(1:screenPatchSize,:);    
    %Phi_rec_right(:,:,r) = Phi_right(:,:,r)/screenImage_left(1:screenPatchSize,:);    
    Phi_rec_left(:,:,r) = Phi_left(:,:,r)*(pinv(screenImage_left')');    
    Phi_rec_right(:,:,r) = Phi_right(:,:,r)*(pinv(screenImage_left')');    
    
end
for r = 1:rank
    sum_Phi_rec_left = sum(Phi_rec_left(:,:,r),2);
    sum_Phi_rec_right = sum(Phi_rec_right(:,:,r),2);
    flip_sign =  sign(sum_Phi_rec_left'*Vr(:,r));
    % figure(r);
    % subplot(211); plot([sum_Phi_rec_left/norm(sum_Phi_rec_left) Vr(:,r)/norm(Vr(:,r))])
    % subplot(212); plot([sum_Phi_rec_right/norm(sum_Phi_rec_right) Vl(:,r)/norm(Vl(:,r))])
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

%% test calibration matrix...
tst = randint(1,1,[1 screenPatchSize])
xTemp = repmat(screenImage_left(1:screenPatchSize,tst),1,screenPatchSize);
It = reshape(sensorImage_left(:,tst),sensorSize(1:2));
It = op_makeSeparable(It);
yTemp = Phi_rec_left(:,:,1)*xTemp*Phi_rec_right(:,:,1)';
figure(1); clf; imagesc([yTemp It It-yTemp]);

saveName = sprintf('%sCalibCompute_downSampled%02d_%s_%s.mat',calibDir,downSamplingFactor,downSamplingMode,calibMode);

save(saveName,'-v7.3');
fprintf('done! \n');

