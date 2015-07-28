% Read saved images, crop and resize them, and save in a mat file

clc
close all
clear

% change working directory
mdir = mfilename('fullpath');
mfile = mfilename;
mdir = mdir(1:end-length(mfile));
cd(mdir);

% prepare calibration data or test data
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
numberOfFrames = 2;

dirPath = sprintf('Data/mask%s_%s%d_14in/',maskType,calibrationSequence,screenPatchSize);

% resize the images...
tmpImage = imread(sprintf('%simage_00000_n_%02d.png',dirPath,numberOfFrames));
cameraSize = size(tmpImage);
downSamplingFactor = 4;
downSamplingMode = 'average';
if length(cameraSize)==3
    sensorSize = [ceil(cameraSize(1:2)/downSamplingFactor) cameraSize(3)];
    numOfChannels = 3;
else
    sensorSize = ceil(cameraSize/downSamplingFactor);
    numOfChannels = 1;
end

msg1 = sprintf('postCalib separable: calibration-%s, screenPatchSize=%d, maskType=%s, downSampling=%02d (%s)',...
    calibrationSequence, screenPatchSize,maskType,downSamplingFactor,downSamplingMode);
fprintf([msg1,'\n']);

%% Read sensor images...
vec = @(x) x(:);
cameraSize = [1036 1384];
sensorSize(1:2) = ceil([1024 1024]/downSamplingFactor);
% cropInd = @(z) z(1:end-1024);
resize_temp = @(z) z([1:1024]+(cameraSize(1)-1024)/2,[1:1024]+(cameraSize(2)-1024)/2,:); 
switch downSamplingMode
    case 'average'
        opReadResize = @(nameImg) vec(imresize(double(resize_temp(imread([dirPath,nameImg])))/numberOfFrames,sensorSize(1:2),'box'));
    case 'subsample'
        downSample = @(z) z(1:downSamplingFactor:end,1:downSamplingFactor:end,:);
        opReadResize = @(nameImg) vec(downSample(double(resize_temp(imread([dirPath,nameImg])))/numberOfFrames));
end


% h = waitbar(0,'Reading images ...');
fprintf('total calibration images = %d \n',N);
msg_old = [];

sensorImageP = zeros(prod(sensorSize),N);
sensorImageN = zeros(prod(sensorSize),N);


% read sensor images; average and resize...
for i = 1:N;
    if calibPM
        % sensorImageP(:,i) = vec(imresize(readImage2(dirPath,sprintf('P_%02d', i),cameraSize,numberOfFrames),cameraPatchSize/2));
        nameImg = sprintf('image_%05d_p_%02d.png', i-1,numberOfFrames);
        sensorImageP(:,i) = opReadResize(nameImg);
        nameImg = sprintf('image_%05d_n_%02d.png', i-1,numberOfFrames);
        sensorImageN(:,i) = opReadResize(nameImg);
        
    else
        % sensorImageP(:,i) = vec(imresize(readImage2(dirPath,sprintf('P_%02d', i),cameraSize,numberOfFrames),cameraPatchSize/2));
        nameImg = sprintf('image_%05d_p_%02d.png', i-1,numberOfFrames);
        sensorImageP(:,i) = opReadResize(nameImg);
    end
    % waitbar(i/totalStep,h);
    msg = sprintf('running ... %d/%d',i,N);
    fprintf([repmat('\b',1,numel(msg_old)),msg]);
    msg_old = msg;
    
    % figure(1); imagesc(reshape([sensorImageP(:,i) sensorImageN(:,i)],sensorSize(1),[])); shg; 
end
fprintf('\n');
% delete(h);

% record zero image for each column
sensorZero = zeros(prod(sensorSize),screenPatchSize);
for i = screenPatchSize-1:screenPatchSize:N;
    nameImg = sprintf('image_%05d_z_%02d.png', i,numberOfFrames);
    sensorZero(:,ceil(i/screenPatchSize)) = opReadResize(nameImg);
end

%% Test data for validation
% 
% testDir = dir(testPath);
% testDir = testDir(3:end);
% screenTest = zeros(N,length(testDir));
% sensorTest = zeros(prod(sensorSize),length(testDir));
% sensorTest_highres = zeros(prod(sensorSize),length(testDir));
% sensorTest_zero = zeros(prod(sensorSize),length(testDir));
% 
% 
% fprintf('total test images = %d \n',length(testDir));
% msg_old = [];
% %
% for i = 1:length(testDir)
%     [pathImage, testImage, extImage] = fileparts(testDir(i).name);
%     
%     nameImg = sprintf(testImage);
%     sensorTest(:,i) = opReadResize(nameImg);
%     screenTest(:,i) = vec(im2double(imresize(imread([testPath testImage extImage]),screenRes,'box')));
%     
%     nameImg = sprintf('highres_%s',testImage);
%     sensorTest_highres(:,i) = opReadResize(nameImg);
%     
%     if exist(sprintf('%stestZero_01_1.bmp',dirPath))
%         nameImg = sprintf('testZero_%02d', i);
%         sensorTest_zero(:,i) = opReadResize(nameImg);
%     end
%     
%     msg = sprintf('running ... %d/%d',i,length(testDir));
%     fprintf([repmat('\b',1,numel(msg_old)),msg]);
%     msg_old = msg;
% end
% fprintf('\n');

%% 
% screenImageP = uint8(screenImageP);
% screenImageN = uint8(screenImageN);
%
% sensorImageP = uint8(sensorImageP);
% sensorImageN = uint8(sensorImageN);
% sensorZero = uint8(sensorZero);
% sensorTest = uint8(sensorTest);
% sensorTest_highres = uint8(sensorTest_highres);
% sensorTest_zero = uint8(sensorTest_zero);

sensorImageP = double(sensorImageP);
sensorImageN = double(sensorImageN);
sensorZero = double(sensorZero);

save(sprintf('%sCalibData_downSampled%02d_%s.mat',dirPath,downSamplingFactor,downSamplingMode),'-v7.3');

% sendmail_lensless({'msalmanasif@gmail.com'},'planar postCalib completed', msg1);
