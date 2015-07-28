% Read saved images, crop and resize them, and save in a mat file

clc
close all
clear

% change working directory
mdir = mfilename('fullpath');
mfile = mfilename;
mdir = mdir(1:end-length(mfile));
cd(mdir);

% prepare test data 
screenPatchSize = 128;
hiRes = 512;
maskType = '255M030C';
numberOfFrames = 2;
isDirectory = true;

testDir = sprintf('Data/mask%s_%dx%d_testImages_14in/',maskType,screenPatchSize,screenPatchSize);
lowrestestDir = sprintf('%slowres%d/',testDir, screenPatchSize);
highrestestDir = sprintf('%shighres%d/',testDir, hiRes);

% resize the images...
tmpImage = imread(sprintf('%simage_00000_p_%02d.png',highrestestDir,numberOfFrames));
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

msg1 = sprintf('postCalib separable for testImages: screenPatchSize=%d, hiRes = %d, maskType=%s, downSampling=%02d (%s)',...
    screenPatchSize,hiRes, maskType,downSamplingFactor,downSamplingMode);
fprintf([msg1,'\n']);

%% Read sensor images
vec = @(x) x(:);
cameraSize = [1036 1384];
sensorSize(1:2) = ceil([1024 1024]/downSamplingFactor);
% cropInd = @(z) z(1:end-1024);
resize_temp = @(z) z([1:1024]+(cameraSize(1)-1024)/2,[1:1024]+(cameraSize(2)-1024)/2,:); 
switch downSamplingMode
    case 'average'
        opReadResize = @(nameImg) vec(imresize(double(resize_temp(imread([nameImg])))/numberOfFrames,sensorSize(1:2),'box'));
    case 'subsample'
        downSample = @(z) z(1:downSamplingFactor:end,1:downSamplingFactor:end);
        opReadResize = @(nameImg) vec(downSample(double(resize_temp(imread([nameImg])))/numberOfFrames));
end  

vec = @(z) z(:);
% Test data for validation
if isDirectory 
    l = dir(highrestestDir); l = l(3:end);
    lowresImages = zeros(screenPatchSize^2*numOfChannels,length(l)/2);
    highresImages = zeros(hiRes^2*numOfChannels,length(l)/2);
    
    l = dir(sprintf('patterns/lowresTestC%02d',screenPatchSize)); 
    
    for i = 1:length(l)/2
        % lowres screen images
        tmp = vec(imread(sprintf('patterns/lowresTestC%02d/pattern_%05d.bmp',screenPatchSize,i)));
        if numel(tmp) < size(lowresImages,1);
            tmp = repmat(tmp,[3 1]);
        end
        lowresImages(:,i) = tmp(:);
        % highres screen images
        tmp = vec(imread(sprintf('patterns/highresTestC%02d/pattern_%05d.bmp',hiRes,i)));
        if numel(tmp) < size(highresImages,1);
            tmp = repmat(tmp,[3 1]);
        end
        highresImages(:,i) = tmp;
    end
    lowresImages = reshape([lowresImages; max(lowresImages(:))-lowresImages],size(lowresImages,1),[]);     
    highresImages = reshape([highresImages; max(highresImages(:))-highresImages],size(highresImages,1),[]);
else
    lowresImages = imread(sprintf('patterns/TestC%02d.bmp',screenPatchSize));
    lowresImages = reshape([lowresImages; max(lowresImages(:))-lowresImages],size(lowresImages,1),[]);
    highresImages = imread(sprintf('patterns/TestC%02d.bmp',hiRes));
    highresImages = reshape([highresImages; max(highresImages(:))-highresImages],size(highresImages,1),[]);
end
% sensorTest_lowres = zeros(prod(sensorSize),size(lowresImages,2));
sensorTest_highres = zeros(prod(sensorSize),size(highresImages,2));

fprintf('total test images = %d \n',size(highresImages,2));
msg_old = [];
%
for i = 1:size(highresImages,2)/2
    
    % lowres images
    nameImg = sprintf('%simage_%05d_p_%02d.png', lowrestestDir, i-1,numberOfFrames);
    sensorTest_lowres(:,2*i-1) = opReadResize(nameImg);
    nameImg = sprintf('%simage_%05d_n_%02d.png', lowrestestDir, i-1,numberOfFrames);
    sensorTest_lowres(:,2*i) = opReadResize(nameImg);

    % highres images
    nameImg = sprintf('%simage_%05d_p_%02d.png', highrestestDir, i-1,numberOfFrames);
    sensorTest_highres(:,2*i-1) = opReadResize(nameImg);
    nameImg = sprintf('%simage_%05d_n_%02d.png', highrestestDir, i-1,numberOfFrames);
    sensorTest_highres(:,2*i) = opReadResize(nameImg);
    
    msg = sprintf('running ... %d/%d',2*i,size(highresImages,2));
    fprintf([repmat('\b',1,numel(msg_old)),msg]);
    msg_old = msg;
end
fprintf('\n');

%%

% for ii=1:size(sensorTest_lowres,2); 
%     figure(1); colormap gray;
%     subplot(221); imagesc(reshape(lowresImages(:,ii),screenPatchSize,[]),[0 255]); axis image; 
%     title(ii)
%     subplot(222); imagesc(reshape(highresImages(:,ii),hiRes,[]),[0 255]); axis image; 
%     subplot(223); imagesc(reshape(sensorTest_lowres(:,ii),sensorSize(1),[]),[0 255]); axis image; 
%     subplot(224); imagesc(reshape(sensorTest_highres(:,ii),sensorSize(1),[]),[0 255]); axis image; 
%     shg; pause(1/60); 
% end
%% 
save(sprintf('%stestImages_downSampled%02d_%s.mat',testDir,downSamplingFactor,downSamplingMode),'-v7.3');

% sendmail_lensless({'msalmanasif@gmail.com'},'planar postCalib completed', msg1);
