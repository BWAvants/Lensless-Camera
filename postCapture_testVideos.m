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
maskType = '255M030';
numberOfFrames = 10;
numberOfChannels = 1;

if numberOfChannels>1
    testDir = sprintf('Data/mask%s_testVideosC/',maskType);
else
    testDir = sprintf('Data/mask%s_testVideos/',maskType);
end

testVideos = {'LenslessCamera','LenslessCameraRT','BirdRT'};%
videoFrames = zeros(length(testVideos),1);

downSamplingFactor = 4;
downSamplingMode = 'average';
sensorSize = round([1024 1280]/downSamplingFactor);

%% Read sensor images...
vec = @(x) x(:);
switch downSamplingMode
    case 'average'
        opReadResize = @(nameImg) vec(imresize(double(imread([nameImg]))/numberOfFrames,sensorSize,'box'));
    case 'subsample'
        downSample = @(z) z(1:downSamplingFactor:end,1:downSamplingFactor:end);
        opReadResize = @(nameImg) vec(downSample(double(imread([nameImg]))/numberOfFrames));
end 

msg1 = sprintf('postCapture separable for realImages: maskType=%s, downSampling=%02d (%s)',...
     maskType,downSamplingFactor,downSamplingMode);
fprintf([msg1,'\n']);
fprintf('total test images = %d \n',length(testVideos));
msg_old = [];


sensorTest = cell(length(testVideos),1);
for k = 1:length(testVideos)
    testVideoDir = [testDir testVideos{k} '/'];
    testVideoDirLS = dir(testVideoDir); testVideoDirLS = testVideoDirLS(3:end);
    if strcmp(testVideos{k}(end-1:end),'RT')
        videoFrames(k) = length(testVideoDirLS);
    else
        videoFrames(k) = length(testVideoDirLS)/2;
    end
    
    tmpImage = imread([testVideoDir testVideoDirLS(1).name]);
    cameraSize = size(tmpImage);
    sensorSize = ceil(cameraSize/downSamplingFactor);
    sensorTest_video = zeros(prod(sensorSize),length(testVideoDirLS));
    
    if strcmp(testVideos{k}(end-1:end),'RT')
        for i = 1:videoFrames(k)
            sensorTest_video(:,i) = opReadResize([testVideoDir testVideoDirLS(i).name]);
            msg = sprintf('running ... %d/%d',i,videoFrames(k));
            fprintf([repmat('\b',1,numel(msg_old)),msg]);
            msg_old = msg;
        end
    else
        for i = 1:videoFrames(k)
            % lowres images
            nameImg = sprintf('%simage_%05d_p_%d.png', testVideoDir, i-1,numberOfFrames);
            sensorTest_video(:,2*i-1) = opReadResize(nameImg);
            nameImg = sprintf('%simage_%05d_n_%d.png', testVideoDir, i-1,numberOfFrames);
            sensorTest_video(:,2*i) = opReadResize(nameImg);

            msg = sprintf('running ... %d/%d',i,videoFrames(k));
            fprintf([repmat('\b',1,numel(msg_old)),msg]);
            msg_old = msg;
        end
    end
    sensorTest(k) = {sensorTest_video};
end
fprintf('\n');

%%
clear sensorTest_video
save(sprintf('%stestVideos_downSampled%02d_%s.mat',testDir,downSamplingFactor,downSamplingMode),'-v7.3');

% sendmail_lensless({'msalmanasif@gmail.com'},'planar postCalib completed', msg1);
