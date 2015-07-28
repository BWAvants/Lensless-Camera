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
maskType = '255M030C';
numberOfFrames = 20;

testDir = sprintf('Data/mask%s_realImages/',maskType);

testImages = {'Angrybird','Animals'};
% testImages = {'Blocks-R','Blocks-G','Blocks-B','Toys-R','Toys-G','Toys-B','Angry Bird-R','Angry Bird-G','Angry Bird-B',...
%     'DSP-R','DSP-G','DSP-B','Magazine-R','Magazine-G','Magazine-B','Books-R','Books-G','Books-B','Fruits-R','Fruits-G','Fruits-B'};%
% testImages = {'angry bird 2ft CBG','angry bird','angry bird WBG','DSP mug','DSP mug WBG','DSP mug2 WBG', 'fruits','fruits WBG','fruits2 WBG','checkerboard2','toys','toys WBG','books','books2','books3','books4','R','1-R','snacks','snacks2'};

%% For mas255M030_hadamardSeparable256
% testImages = {'angry-bird','Blocks2','Blocks3','DSP-Pooh2','Books','Books2','poster3','posterII','posterII2','DSP-Apple','Fruits','MagazineII2','Lovet2'};
%% For mask127M060_hadamardSeparable128
% testImages = {'angry-bird6','books3','DSP-apple','DSP-apple2','DSP-apple-CHKB'...
%     ,'DSP-pooh3','DSP-pooh2 WBG','fruits7','fruits6 WBG','magazine4','poster4','posterII3',...
%     'test-pattern-4in','test-pattern-5in','test-pattern-6in'};
%%


for downSamplingFactor = [4];
    downSamplingMode = 'average';
    numOfChannels = 3; 
    cameraSize = [1036 1384 numOfChannels];
    sensorSize = cameraSize; 
    sensorSize(1:2) = ceil([1024 1024]/downSamplingFactor);
    
    %% Read sensor images...
    vec = @(x) x(:);
    resize_temp = @(z) z([1:1024]+(cameraSize(1)-1024)/2,[1:1024]+(cameraSize(2)-1024)/2,:);    
    switch downSamplingMode
        case 'average'
            opReadResize = @(nameImg) vec(imresize(double(resize_temp(imread([nameImg])))/numberOfFrames,sensorSize(1:2),'box'));
        case 'subsample'
            downSample = @(z) z(1:downSamplingFactor:end,1:downSamplingFactor:end);
            opReadResize = @(nameImg) vec(downSample(double(resize_temp(imread([nameImg])))/numberOfFrames));
    end
    
    msg1 = sprintf('postCapture separable for realImages: maskType=%s, downSampling=%02d (%s)',...
        maskType,downSamplingFactor,downSamplingMode);
    fprintf([msg1,'\n']);
    fprintf('total test images = %d \n',length(testImages));
    msg_old = [];
    sensorTest = zeros(prod(sensorSize),length(testImages));
    for k = 1:length(testImages)
        imageDir = [testDir testImages{k} '/'];
        imageDirLS = dir(imageDir); imageDirLS = imageDirLS(3:end);
        
        for i = 1:length(imageDirLS)
            if strcmp(imageDirLS(i).name(1:4),'zero')
                sensorTest(:,k) = sensorTest(:,k) - opReadResize([imageDir imageDirLS(i).name]);
            else
                sensorTest(:,k) = sensorTest(:,k) + opReadResize([imageDir imageDirLS(i).name]);
            end
        end
        
        msg = sprintf('running ... %d/%d',k,length(testImages));
        fprintf([repmat('\b',1,numel(msg_old)),msg]);
        msg_old = msg;
    end
    fprintf('\n');
    
    %%
    
    save(sprintf('%srealImages_downSampled%02d_%s.mat',testDir,downSamplingFactor,downSamplingMode),'-v7.3');
    
    % sendmail_lensless({'msalmanasif@gmail.com'},'planar postCalib completed', msg1);
end