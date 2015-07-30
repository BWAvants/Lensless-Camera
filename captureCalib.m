% compute calibration patterns, project them onto monitor screen, and save
% multiple frames for each pattern

clc
close all
clear

% change working directory
mdir = mfilename('fullpath');
mfile = mfilename;
mdir = mdir(1:end-length(mfile));
cd(mdir);

addpath misc/

% capture switch (debug)
captureSwitch = true;

% capture parameters
numberOfFrames = 10;
shutterTime = 0.004;
fps = 20;
delayTimeDisplay = 0.2;
%
% (parameters not currently used...)
numSet = 1; % number of sets for separate calibration
indxSet = 1; % indext of the set
%
delayTimeCapture = numberOfFrames*(shutterTime+1/fps);

% calibration sequence parameters
calibrationSequence = 'hadamard';
screenPatchSize = 128;
screenRes = [screenPatchSize screenPatchSize];
monitoScalingFactor = 1;
maskType = '08R10';

% save directory information
tempPath = 'Data/TEMP/';
testPath = 'testImages/';
dirPath = sprintf('Data/mask%s_%s%d/',maskType,calibrationSequence,screenPatchSize);

msg1 = sprintf('captureCalib planar: calibration-%s, screenPatchSize=%d, maskType=%s',...
    calibrationSequence, screenPatchSize,maskType);
fprintf([msg1,'\n']);

% set up figure etc...
mouseLocation = get(0,'PointerLocation');

% make fullscreen image on the left monitor and resize the position
fullscreenFig = figure('units','pixels','outerposition',[-500 200 200 200],'Color','black');
WindowAPI(f,'Position','full');
axes('Units','normalized','position',...
    [(1-monitoScalingFactor)/2 (1-monitoScalingFactor)/2 monitoScalingFactor monitoScalingFactor])
imagesc(double(imread('cameraman.tif'))/255)
axis('off','square');
colormap(gray); drawnow;

N = prod(screenRes);
% impulseVec = [1; zeros(N-1,1)];
fprintf('computing calibration sequences ...');
switch lower(calibrationSequence)
    case 'hadamard'
        % calibMatrix2 = fwht(eye(N)*N,N,'hadamard');
        opFWHT = @(z) fwht(fwht(z,screenPatchSize,'hadamard')',screenPatchSize,'hadamard')';
        tmpMat = reshape(eye(N),screenPatchSize,screenPatchSize,screenPatchSize^2);
        vec = @(z) z(:);
        calibMatrix = zeros(N);
        % calibMatrix = sparse(N,N);
        for i=1:N
            calibMatrix(:,i) = vec(opFWHT(tmpMat(:,:,i)*screenPatchSize));
            % figure(1); imagesc([reshape(calibMatrix(:,i), screenPatchSize,[]) reshape(calibMatrix2(:,i)/screenPatchSize, screenPatchSize,[])]);
        end
        tmpMat = [];
        % opt_H = @(i) reshape(fwht(circshift(impulseVec,i-1)*N,N,'hadamard'),screenRes);
    case 'identity'
        calibMatrix = eye(N);
    case 'dct'
        calibMatrix = dct(eye(N));
    case 'dct2'
        % calibMatrix2 = dct(eye(N));
        opDCT = @(z) dct(dct(z)')';
        tmpMat = reshape(eye(N),screenPatchSize,screenPatchSize,screenPatchSize^2);
        vec = @(z) z(:);
        calibMatrix = zeros(N);
        for i=1:N
            calibMatrix(:,i) = vec(opDCT(tmpMat(:,:,i)));
            % figure(1); imagesc([reshape(calibMatrix(:,i), screenPatchSize,[]) reshape(calibMatrix2(:,i), screenPatchSize,[])]);
        end
        tmpMat = [];
        %         tmpMat = reshape(eye(N),screenPatchSize,screenPatchSize,screenPatchSize^2);
        %         vec = @(z) z(:);
        %         calibMatrix = zeros(N);
        %         for i=1:N
        %             calibMatrix(:,i) = vec(dct2(tmpMat(:,:,i)));
        %         end
end
fprintf('done ! \n');
calibMatrix = calibMatrix/max(abs(calibMatrix(:)));
% (NOT a good idea...) store calibMatrix in 8 or 16 bit signed format
% nBits = 16;
% eval(sprintf('calibMatrix = int%d(calibMatrix*2^%d-1);',nBits,nBits-1));
calibPM = nnz(find(calibMatrix<0));

if captureSwitch
    % save calibration parameters...
    mkdir(dirPath);
    save(sprintf('%sCalibParameters.mat',dirPath),'-v7.3');
end
 
opH = @(i) reshape(calibMatrix(:,i),screenRes);


%% Capture calibration patterns 
totalTime = 2*(N/numSet)*(delayTimeCapture+2*delayTimeDisplay);

if captureSwitch
    opCapture = @(nameImg) captureImage(mouseLocation,delayTimeDisplay,delayTimeCapture,numberOfFrames,nameImg,dirPath,tempPath);
    % sendmail_lensless({'msalmanasif@gmail.com'},'lensless planar calib experiment started', sprintf('%s \n Expected completion time: %s', msg1,datestr(now+totalTime/(60*60*24))));
    sendmail_lensless({'msalmanasif@gmail.com','ccd-l@mailman.rice.edu'},'planar calib experiment started', sprintf('%s \n expected completion time: %s', msg1,datestr(now+totalTime/(60*60*24))));
    % sendmail_lensless({'msalmanasif@gmail.com','ccd-l@mailman.rice.edu'},'test matlab notification', 'matlab will send such emails to notify the start and end of my data collection... ');
else
    opCapture = @(nameImg) [];
end

screenImage = ones(screenRes);
imagesc(screenImage);
axis('off','square');caxis([0 1]); drawnow;
opCapture('one');

tic;
k = 0;

try
    fprintf('total calibration images = %d \n',N);
    msg_old = [];
    for i = (indxSet-1)*(N/numSet)+1:indxSet*(N/numSet);
        
        % record zero image for each column
        if mod(i,screenPatchSize)==1
            screenImage = zeros(screenRes);
            imagesc(screenImage);
            axis('off','square');caxis([0 1]); drawnow;
            nameImg = sprintf('zero_%02d', i);
            opCapture(nameImg);
        end
        
        % for i = 1:(N/numSet);
        screenImage = opH(i);
        if calibPM
            screenImageP = screenImage;
            screenImageP(screenImageP<0) = 0;
            screenImageN = -screenImage;
            screenImageN(screenImageN<0) = 0;
            
            imagesc(screenImageP);
            axis('off','square');caxis([0 1]); drawnow;
            nameImg = sprintf('P_%02d', i);
            opCapture(nameImg);
            %
            % k = 2*i-1; clc;
            % sprintf('progress %3.0f %%, %3.0f mins remaining ...', 100*k/(2*N/numSet), ceil(totalTime*(1-k/(2*N/numSet))/60))
            
            imagesc(screenImageN);
            axis('off','square');caxis([0 1]); drawnow;
            nameImg = sprintf('N_%02d', i);
            opCapture(nameImg);
            %
            % k = 2*i; clc;
            % sprintf('progress %3.0f %%, %3.0f mins remaining ...', 100*k/(2*N/numSet), ceil(totalTime*(1-k/(2*N/numSet))/60))
            
        else
            imagesc(screenImage);
            axis('off','square');caxis([0 1]); drawnow;
            opCapture(sprintf('I_%02d', i));
            %
            % k = 2*i-1; clc;
            % sprintf('progress %3.0f %%, %3.0f mins remaining ...', 100*k/(2*N/numSet), ceil(totalTime*(1-k/(2*N/numSet))/60))
        end
        msg = sprintf('running ... %d/%d',i,N);
        fprintf([repmat('\b',1,numel(msg_old)),msg]);
        msg_old = msg; 
        if toc > 15*60 && captureSwitch
            sendmail_lensless({'msalmanasif@gmail.com'},'lensless planar calib experiment progress...', sprintf('%s \n %s ',msg1,msg));
            tic; 
        end
    end
    fprintf('\n');
    t = toc;
catch err
    sendmail_lensless({'msalmanasif@gmail.com'},'lensless planar calib experiment error', sprintf('%s \n %s \n %s',msg1,err.identifier,err.message));
    % sendmail_lensless({'msalmanasif@gmail.com','ccd-l@mailman.rice.edu'},'planar calib experiment error', sprintf('%s \n %s \n %s',msg1,err.identifier,err.message));
    rethrow(err);
end

%% Test data for validation

% % make fullscreen image on the left monitor and resize the position
% fullscreenFig = figure('units','pixels','outerposition',[-500 200 200 200],'Color','black');
% maximize(fullscreenFig);
% axes('Units','normalized','position',...
%     [(1-monitoScalingFactor)/2 (1-monitoScalingFactor)/2 monitoScalingFactor monitoScalingFactor])
% imagesc(double(imread('cameraman.tif'))/255)
% axis('off','square');
% colormap(gray);

testDir = dir(testPath);
testDir = testDir(3:end);

fprintf('total test images = %d \n',length(testDir));
msg_old = [];

for i = 1:length(testDir)
    [pathImage, testImage, extImage] = fileparts(testDir(i).name);
    
    % low resolution
    screenImage = im2double(imresize(imread([testPath testImage extImage]),screenRes,'box'));
    imagesc(screenImage);
    axis('off','square');caxis([0 1]); drawnow;
    nameImg = sprintf(testImage);
    opCapture(nameImg);
    
    % high resolution
    screenImage = im2double(imread([testPath testImage extImage]));
    imagesc(screenImage);
    axis('off','square');caxis([0 1]); drawnow;
    nameImg = sprintf('highres_%s',testImage);
    opCapture(nameImg);
    
    % zero image
    screenImage = zeros(screenRes);
    imagesc(screenImage);
    axis('off','square');caxis([0 1]); drawnow;
    nameImg = sprintf('testZero_%02d',i);
    opCapture(nameImg);
    
    msg = sprintf('running ... %d/%d',i,length(testDir));
    fprintf([repmat('\b',1,numel(msg_old)),msg]);
    msg_old = msg;
    
end
fprintf('\n');

%%

close(fullscreenFig)

if captureSwitch
    % sendmail_lensless({'msalmanasif@gmail.com'},'lensless planar calib experiment finished', msg1);
    sendmail_lensless({'msalmanasif@gmail.com','ccd-l@mailman.rice.edu'},'planar calib experiment finished', msg1);
end 