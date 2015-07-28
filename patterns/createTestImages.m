% create test images for capture
clc
close all
clear

% change working directory
mdir = mfilename('fullpath');
mfile = mfilename;
mdir = mdir(1:end-length(mfile));
cd(mdir);


screenPatchSize = 128;
screenRes = [screenPatchSize screenPatchSize];
hiRes = 512;
numOfChannels = 3;

testPath = 'testImages/';
testDir = dir(testPath);
testDir = testDir(3:end);

fprintf('total test images = %d \n',length(testDir));
msg_old = [];

lowresImages = zeros(prod(screenRes)*numOfChannels,length(testDir)+1);
highresImages = zeros(hiRes^2*numOfChannels,length(testDir)+1);
vec = @(z) z(:);

for i = 1:length(testDir)
    [pathImage, testImage, extImage] = fileparts(testDir(i).name);
    
    % low resolution
    screenImage = imread([testPath testImage extImage]);
    lowresImages(:,i) = vec(imresize(screenImage,screenRes,'box'));
    
    % high resolution
    highresImages(:,i) = vec(imresize(screenImage,[hiRes hiRes],'box'));
    
    msg = sprintf('running ... %d/%d',i,length(testDir));
    fprintf([repmat('\b',1,numel(msg_old)),msg]);
    msg_old = msg;
    
end
% zero images
lowresImages(:,i+1) = vec(zeros([screenRes numOfChannels]));
highresImages(:,i+1) = vec(zeros([hiRes hiRes numOfChannels]));

if numOfChannels==1;
    imwrite(uint8(lowresImages),sprintf('lowresTest%d.bmp',screenPatchSize));
    imwrite(uint8(highresImages),sprintf('highresTest%d.bmp',hiRes));
else
    imwrite(uint8(lowresImages),sprintf('lowresTestC%d.bmp',screenPatchSize));
    imwrite(uint8(highresImages),sprintf('highresTestC%d.bmp',hiRes));
end
fprintf('\n');

