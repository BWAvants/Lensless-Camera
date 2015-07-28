clc
close all
clear


% addpath('misc');
% change working directory
mdir = mfilename('fullpath');
mfile = mfilename;
mdir = mdir(1:end-length(mfile));
cd(mdir);
cd ..
addpath('utils');

bar = false;
showPattern = true;
screenPatchSize = 256;
patterns = [3];
if showPattern
    fullscreenFig = figure('units','pixels','outerposition',[-500 200 200 200],'Color','black');
    WindowAPI(fullscreenFig,'Position','full');
    for k = 1:length(patterns);
        I = imread(sprintf('patterns/hadamardSeparable%d/pattern_%05d.bmp',screenPatchSize,patterns(k)));
        % I = -(I-1); 
        % I = I';
        imagesc(I/255)
        axis('off','square');
        colormap(gray); drawnow;
        pause
    end
end

if bar
    I = zeros(screenPatchSize,screenPatchSize);
    I(screenPatchSize/2,:) = 1;
    % I(:,screenPatchSize/2) = 1;
    
    fullscreenFig = figure('units','pixels','outerposition',[-500 200 200 200],'Color','black');
    WindowAPI(fullscreenFig,'Position','full');
    imagesc(I)
    axis('off','square');
    colormap(gray); drawnow;
    pause
end
%%
f101 = figure(111);
while(true)
    dirPath = 'Data/TEMP/';
    downSamplingMode = 'average';
    downSamplingFactor = 4;
    numberOfFrames = 1;
    
    resize_temp = @(z) z;
    switch downSamplingMode
        case 'average'
            opReadResize = @(nameImg) imresize(double(resize_temp(imread([dirPath,nameImg])))/numberOfFrames,1/downSamplingFactor,'box');
        case 'subsample'
            downSample = @(z) z(1:downSamplingFactor:end,1:downSamplingFactor:end,:);
            opReadResize = @(nameImg) downSample(double(resize_temp((imread([dirPath,nameImg]))))/numberOfFrames);
    end
    
    dataDir = dir(dirPath);
    dataDir = dataDir(3:end);
    
    sensorSize = size(imread([dirPath dataDir(1).name]));
    theImage = 0;
    for k = 1:length(dataDir)   
    for i = 0
        theImage1 = (opReadResize(dataDir(k).name));
%         theImage = imrotate(theImage1,i,'crop','bilinear');
        theImage = theImage + theImage1;
        for C = 1
            [U S V] = svd(subtract_mean(mean(theImage(:,:,C),3)));
            dS = diag(S); 
            plot(dS(1:10));
            title(sprintf('%s \n S1/S2 = %3.4g, rot. angle = %3.4g, channel-%d',dataDir(k).name, dS(1)/dS(2), i,C),'interpreter','none'); 
            % title(['S1/S2=' num2str(dS(1)/dS(2))]); 
            drawnow
        end
    end
    end
end

%%
img1 = subtract_mean(theImage1);
for i = 0.1:0.01:0.1
theta = -pi/180*i;
rmat = [
cos(theta) sin(theta) 0
-sin(theta) cos(theta) 0
0           0          1];

mx = size(img1,2);
my = size(img1,1);
corners = [
    0  0  1
    mx 0  1
    0  my 1
    mx my 1];
new_c = corners*rmat;

T = maketform('affine', rmat);   %# represents translation
img2 = imtransform(img1, T, ...
    'XData',[min(new_c(:,1)) max(new_c(:,1))],...
    'YData',[min(new_c(:,2)) max(new_c(:,2))]);
figure(222); 
subplot(121); imagesc([img2(150:175,3:end-1)]); colormap gray; shg;
title(sprintf('i = %3.4g',i)); 
img2_center = img2(10:end-10,10:end-10); 
[U S V] = svd(subtract_mean(img2_center)); dS = diag(S); 
subplot(122); plot(dS(1:10)); 
title(sprintf('S1/S2 = %3.4g, rot. angle = %3.4g', dS(1)/dS(2), i));

end