clear; clc; close all

% addpath('misc');
% change working directory
mdir = mfilename('fullpath');
mfile = mfilename;
mdir = mdir(1:end-length(mfile));
cd(mdir);
cd ..
addpath('utils');

monitoScalingFactor = 1;
N = 256;

pattern = 'shifts';


% make fullscreen image on the left monitor and resize the position
fullscreenFig = figure('units','pixels','outerposition',[-500 200 200 200],'Color','black');
WindowAPI(fullscreenFig,'Position','full');

vec = @(z) z(:);

switch pattern
    case 'hadamard' 
        for ii = 1:N^2
            x = fwht(circshift([1;zeros(N^2-1,1)],ii-1)*N^2,N^2,'hadamard');
           
            % imagesc(double(imread('cameraman.tif'))/255)
            imagesc(reshape(x,N,N));
            axis('off','square'); caxis([0 1]);
            colormap(gray);
            pause
        end
    case 'hadamardSeparable'
        for ii = 1:N
            x = fwht(circshift([1;zeros(N-1,1)],ii-1)*N,N,'hadamard');
            x = repmat(x>0,1,N)';
            % imagesc(double(imread('cameraman.tif'))/255)
            imagesc(reshape(x,N,N));
            axis('off','square'); caxis([0 1]);
            colormap(gray);
            pause
        end
    case 'dct'
        tmpMat = reshape(eye(N^2),N,N,N^2);
        vec = @(z) z(:);
        calibMatrix = zeros(N^2);
        for i=1:N^2
            calibMatrix(:,i) = vec(dct2(tmpMat(:,:,i)));
        end
        calibMatrix = calibMatrix/max(abs(calibMatrix(:)));
        for ii =1:N^2
            x = (dct((dct(reshape(circshift([1;zeros(N^2-1,1)],ii-1),N,N)))'))'*N^2;
            x1 = calibMatrix(:,ii);
            % max(abs(x(:)-x1(:)))
            
            imagesc(reshape(x,N,N));
            axis('off','square'); caxis([0 1]);
            colormap(gray);
            pause
        end
    case 'shifts'
        fprintf('\n');
        x1 = circshift([1;zeros(N-1,1)],N/2-1);
        x2 = circshift([1;0; 1; zeros(N-3,1)],N/2-1);
        r = 1; N = N*r; 
        for ii = N/r/2:(N/r)%N/r/2:(N/r)
            fprintf('%d.',ii);
        % for ii = -5:5
            % eval(sprintf('x=x%d;',mod(ii,2)+1));
            
            x = circshift([ones(r,1); zeros(N-r,1)],r*ii-1);
            x = repmat(x(:),1,N)';
            
            % x = zeros(N);
            % x(N/2-r/2:N/2+r/2,:) = 1; 
            % x = imrotate(x,ii,'crop');

            %             x = circshift([1; 0; 0;zeros(N^2-3,1)],ii-1);
            %             x = ones(N^2,1)*(N/ii);
            %             x = zeros(N^2,1)*(N/ii);
            
            % imagesc(double(imread('cameraman.tif'))/255)
            % imagesc(reshape(x,N,N));
%             if (~exist(sprintf('patterns/shiftingBars%d',N),'dir'));
%                 mkdir(sprintf('patterns/shiftingBars%d',N));
%             end
%             imwrite(x,sprintf('patterns/shiftingBars%d/imge_%05d.bmp',N,ii));
%             imwrite(x',sprintf('patterns/shiftingBars%d/imge_%05d.bmp',N,N+ii)); 
            imagesc(x');
            axis('off','square'); caxis([0 1]);
            colormap(gray);
            pause
        end
        fprintf('\n');

    case 'testImages'
        testPath = 'testImages/';
        testDir = dir(testPath);
        testDir = testDir(3:end);
        screenRes = [N N];
        
        fprintf('total test images = %d \n running...',length(testDir));
        
        for i = 1:length(testDir)
            [pathImage, testImage, extImage] = fileparts(testDir(i).name);
            % screenImage = im2double(imresize(imread([testPath testImage extImage]),screenRes,'box'));
            screenImage = im2double((imread([testPath testImage extImage])));
            imagesc(screenImage);
            axis('off','square');caxis([0 1]); drawnow;
            fprintf('\b\b\b %02d',i)
            pause
            
            screenImage = zeros(screenRes);
            imagesc(screenImage);
            axis('off','square');caxis([0 1]); drawnow;
            
        end
        fprintf('\n');
end
