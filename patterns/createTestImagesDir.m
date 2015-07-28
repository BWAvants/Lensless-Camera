clc
clear

% n = 512;
% patterTarget = 'highresTest';
% 
% mkdir(sprintf('%s%d',patterTarget,n));
% I = imread(sprintf('%s%d.bmp',patterTarget,n));
% for i = 1:size(I,2)
%     imwrite(reshape(I(:,i),[n n]),sprintf('%s%d/pattern_%05d.bmp',patterTarget,n,i));
% end


n = 128;
dirTarget = 'testImages';
patternTarget = 'lowresTestC';
l = dir(dirTarget); l = l(3:end);
mkdir(sprintf('%s%d',patternTarget,n));
for i = 1:length(l)
    I = imread(sprintf('%s/%s',dirTarget,l(i).name));
    imwrite(imresize(I,[n n],'nearest'),sprintf('%s%d/pattern_%05d.bmp',patternTarget,n,i));
end


% n = 1024;
% dirTarget = 'testImagesC';
% patternTarget = 'highresTestC';
% l = dir(dirTarget); l = l(3:end);
% mkdir(sprintf('%s%d',patternTarget,n));
% for i = 1:length(l)
%     I = imread(sprintf('%s/%s',dirTarget,l(i).name));
%     imwrite(I,sprintf('%s%d/pattern_%05d.bmp',patternTarget,n,i));
% end