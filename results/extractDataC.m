addpath '../misc/';
addpath '../utils/export_fig/';

% calibration sequence parameters
calibrationSequence = 'hadamardSeparable'; % hadamardSeparable 
screenPatchSize = 256; 
maskType = '255M030C';
downSamplingFactor = 4;
downSamplingMode = 'average';
numOfChannels = 3;
%
figName = sprintf('recovery_calib%s_%dx%d_%s_DS%02d%s_highres_5in', ...
    calibrationSequence, screenPatchSize, screenPatchSize,maskType,downSamplingFactor,downSamplingMode);

dirPath = [figName,'/'];
mkdir(dirPath);

h1 = open([figName,'.fig']);
axesObjs = get(h1, 'Children');  %axes handles
dataObjs = get(axesObjs, 'Children'); %handles to low-level graphics objects in axes

objTypes = get(dataObjs, 'Type');  %type of low-level graphics object   
xdata = get(dataObjs, 'XData');  %data from low-level grahics objects
ydata = get(dataObjs, 'YData');
cdata = get(dataObjs, 'CData');
cdata_cube = reshape(cdata,size(cdata,1),screenPatchSize,[],numOfChannels);
cdata_cube = permute(cdata_cube,[1 2 4 3]);
for ii=1:size(cdata_cube,4);
%    h2 = figure(2); imshow(cdata_cube(size(cdata,1)/2+1:end,:,ii));
    h2 = figure(2); imshow(cdata_cube(end-screenPatchSize+1:end,:,:,ii));
%     figName2 = sprintf('recovery_%dx%d_%s_DS%02d%s_highres_%02d', ...
%     screenPatchSize, screenPatchSize, maskType, downSamplingFactor, downSamplingMode,ii);
    figName2 = sprintf('%s%02d',dirPath,ii);
    
    export_fig(h2, figName2,'-png','-native'); 
    
    if screenPatchSize == 1024
        h3 = figure(3); imshow(cdata_cube(1:size(cdata,1)/2,:,ii));
        figName3 = sprintf('highres_%02d',ii);
        export_fig(h3, figName3,'-png','-native'); 
    end

end