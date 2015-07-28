% simulated 
fprintf('estimating noise-free simulated sensor images ... ');

calibMat = matfile(calibCompute_mat);
% sensorSize = calibMat.sensorSize;
% sensorSize(3) = 3; 
% numOfChannels = 3; 

%% setup inverse problem
script_inverseProblemSetup;

%%
% simulated measurements and reconstruction for test images
scene_test = imresize(double(imread('peppers.png')),[screenPatchSize,screenPatchSize])/255;
sensor_test = op_Phi_rec(scene_test(:));
    
scene_test_rec = iPhi_rec(sensor_test);

% normalize pixel values
% scene_test_rec = bsxfun(@minus, scene_test_rec, min(scene_test_rec,[],1));
% scene_test_rec = bsxfun(@rdivide, scene_test_rec, max(abs(scene_test_rec),[],1));

f881=figure(881); imshow([reshape(scene_test,screenPatchSize,screenPatchSize,numOfChannels); reshape(scene_test_rec,screenPatchSize,screenPatchSize,numOfChannels)]);
title(sprintf([msg1, '\n simulated test images']),'FontSize',24);
colormap gray; axis image; axis off;
%
% maximize(f881);
set(gcf,'Color','w');
% if saveFigs
%     figName = sprintf('recovery_calib%s_%dx%d_%s_DS%02d%s_simulatedImages', ...
%         calibrationSequence, screenPatchSize, screenPatchSize,maskType,downSamplingFactor,downSamplingMode);
%     export_fig(f881,[testDir,figName],'-png','-pdf','-native');
%     saveas(f881,[testDir,figName],'fig');
% end