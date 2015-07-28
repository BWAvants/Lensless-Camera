% real images
tmp_mat = matfile(sprintf('%srealImages_downSampled%02d_%s.mat',testDir,downSamplingFactor,downSamplingMode));
sensorTest = double(tmp_mat.sensorTest);

% sensorSize = tmp_mat.sensorSize;
% numOfChannels = tmp_mat.numOfChannels; 

%%
fprintf('estimating real images ... ');

% test images
if saveFigs
    test_indices = find(sum(abs(sensorTest),1)>0);
else
    test_indices = find(sum(abs(sensorTest),1)>0);
end
% test_indices = find(sum(abs(sensorTest),1)>0);
% test_indices = [13 14 15]; %[2 3 4 5 6 7 9 10 11 12 16];

sensor_test = sensorTest(:,test_indices);
%
% remove faulty pixels
remove_pixels = [];
sensor_test(remove_pixels,:) = [];

% simulated measurements and reconstruction for test images
scene_test_rec = zeros(screenPatchSize^2*numOfChannels,length(test_indices));
for tst = 1:length(test_indices)
    It = reshape(sensor_test(:,tst),sensorSize);
    It = op_makeSeparable(It);
    sensor_test(:,tst) = It(:);
    scene_test_rec(:,tst) = iPhi_rec(sensor_test(:,tst));
    
    f881=figure(881); imagesc([reshape(scene_test_rec(:,tst),screenPatchSize,[])],[0 1]);
    title(sprintf([msg1, '\n real test images']),'FontSize',24); 
    colormap gray; axis image; axis off; drawnow
end
% add white balance correction
% scene_test_rec = bsxfun(@times,reshape(scene_test_rec, screenPatchSize, screenPatchSize, numOfChannels,[]), permute(1./white_balance_correction,[3 1 2]));
% scene_test_rec = reshape(scene_test_rec,[],length(test_indices)); 
% normalize pixel values
% scene_test_rec = bsxfun(@minus, scene_test_rec, min(scene_test_rec,[],1));
% scene_test_rec = bsxfun(@rdivide, scene_test_rec, max(abs(scene_test_rec),[],1));

vec2rgb = @(z) reshape(permute(reshape(z,screenPatchSize,screenPatchSize,numOfChannels,[]),[1 2 4 3]),screenPatchSize,[],numOfChannels);
scene_rec_rgb = vec2rgb(scene_test_rec(:,:)); 

f882=figure(885); imshow([scene_rec_rgb]);
title(sprintf([msg1, '\n real test images']),'FontSize',24);
colormap gray; axis image; axis off;
%
% maximize(f881);
set(gcf,'Color','w');
if saveFigs
    figName = sprintf('recovery_calib%s_%dx%d_%s_DS%02d%s_realImages', ...
        calibrationSequence, screenPatchSize, screenPatchSize,maskType,downSamplingFactor,downSamplingMode);
    export_fig(f881,[testDir,figName],'-png','-pdf','-native');
    saveas(f881,[testDir,figName],'fig');
end