% test images
%
fprintf('reading sensor test images ... ');
% % lowres test images
tmp_mat = matfile(sprintf('%stestImages_downSampled%02d_%s.mat',testDir,downSamplingFactor,downSamplingMode));
lowrestestImages = double(tmp_mat.lowresImages);
% sensorTest_lowres = tmp_mat.sensorTest_lowres;
% highres test images
tmp_mat = matfile(sprintf('%stestImages_downSampled%02d_%s.mat',testDir,downSamplingFactor,downSamplingMode));
highrestestImages = double(tmp_mat.highresImages);
sensorTest_highres = tmp_mat.sensorTest_highres;
%
fprintf('done! \n');

% sensorSize = tmp_mat.sensorSize;
% numOfChannels = tmp_mat.numOfChannels; 

calibdata_mat = matfile(calibCompute_mat);

if saveFigs
    test_indices = find(sum(abs(sensorTest_highres),1)>0);
    test_indices = [17 21 125 127 145 147];    
else
    test_indices = [15 17 21 125:2:129 143 145 147];
    test_indices = [17 21 125 127 145 147];    
end


%% lowres
if exist('sensorTest_lowres','var')
fprintf('estimating lowres test images (%dx%d image from %dx%d measurements) ... ',screenPatchSize,screenPatchSize,sensorSize(1:2));

% test images
% sensorTest_zero = sensorTest_lowres(:,end-1);
sensorTest = sensorTest_lowres(:,1:end);
screenTest = lowrestestImages(:,1:end);

% sensorDC = sensorZero;
sensorDC = mean(calibdata_mat.sensorZero,2);
sensor_test = sensorTest(:,test_indices);%-repmat(mean(sensorDC,2),1,length(test_indices));
%
% remove faulty pixels
% sensor_test(remove_pixels,:) = [];
%
scene_test = screenTest(:,test_indices)/max(screenTest(:));

% simulated measurements and reconstruction for test images
sensor_test_simulated = 0*sensor_test;
scene_test_rec = zeros(screenPatchSize^2*numOfChannels,length(test_indices));
for tst = 1:length(test_indices)
    It = reshape(sensor_test(:,tst),sensorSize);
    It = op_makeSeparable(It);
    sensor_test(:,tst) = It(:);
    sensor_test_simulated(:,tst) = op_Phi_rec(scene_test(:,tst));
    scene_test_rec(:,tst) = iPhi_rec(sensor_test(:,tst));
    
    f881=figure(881); imagesc([reshape(scene_test(:,tst),screenPatchSize,[]); reshape(scene_test_rec(:,tst)/max(scene_test_rec(:,tst)),screenPatchSize,[])],[0 1]);
    title(sprintf([msg1, '\n lowres test images']),'FontSize',24);
    colormap gray; axis image; axis off; drawnow

end
% % add white balance correction
% scene_test_rec = bsxfun(@times,reshape(scene_test_rec, screenPatchSize, screenPatchSize, numOfChannels,[]), permute(1./white_balance_correction,[3 1 2]));
% scene_test_rec = reshape(scene_test_rec,[],length(test_indices));
% % normalize pixel values
% scene_test_rec = bsxfun(@minus, scene_test_rec, min(scene_test_rec,[],1));
% scene_test_rec = bsxfun(@rdivide, scene_test_rec, max(abs(scene_test_rec),[],1));
        
vec2rgb = @(z) reshape(permute(reshape(z,screenPatchSize,screenPatchSize,numOfChannels,length(test_indices)),[1 2 4 3]),screenPatchSize,[],numOfChannels);

figure(f881); imshow([vec2rgb(scene_test); vec2rgb(scene_test_rec)]);
title(sprintf([msg1, '\n lowres test images']),'FontSize',24);
if numOfChannels == 1; colormap gray; end
axis image; axis off;
%
maximize(f881); set(gcf,'Color','w');
if saveFigs
    %     figName = sprintf('recovery_calib%s_%dx%d_%s_misc_DS%02d%s', ...
    %         calibrationSequence, screenPatchSize, screenPatchSize,maskType,downSamplingFactor,downSamplingMode);
    %     export_fig(f771,[testDir,figName],'-png','-pdf');
    % saveas(f771,[testDir,figName],'fig');
    figName = sprintf('recovery_calib%s_%dx%d_%s_DS%02d%s', ...
        calibrationSequence, screenPatchSize, screenPatchSize,maskType,downSamplingFactor,downSamplingMode);
    export_fig(f881,[testDir,figName],'-png','-pdf');
    saveas(f881,[testDir,figName],'fig');
end
fprintf('done! \n');
end

%% highres
fprintf('estimating highres test images ... ');

% test images
% sensorTest_zero = sensorTest_highres(:,end-1);
sensorTest = sensorTest_highres(:,1:end);
screenTest = lowrestestImages(:,1:end);

% sensorDC = sensorZero;
sensorDC = mean(calibdata_mat.sensorZero,2);
sensor_test = sensorTest(:,test_indices);%-repmat(mean(sensorDC,2),1,length(test_indices));
%
% remove faulty pixels
% sensor_test(remove_pixels,:) = [];
%
scene_test = screenTest(:,test_indices)/max(screenTest(:));

% simulated measurements and reconstruction for test images
sensor_test_simulated = 0*sensor_test;
scene_test_rec = zeros(screenPatchSize^2*numOfChannels,length(test_indices));
for tst = 1:length(test_indices)
    It = reshape(sensor_test(:,tst),sensorSize);
    It = op_makeSeparable(It);
    sensor_test(:,tst) = It(:);
    sensor_test_simulated(:,tst) = op_Phi_rec(scene_test(:,tst));
    scene_test_rec(:,tst) = iPhi_rec(sensor_test(:,tst));
    
    f882=figure(882); imagesc([reshape(scene_test(:,tst),screenPatchSize,[]); reshape(scene_test_rec(:,tst),screenPatchSize,[])],[0 1]);
    title(sprintf([msg1, '\n highres test images']),'FontSize',24);
    colormap gray; axis image; axis off; drawnow
end
% add white balance correction
% scene_test_rec = bsxfun(@times,reshape(scene_test_rec, screenPatchSize, screenPatchSize, numOfChannels,[]), permute(1./white_balance_correction,[3 1 2]));
% scene_test_rec = reshape(scene_test_rec,[],length(test_indices));
% normalize pixel values
% scene_test_rec = bsxfun(@minus, scene_test_rec, min(scene_test_rec,[],1));
% scene_test_rec = bsxfun(@rdivide, scene_test_rec, max(abs(scene_test_rec),[],1));
        
vec2rgb = @(z) reshape(permute(reshape(z,screenPatchSize,screenPatchSize,numOfChannels,[]),[1 2 4 3]),screenPatchSize,[],numOfChannels);
scene_rgb = vec2rgb(scene_test(:,:)); 
scene_rec_rgb = vec2rgb(scene_test_rec(:,:)); 

figure(f882); imshow([scene_rgb; scene_rec_rgb]);
title(sprintf([msg1, '\n highres test images']),'FontSize',24);
if numOfChannels == 1; colormap gray; end
axis image; axis off;
%
maximize(f882); set(gcf,'Color','w');

if saveFigs
    %     figName = sprintf('recovery_calib%s_%dx%d_%s_misc_DS%02d%s_highres', ...
    %         calibrationSequence, screenPatchSize, screenPatchSize,maskType,downSamplingFactor,downSamplingMode);
    %     export_fig(f772, [testDir,figName],'-png','-pdf');
    % saveas(f772,[testDir,figName],'fig');
    figName = sprintf('recovery_calib%s_%dx%d_%s_DS%02d%s_highres', ...
        calibrationSequence, screenPatchSize, screenPatchSize,maskType,downSamplingFactor,downSamplingMode);
    export_fig(f882, [testDir,figName],'-png','-pdf');
    saveas(f882,[testDir,figName],'fig');
end
fprintf('done! \n');

break; 
%% compare image projections onto a set of left/rigth singular vectors
[Ul Sl Vl] = svd(Phi_rec_left,'econ'); Sl = diag(Sl);
[Ur Sr Vr] = svd(Phi_rec_right','econ'); Sr = diag(Sr);
Tl=0+[1:256]; Tr = 0+[1:256];
left_proj = @(z) reshape(Vl(:,Tl)*Vl(:,Tl)'*reshape(z,screenPatchSize,[]),screenPatchSize,screenPatchSize,[]);
right_proj = @(z) permute(reshape(reshape(permute(z,[1 3 2]),[],screenPatchSize)*Ur(:,Tr)*Ur(:,Tr)',screenPatchSize,[],screenPatchSize),[1 3 2]);
Itest = reshape(scene_test(:,5),screenPatchSize,screenPatchSize,[]);
figure(1234); imshow([Itest right_proj(left_proj(Itest))]);
