remove_rows = [1:8]; remove_columns = [1:8]; % timestamp pixels in some of the images...
% % remove hot/dark pixels
% experimentID = sprintf('mask%s_%s%d_downSamplingFactor%02d',maskType,calibrationSequence,screenPatchSize,downSamplingFactor); 
% switch experimentID
%     case 'mask16R10_hadamard32_downSamplingFactor04'
%         remove_rows = [1:8]; remove_columns = [1:8]; % 32x32_16R10 test
%     case 'mask16R10_hadamard16_downSamplingFactor08'
%         remove_rows = [1:8]; remove_columns = [1:8]; % 32x32_16R10 test
%     case 'mask16R10_hadamard64_downSamplingFactor04'
%         remove_rows = [1:8]; remove_columns = [1:8]; % 64x64_16R10 test
%     case 'mask16R10_hadamard64_downSamplingFactor08'
%         remove_rows = [1:8]; remove_columns = [1:8]; % 64x64_16R10 test
%     
%     otherwise
%         figure(420); 
%         tmpImages = reshape(sensorImageP(:,1:2:5),[sensorSize,3]);
%         subplot(121); imagesc(reshape(tmpImages,sensorSize(1),[]))
%         subplot(122); imagesc(reshape(permute(tmpImages,[1 3 2]),3*sensorSize(1),[]))
%         error('identify valid rows and columns');
% end

if isempty(remove_rows) || isempty(remove_columns)
    remove_pixels = [];
else
    remove_pixels = sub2ind(sensorSize,vec(remove_rows'*ones(1,numel(remove_columns))),vec(ones(numel(remove_rows),1)*remove_columns));
end
% 
sensorImageP(remove_pixels,:) = [];
sensorImageN(remove_pixels,:) = [];
