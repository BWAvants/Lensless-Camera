%% normalize the columns of left/right matrices to remove the unwanted bright/dark lines in the reconstruction
% I think those lines appear because the columns of Phi_rec_left and
% Phi_rec_right are not normalized properly. 
% e.g., tmp = sum(Phi_rec_left.^2,1); will show spikes in exactly those
% rows that are bright. 
% median filter
med_width = 3;
mfilt = @(z) medfilt1(z,med_width);
% polynomial fit
porder = 16;
pfilt = @(z) polyval(polyfit([1:length(z)]-length(z)/2,z,porder),[1:length(z)]-length(z)/2);

switch screenPatchSize
    case {128}
        nfiltType = 'median';
    case {256,512}
        nfiltType = 'poly';
end
switch nfiltType
    case 'median'
        nfilt = @(z) mfilt(z);
        fprintf('normalizing columns of Phi_rec_left/right using median filter of width %d ...\n',med_width);
    case 'poly'
        nfilt = @(z) pfilt((z)); 
        % figure(1); plot([tmp_left(:) nfilt(tmp_left)])
        fprintf('normalizing columns of Phi_rec_left/right using polyfit with degree %d polynomial...\n',porder);
end

tmp_left = sqrt(sum(Phi_rec_left.^2,1)); 
tmp_left_filt = nfilt(tmp_left);
% normalization_left = (tmp_left_filt./tmp_left);
% normalization_left = (1./tmp_left);
tmp_right = sqrt(sum(Phi_rec_right.^2,1)); 
tmp_right_filt = nfilt(tmp_right); 
% normalization_right = (tmp_right_filt./tmp_right);
% normalization_right = (1./tmp_right); 

% model the intensity decay away from center using a disk...
aperture = imresize(fspecial('disk', screenPatchSize/2), [screenPatchSize,screenPatchSize]);
aperture = aperture/max(aperture(:));
aperture = ones(size(aperture));

includeAperture = 2; 
switch includeAperture
    case 0
        % spatial window (cosine-falloff) and mask patterns kept separate
        Phi_rec_left = Phi_rec_left*diag(1./tmp_left);
        Phi_rec_right = Phi_rec_right*diag(1./tmp_right);
        aperture = (tmp_left_filt(:)*tmp_right_filt(:)');

    case 1
        % spatial window (cosine-falloff) and the mask patterns combined
        Phi_rec_left = Phi_rec_left*diag(tmp_left_filt./tmp_left);
        Phi_rec_right = Phi_rec_right*diag(tmp_right_filt./tmp_right);        
end
aperture = repmat(aperture, [1 1 numOfChannels]); 
