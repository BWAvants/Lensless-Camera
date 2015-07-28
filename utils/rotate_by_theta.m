function img_out = rotate_by_theta(img1,maskType,screenPatchSize)
% rotate an image by small angle and crop the central part of the image... 
 
dS_list = []; 
switch maskType
    case '255M030C'
        switch screenPatchSize
            case 128
                rot_angle_list = [-0.3:0.01:0.1]; 
                crop_center = @(z) z(10:end-10,10:end-10,:); 
            case 256
                rot_angle = 0.06; 
        end
    case '127M060C'
end

for rot_angle = rot_angle_list
theta = -pi/180*rot_angle;

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

img_out = crop_center(img2);

%%
figure(222); 
for cc = 1:size(img1,3)
    subplot(131); imagesc([img1(:,:,cc)]); colormap gray; shg;
    subplot(132); imagesc([img_out(:,:,cc)]); colormap gray; shg;
    title(sprintf('i = %3.4g',i));
    [U S V] = svd(subtract_mean(img_out(:,:,cc))); dS = diag(S);
    subplot(133); plot(dS(1:10));
    title(sprintf('S1/S2 = %3.4g, rot. angle = %3.4g', dS(1)/dS(2), rot_angle));
    dS_list = [dS_list; rot_angle dS(1)/dS(2)];
    pause(1/60); 
end
end
[best_dS best_angle] = max(dS_list(:,2));
fprintf('maskType-%s, screenPatchSize=%d, best dS(1)/dS(2) = %3.4g, best rot_angle = %3.4g. \n ',maskType, screenPatchSize, best_dS, dS_list(best_angle)); 
end