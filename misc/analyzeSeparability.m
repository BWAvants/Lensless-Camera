for row = 1:screenPatchSize;
    UV = Phi_rec_left(:,row)*sum(Phi_rec_right,2)';
    UVt = UV+mean(UV(:)); % mean(UV,2)*ones(1,size(UV,2));
    [Ut St Vt] = svd(UVt); St = diag(St);
    figure(1);
    subplot(121)
    imagesc(UVt);
    axis image;
    title(sprintf('row=%d, col=%d -- mean=%3.4g.',row,col,mean(UVt(:))));
    colormap gray; colorbar; shg;
    subplot(122)
    imagesc([Ut(:,1)*St(1)*Vt(:,1)' Ut(:,2)*St(2)*Vt(:,2)']);
    axis image;
end;
    