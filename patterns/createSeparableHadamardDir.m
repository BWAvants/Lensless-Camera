n = 128; 
% 2D hadamard for nxn image
% H = fwht(eye(n^2)*n^2,n^2,'hadamard'); % reshape
% % separable hadamard
Hrow = fwht(eye(n)*n,n,'hadamard');
% Hsep = [repmat(Hrow,n,1) reshape(repmat(permute(Hrow',[3 2 1]),[n 1 1]),n^2,n)];
% write to bmp images as 0/1

% mkdir(sprintf('hadamard%d',n));
% for i = 1:size(H,2)
%     imwrite(reshape(H(:,i),[n n])>0,sprintf('hadamard%d/pattern_%05d.bmp',n,i));
% end

% Hsep = imread(sprintf('hadamardSeparable%d.bmp',n));
%% separable patterns
dirName = sprintf('hadamardSeparable%d',n);
mkdir(dirName);
for i = 1:size(Hrow,2)
    imwrite(repmat(Hrow(:,i),1,n)>0,sprintf('%s/pattern_%05d.bmp',dirName,i));
end
for i = 1+size(Hrow,2):2*size(Hrow,2)
    imwrite(repmat(Hrow(:,i-size(Hrow,2))',n,1)>0,sprintf('%s/pattern_%05d.bmp',dirName,i));
end

break; 
%% rotated monitor pattern
dirName = sprintf('hadamardSeparable%d_rotated',n);
mkdir(dirName);
for i = 1:size(Hrow,2)
    Hsep = repmat(Hrow(:,i),1,n)>0;
    imwrite(Hsep',sprintf('%s/pattern_%05d.bmp',dirName,i));
    imwrite(Hsep,sprintf('%s/pattern_%05d.bmp',dirName,i+size(Hrow,2)));
end

% 
% % Hsep = imread(sprintf('hadamardSeparable%d.bmp',n));
% mkdir(sprintf('hadamardSeparable%d',n));
% for i = 1:size(Hsep,2)
%     imwrite(reshape(Hsep(:,i),[n n])>0,sprintf('hadamardSeparable%d/pattern_%05d.bmp',n,i));
% end