n = 128; 
% 2D hadamard for nxn image
% H = fwht(eye(n^2)*n^2,n^2,'hadamard'); % reshape
% separable hadamard
Hrow = fwht(eye(n)*n,n,'hadamard');

vec = @(z) z(:);

Hsep = zeros(n^2,2*n);
for i = 1:n
    Hsep(:,i) = vec(repmat(Hrow(:,i),1,n)>0);
end

for i = 1+n:2*n
    Hsep(:,i) = vec(repmat(Hrow(:,i-n)',n,1)>0);
end


% Hsep = [repmat(Hrow,n,1) reshape(repmat(permute(Hrow',[3 2 1]),[n 1 1]),n^2,n)];
% write to bmp images as 0/1
% imwrite(H>0,sprintf('hadamard%d.bmp',n));
imwrite(Hsep>0,sprintf('hadamardSeparable%d.bmp',n));
