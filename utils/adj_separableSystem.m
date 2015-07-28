function x = adj_separableSystem(y, Phi_rec_left, Phi_rec_right,screenPatchSize,sensorSize,varargin) 
% adjoint operator for a multi-rank camera system... 
  
numOfChannels = 1; 
if length(sensorSize) > 2; numOfChannels = 3; end

if nargin > 5 
    lambda = varargin{1}; 
else 
    lambda = 0;
end
if lambda > 0 
    Y = reshape(y(1:prod(sensorSize)),sensorSize);
    X = reshape(lambda*y(prod(sensorSize)+1:end),screenPatchSize,[]);
else
    Y = reshape(y,sensorSize);
    X = zeros(screenPatchSize,screenPatchSize, numOfChannels);
end

X = reshape(X,screenPatchSize,screenPatchSize, numOfChannels); 
for r = 1:size(Phi_rec_left,3)
    for c = 1:numOfChannels
        X(:,:,c) = X(:,:,c)+Phi_rec_left(:,:,r)'*Y(:,:,c)*Phi_rec_right(:,:,r);
    end
end
x = X(:);