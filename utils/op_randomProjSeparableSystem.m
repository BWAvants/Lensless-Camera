function out = op_randomProjSeparableSystem(in,mode, Lr, Rr, varargin) 

global m n T maxrank

if sum(mode == 1) || strcmp(mode,'notransp'); 
    A_cube = reshape(Lr,m,n,maxrank);
    B_cube = reshape(Rr,m,n,maxrank);
    
    Y = zeros(m^2,T);
    f_h = @(z) op_separableSystem(z, A_cube, B_cube,n,[m m]);
    X = reshape(in,n^2,T);
    for t = 1:T
    	Y(:,t) = f_h(X(:,t));
    end
    out = Y(:); 
    
elseif sum(mode == 0) 
   % y = ft_h(f_h(x));
else
    A_cube = reshape(Lr,m,n,maxrank);
    B_cube = reshape(Rr,m,n,maxrank);
    
    X = reshape(varargin{1},n,n,T);
    yhat = reshape(in,m,m,T);
    
    adjoint_times_L = zeros(m,n,maxrank);
    adjoint_times_R = zeros(m,n,maxrank);
    for r = 1:maxrank
        for t = 1:T
            adjoint_times_R(:,:,r) = adjoint_times_R(:,:,r)+yhat(:,:,t)*B_cube(:,:,r)*X(:,:,t)';
            adjoint_times_L(:,:,r) = adjoint_times_L(:,:,r)+yhat(:,:,t)'*A_cube(:,:,r)*X(:,:,t);
        end
    end    
    adjoint_times_R = reshape(adjoint_times_R,[],maxrank);
    adjoint_times_L = reshape(adjoint_times_L,[],maxrank);
    GradL = 2*(Lr - adjoint_times_R); 
    GradR = 2*(Rr - adjoint_times_L);
    out = [GradL(:);GradR(:)];
end