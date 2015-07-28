%% setup inverse problem
vec = @(z) z(:);
calibMat = matfile(calibCompute_mat);

Phi_rec_left = calibMat.Phi_rec_left;
Phi_rec_right = calibMat.Phi_rec_right;
% white_balance_correction = calibMat.white_balance_correction;
% screenImage_left = calibMat.screenImage_left;
% screenImage_right = calibMat.screenImage_right;
sensorSize = calibMat.sensorSize;
numOfChannels = calibMat.numOfChannels;

fprintf('Calib matrices rank = %d \n',size(Phi_rec_left,3));

% normalize columns
normalize_Phi_columns;


fprintf('solving inverse problem using %s ...  ',solver);

vec = @(z) z(:);
rshp2D_sensor = @(z) reshape(z,sensorSize);
rshp2D_scene = @(z) reshape(z,screenPatchSize,screenPatchSize,numOfChannels);

% measurement operators
op_Phi_rec = @(z) op_separableSystem(aperture(:).*z,1, Phi_rec_left, Phi_rec_right, screenPatchSize,sensorSize);
A_h = @(z) op_separableSystem(z,1, Phi_rec_left, Phi_rec_right, screenPatchSize,sensorSize);
At_h = @(z) op_separableSystem(z, 2, Phi_rec_left, Phi_rec_right, screenPatchSize,sensorSize);

switch solver
    case {'lasso'}
        % Wavelet parameters
        global QMF waveletLevel
        waveletType = 'Daubechies';
        waveletPar = 8;
        waveletLevel = 5;
        QMF = MakeONFilter(waveletType,waveletPar);
        psiT = @(z) vec(FWT2_PO(rshp2D_scene(z), waveletLevel, QMF));
        psi = @(z) vec(IWT2_PO(rshp2D_scene(z),waveletLevel,QMF));
end

switch solver
    case 'svd'
        % svd for left- and right-side matrices
        [Ul Sl Vl] = svd(Phi_rec_left,'econ'); Sl = diag(Sl);
        [Ur Sr Vr] = svd(Phi_rec_right','econ'); Sr = diag(Sr);
        
        T = round(screenPatchSize/2); lambda = 1;
        iSl = speye(min(screenPatchSize,size(Ul,2))); iSl = spdiags(1./(Sl+lambda),0,iSl);
        iSr = speye(min(screenPatchSize,size(Vr,2))); iSr = spdiags(1./(Sr+lambda),0,iSr);
        % applying left/right inverse matrices on rgb channels in a single operation
        iPhi_rec_left = @(z) reshape(Vl(:,1:T)*(iSl(1:T,1:T)*(Ul(:,1:T)'*reshape(z,sensorSize(1),[]))),screenPatchSize,sensorSize(2),[]);
        iPhi_rec_right = @(z) permute(reshape([reshape(permute(z,[1 3 2]),[],sensorSize(2))*Vr(:,1:T)*iSr(1:T,1:T)*Ur(:,1:T)'],screenPatchSize,[],screenPatchSize),[1 3 2]);
        % including white-balance correction
        iPhi_rec = @(z) vec(iPhi_rec_right(iPhi_rec_left(z)));
         
    case 'lsqr'
        % lambda = sqrt(1);
        % f_h = @(z) op_separableSystem(z, Phi_rec_left, Phi_rec_right, screenPatchSize,sensorSize,lambda);
        % ft_h = @(z) adj_separableSystem(z, Phi_rec_left, Phi_rec_right, screenPatchSize,sensorSize,lambda);
        %
        % fhandle = @(z,mode) fhandle_mode(z,mode, f_h, ft_h);
        % tol = 1e-6;
        % maxit = 1000;
        % iPhi_rec = @(z) lsqr(fhandle,[z;zeros(screenPatchSize^2*(lambda>0),1)],tol,maxit);
        
        
        damp = 1;
        atol = 1e-6;
        btol = 1e-6;
        conlim = 1e6;
        itnlim = 1000;
        show = 0; 
        f_h = @(z) A_h(aperture(:).*z); ft_h = @(z) aperture(:).*(At_h(z));
        fhandle = @(z,mode) fhandle_mode(z,mode, f_h, ft_h);
        iPhi_rec = @(z) lsqrSOL(numel(z), screenPatchSize^2*numOfChannels, fhandle, z, damp, atol, btol, conlim, itnlim, show );
   
    case 'tvdn'
        param.verbose = 1;  % Print log or not
        param.gamma = 1e-4; % Converge parameter
        param.tol = 1e-4;   % Stopping criterion for the TVDN problem
        param.maxit = 200;  % Max. number of iterations for the TVDN problem
        % param.maxit_tv = 100; % Max. nb. of iter. for the sub-problem (proximal TV operator)
        param.nu_b2 = 1e3;    % Bound on the norm of the operator A
        param.tol_b2 = 1e-4; % Tolerance for the projection onto the L2-ball
        param.tight_b2 = 0; % Indicate if A is a tight frame (1) or not (0)
        % param.maxit_b2 = 500;
        
        SIGMA = 100;%norm(vec(Y-A*sceneImgRecPI*A'));
        optTVDN = @(z) rshp2D_sensor(A_h(z));
        optTVDNt = @(z) rshp2D_scene(At_h(z)); 
        
        iPhi_rec = @(z) vec(solve_tvdn(rshp2D_sensor(z),SIGMA,optTVDN,optTVDNt,param));
        
    case 'tval3'
        
        opts = []; 
        opts.mu = 2^10;              % [2^8] (primary penalty parameter)
        opts.beta = 2^4;            % [2^5] (secondary penalty parameter)
        % opts.mu0 = 2^4;             % opts.mu (initial mu for continuation) % trigger continuation shceme
        % opts.beta0 = 2^4;           % opts.beta (initial beta for continuation) % trigger continuation shceme
        opts.tol = 1e-6;            % [1.e-6] (outer stopping tolerance)
        opts.tol_inn = 1e-6;        % [1.e-3] (inner stopping tolerance)
        opts.maxit = 1000;          % [1025] (maximum total iterations)
        opts.maxcnt = 20;           % [10] (maximum outer iterations)
        opts.TVnorm = 2;            % [2] (isotropic or anisotropic TV)
        opts.nonneg = true;         % [false] (switch for nonnegative models)
        opts.TVL2 = true;           % [false] (switch for TV/L2 models)
        opts.isreal = true;         % [false] (switch for real signals/images)
        % opts.scale_A = true;       % [true] (switch for scaling A)
        % opts.scale_b = true;       % [true] (switch for scaling b)
        opts.disp = false;          % [false] (switch for iteration info printout)
        opts.init = 0;              % [1] (initial guess)
 
        f_h = @(z) A_h(aperture(:).*z); ft_h = @(z) aperture(:).*(At_h(z));
        fhandle = @(z,mode) fhandle_mode(z,mode, f_h, ft_h);
        iPhi_rec = @(z) vec(TVAL3(fhandle,z(:),screenPatchSize,screenPatchSize,opts)); 
        
        %     case 'tv-unlock'
        %         epsilon = 3000;
        %         param.maxit = 6;
        %         iPhi_rec = @(z) vec(solve_tvdn(reshape(sensor_test(:,tst),sensorSize), epsilon, op_Phi_rec_TV, opt_Phi_rec_TV,param));
        %

    case 'bpdn'
        opts = [];
        opts.tol = 1e-5;  % Use super strict tolerance
        opts.recordObjective = true; %  Record the objective function so we can plot it
        opts.verbosity = 0;
        opts.iterations = 1000;
        f_h = @(z) A_h(psi(z));
        ft_h = @(z) psiT(At_h(z));
        fhandle = @(z,mode) fhandle_mode(z,mode, f_h, ft_h);
        SIGMA = 1e2;% norm(Y(:)-A_h(sceneImgRecPI(:)))*2;
        iPhi_rec = @(z) psi(spg_bpdn(fhandle,z(:),SIGMA, opts));
    case 'lasso'
        opts = [];
        opts.tol = 1e-5;  % Use super strict tolerance
        opts.recordObjective = true; %  Record the objective function so we can plot it
        opts.verbose = true;
        opts.maxIters = 20; 
         
        %         if rectifyCosFactor
        %             F_h = @(z) A_h(CosRectFactor(:).*psi(z));
        %             Ft_h = @(z) psiT(CosRectFactor(:).*At_h(z));
        %             sceneImgRecPI = sceneImgRecPI./CosRectFactor;
        %         else
        %             F_h = @(z) A_h(psi(z));
        %             Ft_h = @(z) psiT(At_h(z));
        %         end
        
        F_h = @(z) A_h(psi(z));
        Ft_h = @(z) psiT(At_h(z));
        % tau = 0.9*sum(abs(psiT(sceneImgRecPI(:))));
        tau = screenPatchSize^2*3;
        
        iPhi_rec = @(z) psi(fasta_lasso(F_h,Ft_h,z,tau, zeros(screenPatchSize^2*numOfChannels,1), opts)); 
        
    case 'nesta-tv'

        % Solver parameters
        muf = 1e-6;
        NESTAmaxiter = 200;
        
        opts = [];
        opts.MaxIntIter = 5;
        opts.Verbose = 0;
        opts.maxiter = NESTAmaxiter;
        opts.TOlVar = 1e-6;
        opts.stoptest = 1;
        opts.ROW = screenPatchSize;
        opts.COL = screenPatchSize; 
        opts.T_frames = 1;

        opts.typemin = 'tv';
        
        % CG for AAtinv
        delta = 0;
        A_function = @(x) A_h(At_h(x));
        cg_tol = 1e-6; cg_maxit = 20;
        CGwrapper(); % (first, zero-out the CGwrapper counters)
        opts.AAtinv = @(b) CGwrapper(A_function,b,cg_tol,cg_maxit);
        % opts.errFcn = @(x) norm( x - I(:)) / norm(I(:));
        % opts.outFcn = @(x) [norm( x - I(:), 'inf' ), norm( x - I(:)) / norm(I(:))];
          
        % NESTA constrained solver
        iPhi_rec = @(z) NESTA(A_h,At_h,z,muf,delta,opts);
        % [xk_NESTA,niter,resid,outputData] = NESTA(A_h,At_h,y,muf,delta,opts);
        % NESTA unconstrained solver
    case 'nesta-L1'        
        
        % Solver parameters
        muf = 1e-6;
        NESTAmaxiter = 200;
        
        opts = [];
        opts.MaxIntIter = 5;
        opts.Verbose = 0;
        opts.maxiter = NESTAmaxiter;
        opts.TOlVar = 1e-6;
        opts.stoptest = 1;
        opts.ROW = screenPatchSize;
        opts.COL = screenPatchSize; 
        opts.T_frames = 1; 
        
        opts.typemin = 'L1';
        opts.U = psiT;
        opts.Ut = psi;
        
        
        % CG for AAtinv
        delta = 0;
        A_function = @(x) A_h(At_h(x));
        cg_tol = 1e-6; cg_maxit = 20;
        CGwrapper(); % (first, zero-out the CGwrapper counters)
        opts.AAtinv = @(b) CGwrapper(A_function,b,cg_tol,cg_maxit);
        opts.errFcn = @(x) norm( x - I(:)) / norm(I(:));
        opts.outFcn = @(x) [norm( x - I(:), 'inf' ), norm( x - I(:)) / norm(I(:))];
         
        % NESTA constrained solver
        iPhi_rec = @(z) NESTA(A_h,At_h,z,muf,delta,opts);
        % [xk_NESTA,niter,resid,outputData] = NESTA(A_h,At_h,y,muf,delta,opts);
        % NESTA unconstrained solver
end

fprintf('\n');