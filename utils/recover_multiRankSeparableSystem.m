
function [L,R] = recover_multiRankSeparableSystem(sensorImages,screenImages,LInit,RInit)
%
%   Uses the method of multipliers to solve the algorithm
%
%   minimize sum_r ||Lr||^2 + ||Rr||^2
%   subject to Y_i = sum_i sum_r (Lr*X_i*Rr')
%
%   based on the approach described in [1].  The inner iterations descend
%   the extended Lagrangian using minFunc.
%
%   Inputs LInit and RInit are optional.  They give starting points for the
%   algorithm.
%
%
%   References:
%
%   [1] Samuel Burer and Renato D. C. Monteiro.  "A nonlinear programming
%   algorithm for solving semidefinite programs via low-rank
%   factorization." Mathematical Programming (Series B), Vol. 95, 
%   2003. pp 329-357.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The parameters are
%
%   maxOutIter = maximum number of outer iterations
%
%   rmseTol = the root mean square of the errors must drop below
%   rmseTol before termination
%
%   sigmaInit = starting value for the Augmented Lagrangian parameter
%
%   LR1 = the feasibility must drop below LR1 times the previous
%   feasibility before updating the Lagrange multiplier.
%   Otherwise we update sigma.
%
%   LR2 = the amount we update sigma by if feasibility is not improved
%
%   progTol = like LR1, but much smaller.  If feasibility is not
%   improved by progtol for numbaditers iterations, then we decide
%   that further progress cannot be made and quit
%
%   numbaditers = the number of iterations where trivial progress
%   is made improving the feasibility of L and R before algorithm
%   termination.    

global m n T maxrank

pars.maxOutIter = 25;
pars.rmseTol = 1e-8;
pars.sigmaInit = 1e4;
pars.LR1 = 0.25;
pars.LR2 = 10;
pars.progTol = 1e-3;
pars.numbaditers = 6;

% options for minFunc
options = [];
options.display = 'none';
options.maxFunEvals = 50000;
options.Method = 'lbfgs';


% Derived constants:
siglen = length(sensorImages(:)); 

% If starting points are specified, use them. otherwise use the default.
% A better default might speed things up.  Not sure.
if nargin>3 & ~isempty(LInit),
    L = LInit;
else
    L = 1e-2*rand(m*n,maxrank);
    % Z = repmat(1e-2*randn(n1,1),1,maxrank);
end

if nargin>4 & ~isempty(RInit),
    R = RInit;
else
    R = 1e-2*rand(m*n,maxrank);
    % H = repmat(1e-2*randn(n1,1),1,maxrank);
end

% y is the Lagrange multiplier
y = zeros(siglen,1);
% sigma is the penalty parameter
sigma = pars.sigmaInit;

% compute initial infeasibility
dev = op_randomProjSeparableSystem(screenImages,1,L,R) - sensorImages;
vOld = norm(dev(:))^2;


  
v = vOld;
badcnt = 0;
T0 = clock;
  
fprintf('|      |          |          | iter  | tot   |\n');
fprintf('| iter |  rmse    |  sigma   | time  | time  |\n');
for j=1:46, fprintf('-'); end
fprintf('\n');
  
iterCount = 0;

for outIter=1:pars.maxOutIter,
    
    T1 = clock;

    % minimize the Augmented Lagrangian using BFGS
    [x,mval] = minFunc(@subproblem_cost,[L(:);R(:)],options,screenImages,sensorImages,y,sigma);
    L = reshape(x(1:(m*n*maxrank)),m*n,maxrank);
    R = reshape(x((m*n*maxrank) + (1:(m*n*maxrank))),m*n,maxrank);
    
    % compute the equation error 
    dev = op_randomProjSeparableSystem(screenImages,1,L,R) - sensorImages;

    vLast = v;
    v = norm(dev(:))^2; % v is sum of the squares of the errors
    
    % if unable to improve feasibility for several iterations, quit.
    if abs(vLast-v)/vLast<pars.progTol,
      badcnt = badcnt+1;
      if badcnt>pars.numbaditers,
        fprintf('\nunable to make progress. terminating\n');
        break;
      end
    else
      badcnt = 0;
    end
    
    % print diagnostics
    fprintf('| %2d   | %.2e | %.2e |  %3.0f  |  %3.0f  |\n',...
            outIter,sqrt(v/siglen),sigma,etime(clock,T1),etime(clock,T0));

    % if solution is feasible to specified tolerance, we're done
    if sqrt(v/siglen)<pars.rmseTol,
      break;
    % if normed feasibility is greatly reduced, update Lagrange multipliers
    elseif v < pars.LR1*vOld
      y = y - sigma*dev;
      vOld = v;
    % if normed feasibility is not reduced, increase penalty term
    else
      sigma = pars.LR2*sigma;    
    end
    
  end
  % print final diagnostics
  fprintf('elapsed time: %.0f seconds\n',etime(clock,T0));
  
function [mval,g]=subproblem_cost(x,screenImages,sensorImages,y,sigma)
% This helper function computes the value and gradient of the augmented
% Lagrangian     
    
    global m n maxrank T 
    siglen = numel(sensorImages);
      
    L = reshape(x(1:(m*n*maxrank)),m*n,maxrank);
    R = reshape(x((m*n*maxrank) + (1:(m*n*maxrank))),m*n,maxrank);
    
    % compute equation error
    dev = op_randomProjSeparableSystem(screenImages,1,L,R) - sensorImages;
    
    % compute the cost of the extended Lagrangian penalty function
    mval = norm(L,'fro')^2+norm(R,'fro')^2 - 2*real(y'*dev(:)) + sigma*norm(dev(:),'fro')^2; 
    
    % compute the gradient
    yhat = y - sigma*dev(:);

    % Lagrange multiplier and L2 norm terms should be considered jointly 
    %
    g = op_randomProjSeparableSystem(yhat,2,L,R,screenImages);

