function [X,err,iter,veck] = lrtc_tnn_fastNEW(M,omega,originX,opts)

tol = 1e-6;
max_iter = 500;
rho = 1.1;
mu = 1e-4;
max_mu = 1e14;
DEBUG = 0;
veck=130*ones(1,size(M,3));
rank_min = opts.rank_min;
rank_max = opts.rank_max;
if ~exist('opts', 'var')
    opts = [];
end
if isfield(opts, 'veck');         veck = opts.veck;              end
if isfield(opts, 'tol');         tol = opts.tol;              end
if isfield(opts, 'max_iter');    max_iter = opts.max_iter;    end
if isfield(opts, 'rho');         rho = opts.rho;              end
if isfield(opts, 'mu');          mu = opts.mu;                end
if isfield(opts, 'max_mu');      max_mu = opts.max_mu;        end
if isfield(opts, 'DEBUG');       DEBUG = opts.DEBUG;          end
lp=opts.lp;
dim = size(M);
k = length(dim);
% omegac = setdiff(1:prod(dim),omega);

X = zeros(dim);
X(omega) = M(omega);
E = zeros(dim);
Y = E;
[A,B]=ini_FactorizationTensor(X,veck);

for iter = 1 : max_iter
    Xk = X;
    Ek = E;
    %update A and B
    C=-E+M+Y/mu;
    [A,B,Bsq]=update_FactorizationTensor(C,A,B,veck);
    [A,B,veck]=EstRankAdjustAB_increase(A,B,X,C,veck,dim(3),rank_max);
    
    [X,~,~,Bsq] =prox_Gfun_tnnFast(A,B,1/mu,@Generalized_Soft_Thresholding,lp);%0.8
    
    % update E
    E = M-X+Y/mu;
    E(omega) = 0;
    
    dY = M-X-E;
    chgX = max(abs(Xk(:)-X(:)));
    chgE = max(abs(Ek(:)-E(:)));
    chg = max([chgX chgE max(abs(dY(:)))]);
    
    if DEBUG
        if iter == 1 || mod(iter, 10) == 0
            err = norm(dY(:));
            disp(['iter ' num2str(iter) ', mu=' num2str(mu) ...
                ', err=' num2str(err)]);
        end
    end
    
    if chg < tol
        break;
    end
    Y = Y + mu*dY;
    mu = min(rho*mu,max_mu);
    
    [A,B,veck]=EstRankAdjustAB_decreaseNew(A,B,Bsq,veck,0.999,-1*ones(1,dim(3)),1.01,dim(3),rank_min);
    
end
err = norm(dY(:));

