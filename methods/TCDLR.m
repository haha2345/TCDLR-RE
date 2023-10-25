function [Xhat,time,psnr_result,ssim_result,fsim_result,err,iter] = TCDLR(X,p,num,option)
disp(['#######TCDLR########' num2str(num)]);

[n1,n2,n3]=size(X);
opts.mu = 1e-5;
opts.rho = 1.3;
opts.max_iter = 50;
opts.DEBUG = 1;
opts.veck=100*ones(n3);

if isfield(option, 'mu');         opts.mu = option.mu;              end
if isfield(option, 'rho');        opts.rho = option.rho;              end
if isfield(option, 'tol');         opts.tol = option.tol;              end
if isfield(option, 'max_iter');  opts.max_iter = option.max_iter;              end
if isfield(option, 'r')         
    opts.veck=round(n1*0.5)*ones(n3);    
end

opts.rank_min = 25*ones(n3);
opts.lp=option.lp;
rate=option.rate;
opts.rank_max = floor(rate*min(size(X,1),size(X,2))*ones(n3));
maxP = max(X(:));
Omega = find(rand(n1*n2*n3,1)<p);
data=X(Omega);
Xn=zeros(size(X));
Xn(Omega)=data;
tic
[Xhat,err,iter,~] = lrtc_tnn_fastNEW_fixed(Xn,Omega,X,opts);
time=toc;
psnr_result=PSNR(X,Xhat,maxP);
ssim_result=ssim(Xhat,X);
fsim_result=FeatureSIM(X,Xhat);
