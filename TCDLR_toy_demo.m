addpath(genpath(cd))
clear all;
clc;
% TCDLR in toy data with fixed rank
method={'TCDLR'};
short_name={'TCDLR'};

lp=0.8;
rate=0.5;
rho=1.5;
option.lp=lp;
option.rate=rate;
option.mu = 1e-5;
option.rho = rho;
option.tol=1e-9;
option.max_iter = 100;
option.k_max=1.5;

for n1=[1000,2000,3000,4000]
    r=floor(n1*0.1);
    X1=randn(n1,r,3);
    X2=randn(r,n1,3);
    I=tprod(X1,X2);
    option.r=r;
    p=0.3; %sampling rate
    j=1; % methods
    for i=1:5 % times
        X=double(I);
        funcc=str2func(method{j});
        [Xhat,time] = funcc(X,p,i,option);
        TIME_(i)=time;
        RSE(i)=norm(Xhat(:)-X(:),'fro')/norm(X(:),'fro');
    end
    p1=round(p*10);
    time_list(p1,r)=mean(TIME_);
    rse_list(p1,r)=mean(RSE);
end


