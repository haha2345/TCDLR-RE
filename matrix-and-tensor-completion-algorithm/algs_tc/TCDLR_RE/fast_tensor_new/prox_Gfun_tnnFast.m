function [X,tnn,trank,Bsq] = prox_Gfun_tnnFast(A,B,rho,fun,normpar)
n1=size(A{1},1);
n2=size(B{1},2);
n3=length(A);
X = zeros(n1,n2,n3);
Bsq=cell(1,n3);
tnn = 0;
trank = 0;
% first frontal slice
[U,S,V] = svd(A{1},'econ');
S = diag(S);
S0=fun(S,rho,normpar);
Bsq{1}=S0;
r = length(find(S0>0));
if r>=1
    X(:,:,1) = U(:,1:r)*diag(S0(1:r))*V(:,1:r)'*B{1};
    tnn = tnn+sum(S0);
    trank = max(trank,r);
end
% i=2,...,halfn3
halfn3 = round(n3/2);
for i = 2 : halfn3
    [U,S,V] = svd(A{i},'econ');
    S = diag(S);
    S0=fun(S,rho,normpar);
    Bsq{i}=S0;
    Bsq{n3+2-i}=S0;
    r = length(find(S0>0));
    if r>=1
        X(:,:,i) = U(:,1:r)*diag(S0(1:r))*V(:,1:r)'*B{i};
        tnn = tnn+sum(S0)*2;
        trank = max(trank,r);
    end
    X(:,:,n3+2-i) = conj(X(:,:,i));
end

% if n3 is even

if mod(n3,2) == 0
    i = halfn3+1;
    [U,S,V] = svd(A{i},'econ');
    S = diag(S);
    S0=fun(S,rho,normpar);
    Bsq{i}=S0;
    r = length(find(S0>0));
    if r>=1
        X(:,:,i) = U(:,1:r)*diag(S0(1:r))*V(:,1:r)'*B{i};
        tnn = tnn+sum(S0);
        trank = max(trank,r);
    end
end

tnn = tnn/n3;
X = ifft(X,[],3);
