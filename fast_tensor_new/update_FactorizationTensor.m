function [newA,newB,Bsq]=update_FactorizationTensor(Y,A,B,veck)
[n1,n2,n3]=size(Y);
newA=cell(n3);
newB=cell(n3);
Bsq=cell(n3);
Y= fft(Y,[],3);

% first frontal slice
newA{1}=Y(:,:,1)*B{1}'*pinv(B{1}*B{1}');
Bsq{1}=newA{1}'*newA{1};
newB{1}=pinv(Bsq{1})*newA{1}'*Y(:,:,1);
%提取
[Q,R]=qr(newB{1}',0);
newA{1}=newA{1}*R';
newB{1}=Q';
halfn3 = round(n3/2);
for i = 2 : halfn3
    newA{i}=Y(:,:,i)*B{i}'*pinv(B{i}*B{i}');
    Bsq{i}=newA{i}'*newA{i};
    newB{i}=pinv(Bsq{i})*newA{i}'*Y(:,:,i);
    [Q,R]=qr(newB{i}',0);
    newA{i}=newA{i}*R';
    newB{i}=Q';
    
    newA{n3+2-i}=conj(newA{i});
    Bsq{n3+2-i}=conj(Bsq{i});
    newB{n3+2-i}=conj(newB{i});
end
% if n3 is even
if mod(n3,2) == 0
    i = halfn3+1;
    newA{i}=Y(:,:,i)*B{i}'*pinv(B{i}*B{i}');
    Bsq{i}=newA{i}'*newA{i};
    newB{i}=pinv(Bsq{i})*newA{i}'*Y(:,:,i);
    [Q,R]=qr(newB{i}',0);
    newA{i}=newA{i}*R';
    newB{i}=Q';
end
