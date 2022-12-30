function [A,B,veck]=EstRankAdjustAB_increase(A,B,X,C,veck,n3,rank_max)
l=floor(min(size(X,1),size(X,2))/50);%2;%
C=fft(C,[],3);
X=fft(X,[],3);
[n1,n2,n3]=size(X);
% first frontal slice
for i=1
    if l+veck(i)<rank_max(i)
        % diff=C(:,:,i)-X(:,:,i);
        diff=C(:,:,i)-A{i}*B{i};
        d=diff(1:20,1:20);
        sigma=std(d(:));%norm(diff(1:20,1:20),'fro')/399;%标准差估计
        H=randn(1,n1)*diff;
        [Q,R]=qr(H',0);
        if  abs(R(1))>sigma*(sqrt(size(X,1))+sqrt(size(X,2)))+1;
            
            H=randn(l,n1)*diff;
            [Q,R]=qr(H',0);
  
            A{i}=[A{i} zeros(n1,l)];
            B{i}=[B{i};Q'];
            veck(i)=veck(i)+l;
            
        end
    end
    
end
halfn3 = round(n3/2);
for i=2:halfn3
    if l+veck(i)<rank_max(i)
        %diff=C(:,:,i)-X(:,:,i);
        diff=C(:,:,i)-A{i}*B{i};
        sigma=norm(diff(1:20,1:20),'fro')/399;
        H=randn(1,n1)*diff;
        [Q,R]=qr(H',0);
        if  abs(R(1))>sigma*(sqrt(size(X,1))+sqrt(size(X,2)))+1;
            
            H=randn(l,n1)*diff;
            [Q,R]=qr(H',0);
            A{i}=[A{i} zeros(n1,l)];
            B{i}=[B{i};Q'];
            veck(i)=veck(i)+l;
            A{n3+2-i}=conj(A{i});
            B{n3+2-i}= conj(B{i});
            veck(n3+2-i)=veck(i);
        end
    end
    
end

% if n3 is even
if mod(n3,2) == 0
    i = halfn3+1;
    if l+veck(i)<rank_max(i)
        % diff=C(:,:,i)-X(:,:,i);
        diff=C(:,:,i)-A{i}*B{i};
        sigma=norm(diff(1:20,1:20),'fro')/399;
        H=randn(1,n1)*diff;
        [Q,R]=qr(H',0);
        if  abs(R(1))>sigma*(sqrt(size(X,1))+sqrt(size(X,2)))+1;
            
            H=randn(l,n1)*diff;
            [Q,R]=qr(H',0);
            A{i}=[A{i} zeros(n1,l)];
            B{i}=[B{i};Q'];
            veck(i)=veck(i)+l;
        end
    end
    
end