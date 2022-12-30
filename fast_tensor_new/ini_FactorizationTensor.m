function [A,B]=ini_FactorizationTensor(Y,veck)
[n1,n2,n3]=size(Y);
A=cell(n3);
B=cell(n3);
Y= fft(Y,[],3);
%  k=max(veck);
%  A=zeros(n1,k,n3);
%  B=zeros(k,n2,n3);
%  % first frontal slice
%   [Uk,Sigmak,Vk]=svds(Y(:,:,1),k);
%      A(:,:,1)=Uk*Sigmak;
%      B(:,:,1)=Vk'; 
%  halfn3 = round(n3/2);
%  for i = 2 : halfn3
%      [Uk,Sigmak,Vk]=svds(Y(:,:,i),k);
%      A(:,:,i)=Uk*Sigmak;
%      B(:,:,i)=Vk'; 
%      A(:,:,n3+2-i)=conj(A(:,:,i));
%      B(:,:,n3+2-i)=conj(B(:,:,i));
%  end
%  % if n3 is even
% if mod(n3,2) == 0
%     i = halfn3+1;
%     [Uk,Sigmak,Vk]=svds(Y(:,:,i),k);
%      A(:,:,i)=Uk*Sigmak;
%      B(:,:,i)=Vk'; 
% end
% %   for i = 1 :n3
% %      [Uk,Sigmak,Vk]=svds(Y(:,:,i),k);
% %      A(:,:,i)=Uk*Sigmak;
% %      B(:,:,i)=Vk'; 
% %   end
%  A=ifft(A,[],3);
%  B=ifft(B,[],3);
  %  for i = 1 :n3
%      [Uk,Sigmak,Vk]=svds(Y(:,:,i),veck(i));
%      A{i}=Uk*Sigmak;
%      B{i}=Vk';
%  end
% first frontal slice
 [Uk,Sigmak,Vk]=svds(Y(:,:,1),veck(1));
 A{1}=Uk*Sigmak;
 B{1}=Vk';
 halfn3 = round(n3/2);
 for i = 2 : halfn3
     [Uk,Sigmak,Vk]=svds(Y(:,:,i),veck(i));
     A{i}=Uk*Sigmak;
     B{i}=Vk';
     A{n3+2-i}=conj(Uk*Sigmak);
     B{n3+2-i}=conj(Vk');
 end
 % if n3 is even
if mod(n3,2) == 0
    i = halfn3+1;
    [Uk,Sigmak,Vk]=svds(Y(:,:,i),veck(i));
    A{i}=Uk*Sigmak;
    B{i}=Vk';
end
