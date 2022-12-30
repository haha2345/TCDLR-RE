function [X,Y,coreNway]=EstRankAdjustAB_decreaseNew(X,Y,Xsq,coreNway,rho,rank_adj,mu,n3,rank_min)
%inc=zeros(1,n3);
% first frontal slice
n=1;
ppp=10;
if rank_adj(n) == -1 && rho<1
    max_k=coreNway(n);
    sum_k=coreNway(n);
    sigmas=zeros(max_k,1);
    s = Xsq{n};%svd(Xsq{i},'econ');%
    sigmas(1:length(s))=s;
    [dR,id]=sort(sigmas,'descend');
    drops = dR(1:sum_k-1)./dR(2:sum_k);%lamda_i/lambda_{i+1}
    [dmx,imx] = max(drops);%dmx:hat lambda_{tk} imx:tk
    rel_drp = (sum_k-1)*dmx/(sum(drops)-dmx);%tua^k ����ط������Ĳ�һ�£��Ա��������ĺ����Ӧ���������ı���
    % find��������������
    if rel_drp>ppp
        thold=rho*sum(dR);
        iidx=0;ss=0;
        len=length(dR);
        for i=1:len
            ss=ss+dR(i);
            if(ss>thold)
                iidx=i;
                break;
            end
        end
        if(iidx<coreNway(n))
            if iidx>rank_min(n)
                coreNway(n) =iidx;
            else
                coreNway(n) = rank_min(n);
            end
            [Qx,Rx] = qr(X{n},0);
            [Qy,Ry] = qr(Y{n}',0);
            [U,S,V] = svd(Rx*Ry');
            sigv = diag(S);
            X{n} = Qx*U(:,1:coreNway(n))*spdiags(sigv(1:coreNway(n)),0,coreNway(n),coreNway(n));
            Y{n} = (Qy*V(:,1:coreNway(n)))';
            
        end
        rho=rho*mu;
    end
end


% i=2,...,halfn3
halfn3 = round(n3/2);
for n=2 : halfn3
    if rank_adj(n) == -1 && rho<1
        max_k=coreNway(n);
        sum_k=coreNway(n);
        sigmas=zeros(max_k,1);
        s = Xsq{n};%svd(Xsq{i},'econ');%
        sigmas(1:length(s))=s;%+eps*ones(length(s),1)
        [dR,id]=sort(sigmas,'descend');
        drops = dR(1:sum_k-1)./dR(2:sum_k);%lamda_i/lambda_{i+1}
        [dmx,imx] = max(drops);%dmx:hat lambda_{tk} imx:tk
        rel_drp = (sum_k-1)*dmx/(sum(drops)-dmx);%tua^k ����ط������Ĳ�һ�£��Ա��������ĺ����Ӧ���������ı���
        % find��������������
        
        if rel_drp>ppp
            thold=rho*sum(dR);
            iidx=0;ss=0;
            len=length(dR);
            for i=1:len
                ss=ss+dR(i);
                if(ss>thold)
                    iidx=i;
                    break;
                end
            end
            if(iidx<coreNway(n))
                if iidx>rank_min(n)
                    coreNway(n) =iidx;
                else
                    coreNway(n) = rank_min(n);
                end
                [Qx,Rx] = qr(X{n},0);
                [Qy,Ry] = qr(Y{n}',0);
                [U,S,V] = svd(Rx*Ry');
                sigv = diag(S);
                X{n} = Qx*U(:,1:coreNway(n))*spdiags(sigv(1:coreNway(n)),0,coreNway(n),coreNway(n));
                Y{n} = (Qy*V(:,1:coreNway(n)))';
            end
            rho=rho*mu;
        end
    end
    X{n3+2-n} =conj(X{n});
    Y{n3+2-n} =conj(Y{n});
    coreNway(n3+2-n)=coreNway(n);
end

% if n3 is even
if mod(n3,2) == 0
    n = halfn3+1;
    if rank_adj(n) == -1 && rho<1
        max_k=coreNway(n);
        sum_k=coreNway(n);
        sigmas=zeros(max_k,1);
        s = Xsq{n};%svd(Xsq{i},'econ');%
        sigmas(1:length(s))=s;
        [dR,id]=sort(sigmas,'descend');
        drops = dR(1:sum_k-1)./dR(2:sum_k);%lamda_i/lambda_{i+1}
        [dmx,imx] = max(drops);%dmx:hat lambda_{tk} imx:tk
        rel_drp = (sum_k-1)*dmx/(sum(drops)-dmx);%tua^k ����ط������Ĳ�һ�£��Ա��������ĺ����Ӧ���������ı���
        % find��������������
        if rel_drp>ppp
            thold=rho*sum(dR);
            iidx=0;ss=0;
            len=length(dR);
            for i=1:len
                ss=ss+dR(i);
                if(ss>thold)
                    iidx=i;
                    break;
                end
            end
            if(iidx<coreNway(n))%ԭ����if(iidx>sum(rank_min(n3)))���������⣿Ӧ����if(iidx>sum(rank_min))��
                %idx=floor((id(iidx+1:sum_k)-1)/max_k);%�ҵ����ں�������ֵ������ǰ����Ƭ
                %for n=1:n3
                % num=length(find(idx==n-1));%��ߣ�ޣ�%����㷨��ȱ����ֻ����С���ƣ������ʼֵ�Ƚ�С��ʱ��Ч����Ƚϲ
                %if(num>0)
                if iidx>rank_min(n)
                    coreNway(n) =iidx;
                else
                    coreNway(n) = rank_min(n);
                end
                [Qx,Rx] = qr(X{n},0);
                [Qy,Ry] = qr(Y{n}',0);
                [U,S,V] = svd(Rx*Ry');
                sigv = diag(S);
                X{n} = Qx*U(:,1:coreNway(n))*spdiags(sigv(1:coreNway(n)),0,coreNway(n),coreNway(n));
                Y{n} = (Qy*V(:,1:coreNway(n)))';
                % end
                % end
            end
            rho=rho*mu;
        end
    end
    
end

