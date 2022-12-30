function [X,Y,coreNway]=EstRankAdjustAB_decrease(X,Y,Xsq,coreNway,rho,rank_adj,mu,n3,rank_min)

if rank_adj(n3) == -1 && rho<1
        max_k=max(coreNway);
        sum_k=sum(coreNway);
        sigmas=zeros(max_k*n3,1);
        for i=1:n3
            s =svd(Xsq{i},'econ');% Xsq{i};%
            sigmas((i-1)*max_k+1:(i-1)*max_k+length(s))=s;
        end
        [dR,id]=sort(sigmas,'descend');
        drops = dR(1:sum_k-1)./dR(2:sum_k);%lamda_i/lambda_{i+1}
        [dmx,imx] = max(drops);%dmx:hat lambda_{tk} imx:tk
        rel_drp = (sum_k-1)*dmx/(sum(drops)-dmx);%tua^k ����ط������Ĳ�һ�£��Ա��������ĺ����Ӧ���������ı���
        % find��������������
        rel_drp
        if rel_drp>100
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
            if(iidx>sum(rank_min))%ԭ����if(iidx>sum(rank_min(n3)))���������⣿Ӧ����if(iidx>sum(rank_min))��
                idx=floor((id(iidx+1:sum_k)-1)/max_k);%�ҵ����ں�������ֵ������ǰ����Ƭ��
                for n=1:n3
                    num=length(find(idx==n-1));%��ߣ�ޣ�%����㷨��ȱ����ֻ����С���ƣ������ʼֵ�Ƚ�С��ʱ��Ч����Ƚϲ
                    if(num>0)
                        if coreNway(n)-num>rank_min(n)
                            coreNway(n) = coreNway(n)-num;
                        else
                            coreNway(n) = rank_min(n);
                        end
                        [Qx,Rx] = qr(X{n},0);
                        [Qy,Ry] = qr(Y{n}',0);
                        [U,S,V] = svd(Rx*Ry');
                        sigv = diag(S);
                        X{n} = Qx*U(:,1:coreNway(n))*spdiags(sigv(1:coreNway(n)),0,coreNway(n),coreNway(n));
                        Y{n} = (Qy*V(:,1:coreNway(n)))';
                       % C(:,:,n)=X{n}*Y{n};
                    end
                    
                    
                    
                    
                end
            end
            rho=rho*mu;
        end
end
%if rank_adj(n3) == 1



