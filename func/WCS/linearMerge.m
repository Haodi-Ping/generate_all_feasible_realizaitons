function [ final_Pos ] =  linearMerge (points,initPos,neighbors, patchCoordMatrix,MOCC,MOCCcoord,radio,range)
% points            n*2             Ground truth,转换最终结果用
% initPos           n*2             初始化目标函数的P，在这里是直接用 ARAP 给的
% neighbors         n*(maxNeiNum+1)	第一列是邻居数，后面是邻居ID，空余的填0
% patchCoordMatrix  n*n*2           2个n*n的矩阵 ?
% MOCC              subNum*n        存储 MOCC 子图
% MOCCcoord         (2subNum)*n     MOCC的ARAP实现结果
% radio             1*subNum        每个子图的冗余比
% range             常数            算weight用

rmin = min(radio);
rmax = max(radio);
[moccnum,npoints]=size(MOCC);
if rmin == rmax
    weight = ones(1,moccnum)*range;
else
    weight = 1+(radio-rmin)/(rmax - rmin)*range;
end
%weight = radio;
IterNum=300;
tol=1e-9;

I=0;
J=0;
S=0;
row=0;
cnt=0;
nnz=0;

% Constraint
kc=1;
for kc=1:npoints
    if neighbors(kc,1)>2
        break;
    end
end

% Prefactor coefficient matrix
for i=1:npoints
    ind=neighbors(i,:);
    nbrNum=ind(1);
    s=cat(2,[i],ind(2:nbrNum+1));
    s=sort(s);
    num=nbrNum+1;
    
    d=0;
    if i==kc
        d=1.0/2;
    end
    
    row=row+1;
    for j=1:num
        cnt=cnt+1;
        I(cnt)=row;
        J(cnt)=s(j);
        if s(j)==i
            S(cnt)=nbrNum+d;
        else
            S(cnt)=-1;
        end
    end        
    nnz=nnz+num;
end

m=row;
n=npoints;
A=sparse(I',J',S',m,n,nnz);
[L,U]=lu(A);    



%ARAP迭代过程
pos=initPos;
toIter=1;
iter=0;
while toIter
    initPos=pos;
    R=zeros(npoints, 4);
    for i=1:npoints
        ind=neighbors(i,:);
        X=patchCoordMatrix(i,:,1);
        Y=patchCoordMatrix(i,:,2);
        nbrNum=ind(1);
        P=zeros(nbrNum,2);
        Q=zeros(nbrNum,2);
        for j=2:nbrNum+1
            Q(j-1,:)=[X(ind(j))-X(i),Y(ind(j))-Y(i)];
            P(j-1,:)=pos(ind(j),:)-pos(i,:);
        end
        [E,S,F]=svd(Q'*P);
        H=E*F';
        R(i,:)=reshape(H,1,4);
    end
    
    r=zeros(npoints,2);
    ind=neighbors(kc,:);
    Xi=patchCoordMatrix(kc,:,1);
    Yi=patchCoordMatrix(kc,:,2);
    r(kc,:)=[Xi(kc),Yi(kc)];
    nbrNum=ind(1);
    for j=2:nbrNum+1
        r(ind(j),:)=[Xi(ind(j)),Yi(ind(j))]*0.5;
    end
    
    for i=1:npoints
        ind=neighbors(i,:);
        Xi=patchCoordMatrix(i,:,1);
        Yi=patchCoordMatrix(i,:,2);
        nbrNum=ind(1);
        Ri=reshape(R(i,:),2,2);
        
        t=[0,0];
        
        for j=2:nbrNum+1
            Xj=patchCoordMatrix(ind(j),:,1);
            Yj=patchCoordMatrix(ind(j),:,2);
            Rj=reshape(R(ind(j),:), 2, 2);
            totalweight = 0;
            pp = [0,0];
            for k = 1:moccnum
                if(MOCC(k,i)*MOCC(k,ind(j)))
                    totalweight = totalweight + weight(k);
                    pp = pp + weight(k)*[MOCCcoord(2*k-1,i)-MOCCcoord(2*k-1,ind(j)),MOCCcoord(2*k,i)-MOCCcoord(2*k,ind(j))];
                end
            end
            t=t+([Xi(i)-Xi(ind(j)), Yi(i)-Yi(ind(j))]*Ri+[Xj(i)-Xj(ind(j)), Yj(i)-Yj(ind(j))]*Rj+pp)/(totalweight+2);
        end
        %r(i,:)=r(i,:)+t;
        r(i,:)=t;
    end
    
    pos=U\(L\r);
    if(sum(sum(isnan(pos)))>0)
        pos
    end
    [D,Z,T]=procrustes(points, pos, 'Scaling', false);
    pos=Z;
    
    iter=iter+1;
    d=norm(initPos-pos, 'fro');
    
    if iter>IterNum || d<tol
        toIter=0;
    end
end

final_Pos=pos;





end

