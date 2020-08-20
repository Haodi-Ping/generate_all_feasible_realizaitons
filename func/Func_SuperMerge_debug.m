function [ final_Pos ] =  Func_SuperMerge_debug...
    (points,initPos,neighbors, patchCoordMatrix,MOCC,MOCCcoord,radio)
% points            n*2             Ground truth,转换最终结果用
% initPos           n*2             初始化目标函数的P，在这里是直接用 ARAP 给的
% neighbors         n*(maxNeiNum+1)	第一列是邻居数，后面是邻居ID，空余的填0
% patchCoordMatrix  n*n*2           2个n*n的矩阵 ?
% MOCC              subNum*n        存储 MOCC 子图
% MOCCcoord         (2subNum)*n     MOCC的ARAP实现结果
% radio             1*subNum        每个子图的冗余比
% range             常数            算weight用
show_convergence = 0;
if show_convergence == 1
    mat_path = '.\mats\iteration\';
    delete ./mats/iteration/*.mat
end

range = 6;
[moccnum,npoints]=size(MOCC);
rmin = min(radio);
rmax = max(radio);
if rmin==rmax
    weight = radio;
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
if show_convergence == 1
    this_pos = strcat('pos_',num2str(0));
    eval([this_pos '=pos;']);
    save([mat_path,'pos_',num2str(0),'.mat'],this_pos)
    times = [];
end
toIter=1;
iter=0;

while toIter
    tic
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
        r(i,:) = t;
    end
    
    pos=U\(L\r);
    if(sum(sum(isnan(pos)))>0)
        pos
    end
    [D,Z,T]=procrustes(points, pos, 'Scaling', false);
    pos=Z;
    
    iter=iter+1;
    d=norm(initPos-pos, 'fro');
    
    %     if iter>IterNum || d<tol
    %         toIter=0;
    %     end
    if show_convergence == 1
        
        if iter>IterNum
            toIter=0;
            disp(['max iteration number reached --- ', num2str(IterNum)]);
        end
        if d<tol
            toIter=0;
            disp(['diff between two iterations is within the threshould in --- ',num2str(iter)]);
        end
        this_time = toc;
        times = [times this_time];
        this_pos = strcat('pos_',num2str(iter));
        eval([this_pos '=pos;']);
        save([mat_path,'pos_',num2str(iter),'.mat'],this_pos)
    end
end
if show_convergence == 1
    
    save([mat_path,'times.mat'],'times')
    final_Pos=pos;
    save([mat_path,'points.mat'],'points')
    
end
end