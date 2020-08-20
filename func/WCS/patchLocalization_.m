function [patchCoordMatrix]=patchLocalization...
    (npoints, neighbors, distMatrix, ConnectivityM)
% 输入：
    % npoints 节点数
    % neighbors 邻居列表阵
    % distMatrix 和 ConnectivityM
% 输出：
    % patchCoordMatrix 每个节点的patch定位结果
    % npoints 个 maxNbrNum+1 × 2 矩阵，即每个 patch 的局部实现结果，
    % 此处为 SMACOF
    
IterNum=50;         % maximum iteration of SMACOF
tol=1e-7;           % tolerance for SMACOF iteration

maxNbrNum=max(neighbors(:,1)); % 邻居数列表
patchCoordMatrix=zeros(npoints, maxNbrNum+1, 2); % 加1是因为要算上自己
%size(patchCoordMatrix)
% 遍历每个节点
for i=1:npoints

    nbrNum=neighbors(i,1); % nbrNum 是当前节点的邻居数
    ind=[neighbors(i,2:nbrNum+1), [i]]; % ind 是当前节点的邻居列表
    pairDistMatrix=zeros(nbrNum+1, nbrNum+1);
    num=nbrNum+1;
    % pairDistMatrix 是当前 patch 的邻接矩阵，对角为0，有边为1，无边为 NaN
    for j=1:num
        for k=j:num
            if j==k
                pairDistMatrix(j,k)=0;
            elseif ConnectivityM(ind(j),ind(k))==1
                pairDistMatrix(j,k)=distMatrix(ind(j),ind(k));
            else
                pairDistMatrix(j,k)=NaN;
            end
        end
    end
    % 猜边，得到完全图   
    for j=1:num
        for k=j:num
            if isnan(pairDistMatrix(j,k)) % 如果 j k 间无边
                maxDist_j=max(distMatrix(ind(j),:)); % maxDist_j 是 j 的最远邻居距离
                maxDist_k=max(distMatrix(ind(k),:)); % maxDist_k 是 k 的最远邻居距离
                lowDist=max(maxDist_j, maxDist_k); % lowDist 是其中更远的那一个
                upDist=Inf; % upDist 是无穷大
                for p=1:npoints
                    if ConnectivityM(ind(j),p)==1 && ConnectivityM(ind(k),p)==1
                        len=distMatrix(ind(j),p)+distMatrix(ind(k),p);
                        upDist=min(upDist, len);
                    end
                end
                pairDistMatrix(j,k)=(lowDist+upDist)*0.5;
            end
        end
    end
    
    % Make matrix symmetry
    for j=1:num
        for k=j:num
            pairDistMatrix(k,j)=pairDistMatrix(j,k);
        end
    end
    
    % MDS
    if(size(pairDistMatrix,1)==1)
        pos=[1,1];
    else
        [Y, e] = mdscale(pairDistMatrix,2);%metricstress
        if(size(Y,2)==1)
            Y(:,2)=zeros()';
        end
        pos=Y(:,1:2);

        % SMACOF
        toIter=1;
        iter=0;
        while toIter
            upos=pos;
            for k=1:num
                nbr=0;
                t=[0,0];
                for p=1:num
                    if ConnectivityM(ind(k),ind(p))==1
                        nbr=nbr+1;
                        z=upos(k,:)-upos(p,:);
                        l=norm(z);
                        if l>0
                            l=1.0/l;
                        end
                        t=t+upos(p,:)+z*l*distMatrix(ind(k),ind(p));
                    end
                end
                if nbr>0
                    t=t*1.0/nbr;
                    pos(k,:)=t;
                end
            end
            iter=iter+1;
            d=norm(pos-upos,'fro');
            if iter>IterNum || d<tol
                toIter=0;
            end
        end
    end
    
    for k=1:nbrNum+1
        patchCoordMatrix(i,ind(k),:) = pos(k,:);
    end
end
%size(patchCoordMatrix)