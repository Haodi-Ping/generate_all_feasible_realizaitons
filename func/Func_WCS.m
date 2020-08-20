function [pos] = Func_WCS(points,distMatrix,ConnectivityM)
%UNTITLED5 此处显示有关此函数的摘要
%   此处显示详细说明
[npoints,~]=size(ConnectivityM);     %nodes number
Boxscale = 100;   %range

nbrNums=sum(ConnectivityM, 2);
maxNbrNum = max(nbrNums);
neighbors=zeros(npoints, maxNbrNum);% neighbors 行数为：节点数，列数为：最大邻居数
neighbors(:,1)=nbrNums;% neighbors 第一列代表当前节点邻居数
for i=1:npoints
    cnt = 1;
    for j=1:npoints
        if ConnectivityM(i,j)==1
            cnt=cnt+1;
            neighbors(i,cnt)=j;
        end
    end
end
% neighbors 其它列是当前节点邻居的ID
%% ARAP，WCS会用到
patchCoordMatrix = patchLocalization(npoints, neighbors, distMatrix, ConnectivityM);
[AAAP_Pos]=aaap(npoints, neighbors, patchCoordMatrix);
ARAP_Pos=arap(npoints, AAAP_Pos, neighbors, patchCoordMatrix);
[~,Z,~]=procrustes(points, ARAP_Pos, 'Scaling', false);
ARAP_Pos=Z;
%% WCS
G=graph(ConnectivityM);% 返回一个G格式【Edges，Nodes】的图，graph为MATLAB toolbox函数
m =sum(sum(ConnectivityM))/2;% 总边数
[~,subgraph]=ordermatch(G, ConnectivityM, m, npoints,points'*Boxscale);% 得到全部 BRC
%nodeSub = tranGraph(G, ConnectivityM,npoints,subgraph);% 好像没啥用，注释掉也能出结果
[MOCC,ratio] = merge(subgraph,G,ConnectivityM);% 合并 BRC，得到 MOCC 与对应 RR
[a,~]=size(MOCC);
MOCCcoord = zeros(2*a,npoints);% 用ARAP给出 MOCC 的初始坐标
for j = 1:a
    id = find(MOCC(j,:));
    Subdist = distMatrix(id,id);
    SubCon = ConnectivityM(id,id);
    Subpoints = points(id,:);
    coord = runarap( Subpoints,Subdist,SubCon);
    MOCCcoord(2*j-1:2*j,id)  = coord';
    na = sum(MOCC(j,:));
end
range = 6; % 
% ARAP_Pos
% neighbors
% patchCoordMatrix
% MOCC
% MOCCcoord
% ratio
pos = linearMerge (points,ARAP_Pos,neighbors, patchCoordMatrix,MOCC,MOCCcoord,ratio,range);


end

