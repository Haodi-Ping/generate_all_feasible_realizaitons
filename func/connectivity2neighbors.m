function [neighbors] = connectivity2neighbors(ConnectivityM)
%CONNECTIVITY2NEIGHBORS 此处显示有关此函数的摘要
%   此处显示详细说明
nClean = size(ConnectivityM,1);
nbrNums = sum(ConnectivityM,2);
maxNbrNum = max(nbrNums);
neighbors = zeros(nClean, maxNbrNum);% neighbors 行数为：节点数，列数为：最大邻居数
neighbors(:,1)=nbrNums;% neighbors 第一列代表当前节点邻居数
for i=1:nClean
    cnt = 1;
    for j=1:nClean
        if ConnectivityM(i,j)==1
            cnt=cnt+1;
            neighbors(i,cnt)=j;
        end
    end
end

end

