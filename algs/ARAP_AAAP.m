function [pos_AAAP,pos_ARAP] = ARAP_AAAP(ConnectivityM,distMatrix,points)
%UNTITLED5 此处显示有关此函数的摘要
%   此处显示详细说明
[npoints,~]=size(points);

% if npoints ==3 && (sum(sum(ConnectivityM)) == 6)
%     pos_ARAP = points;
%     return
% end

%% patch compute
nbrNums=sum(ConnectivityM, 2);
maxNbrNum=max(nbrNums);
neighbors=zeros(npoints, maxNbrNum);
neighbors(:,1)=nbrNums;
for i=1:npoints
    cnt = 1;
    for j=1:npoints
        if ConnectivityM(i,j)==1
            cnt=cnt+1;
            neighbors(i,cnt)=j;
        end
    end
end
patchCoordMatrix=patchLocalization(npoints, neighbors, distMatrix, ConnectivityM);
% if isempty(patchCoordMatrix)
%     pos_ARAP = points;
%     return;
% end
%% stitching

[AAAP_Pos]=aaap(npoints, neighbors, patchCoordMatrix);
%[D,Z,T]=procrustes(points, AAAP_Pos, 'Scaling', false);
%registAAAPPos=Z;

ARAP_Pos=arap(npoints, AAAP_Pos, neighbors, patchCoordMatrix);
[D,Z,T]=procrustes(points, ARAP_Pos, 'Scaling', false);
%registARAPPos=Z;
pos_ARAP=Z;
[D,pos_AAAP,T]=procrustes(points, AAAP_Pos, 'Scaling', false);

%output_args=(Z+points)/2;
end

