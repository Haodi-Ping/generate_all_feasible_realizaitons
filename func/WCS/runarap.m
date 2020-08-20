function [ output_args ] = runarap( points,distMatrix,ConnectivityM)

[npoints,opts]=size(points);

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
if isempty(patchCoordMatrix)
    output_args = points;
    return;
end
%% stitching

[AAAP_Pos]=aaap(npoints, neighbors, patchCoordMatrix);
[D,Z,T]=procrustes(points, AAAP_Pos, 'Scaling', false);
registAAAPPos=Z;

ARAP_Pos=arap(npoints, AAAP_Pos, neighbors, patchCoordMatrix);
[D,Z,T]=procrustes(points, ARAP_Pos, 'Scaling', false);
registARAPPos=Z;
output_args=Z;
%output_args=(Z+points)/2;
end

