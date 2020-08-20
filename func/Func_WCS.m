function [pos] = Func_WCS(points,distMatrix,ConnectivityM)
%UNTITLED5 �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
[npoints,~]=size(ConnectivityM);     %nodes number
Boxscale = 100;   %range

nbrNums=sum(ConnectivityM, 2);
maxNbrNum = max(nbrNums);
neighbors=zeros(npoints, maxNbrNum);% neighbors ����Ϊ���ڵ���������Ϊ������ھ���
neighbors(:,1)=nbrNums;% neighbors ��һ�д���ǰ�ڵ��ھ���
for i=1:npoints
    cnt = 1;
    for j=1:npoints
        if ConnectivityM(i,j)==1
            cnt=cnt+1;
            neighbors(i,cnt)=j;
        end
    end
end
% neighbors �������ǵ�ǰ�ڵ��ھӵ�ID
%% ARAP��WCS���õ�
patchCoordMatrix = patchLocalization(npoints, neighbors, distMatrix, ConnectivityM);
[AAAP_Pos]=aaap(npoints, neighbors, patchCoordMatrix);
ARAP_Pos=arap(npoints, AAAP_Pos, neighbors, patchCoordMatrix);
[~,Z,~]=procrustes(points, ARAP_Pos, 'Scaling', false);
ARAP_Pos=Z;
%% WCS
G=graph(ConnectivityM);% ����һ��G��ʽ��Edges��Nodes����ͼ��graphΪMATLAB toolbox����
m =sum(sum(ConnectivityM))/2;% �ܱ���
[~,subgraph]=ordermatch(G, ConnectivityM, m, npoints,points'*Boxscale);% �õ�ȫ�� BRC
%nodeSub = tranGraph(G, ConnectivityM,npoints,subgraph);% ����ûɶ�ã�ע�͵�Ҳ�ܳ����
[MOCC,ratio] = merge(subgraph,G,ConnectivityM);% �ϲ� BRC���õ� MOCC ���Ӧ RR
[a,~]=size(MOCC);
MOCCcoord = zeros(2*a,npoints);% ��ARAP���� MOCC �ĳ�ʼ����
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

