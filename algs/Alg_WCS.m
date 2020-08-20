function [pos_WCS] = Alg_WCS...
    (npoints,points,pos_ARAP,neighbors,ConnectivityM,distMatrix,Boxscale)
%ALG_WCS �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
patchCoordMatrix = patchLocalization(npoints, neighbors, distMatrix, ConnectivityM);
G=graph(ConnectivityM);% ����һ��G��ʽ��Edges��Nodes����ͼ��graphΪMATLAB toolbox����
m =sum(sum(ConnectivityM))/2;% �ܱ���
[matching,subgraph]=ordermatch(G, ConnectivityM, m, npoints,points'*Boxscale);
nodeSub = tranGraph(G, ConnectivityM,npoints,subgraph);
[MOCC,ratio] = merge(subgraph,G,ConnectivityM);
[a,b]=size(MOCC);
if a == 1
    pos_WCS = pos_ARAP;
    return
end
MOCCcoord = zeros(2*a,npoints);
for j = 1:a
    id = find(MOCC(j,:));
    this_num = sum(MOCC(j,:));
    Subdist = distMatrix(id,id);
    SubCon = ConnectivityM(id,id);
    Subpoints = points(id,:);
    %coord = runarap( Subpoints,Subdist,SubCon);
    [coord] = Alg_ARAP(SubCon,Subdist,Subpoints);
    MOCCcoord(2*j-1:2*j,id)  = coord';
    na = sum(MOCC(j,:));
end
range = 6;
pos_WCS = linearMerge...
    (points,pos_ARAP,neighbors, patchCoordMatrix,MOCC,MOCCcoord,ratio,range);
end

