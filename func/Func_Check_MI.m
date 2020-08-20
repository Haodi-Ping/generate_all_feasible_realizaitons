function [distMatrix,ConnectivityM,flag_MI] = Func_Check_MI...
    (cut_set,subgraph,radius,points,distMatrix,ConnectivityM)
%Func_CheckMI 输入pos_1,pos_2检查这两个模块中是否违背MI
%   此处显示详细说明
disp('Func_Check_MI Begin');
%cut_set
flag_MI = 0; % 为0表示不违背，此2种实现可行
%% 参数初始化
graph1 = subgraph(1,:);
graph2 = subgraph(2,:);
ID1 = find(graph1>0);
subPP1 = points(ID1,:);
subDis1 = distMatrix(ID1,ID1);
subCon1 = ConnectivityM(ID1,ID1);
ID2 = find(graph2>0);
subPP2 = points(ID2,:);
subDis2 = distMatrix(ID2,ID2);
subCon2 = ConnectivityM(ID2,ID2);
%try
%% 局部实现
disp('开始局部实现');
try
    pos_1 = Func_WCS(subPP1,subDis1, subCon1);
catch
    pos_1 = Alg_ARAP(subPP1,subDis1, subCon1);
end
try
    pos_2 = Func_WCS(subPP2,subDis2, subCon2);
catch
    pos_2 = Alg_ARAP(subPP2,subDis2, subCon2);
end
[~,Z1,~]=procrustes(subPP1, pos_1, 'Scaling', false);
[~,Z2,~]=procrustes(subPP2, pos_2, 'Scaling', false);
pos_1 = Z1;
pos_2 = Z2;
%----- plot 出局部模块的实现结果
% plotpositionsaa(subPP1',pos_1',subCon1,100,'G1的实现结果');
% hold on;
% plotsubgraph(points'*100,ConnectivityM,subgraph);
%
% plotpositionsaa(subPP2',pos_2',subCon2,100,'G2的实现结果');
% hold on;
% plotsubgraph(points'*100,ConnectivityM,subgraph);


% Y1=pdist(pos_1);
% dist1=squareform(Y1);
% Y2=pdist(pos_2);
% dist2=squareform(Y2);
% figure;
% plotgraph([],pos_1'*100,subCon1);%绘00制拓扑图
% %hold on
% figure;
% plotgraph([],pos_2'*100,subCon2);%绘00制拓扑图

cut1_ID_in_graph1 = find (ID1 == cut_set(1));
i_l = pos_1(cut1_ID_in_graph1,:);
cut2_ID_in_graph1 = find (ID1 == cut_set(2));
j_l = pos_1(cut2_ID_in_graph1,:);
cut1_ID_in_graph2 = find (ID2 == cut_set(1));
i_r = pos_2(cut1_ID_in_graph2,:);
cut2_ID_in_graph2 = find (ID2 == cut_set(2));
j_r = pos_2(cut2_ID_in_graph2,:);
%% 如果刚性变换时 {vi,vj} 对不齐，这种情况，将点多的实现结果赋给点少的，然后重新实现-----对不齐不用管
%
% if(ConnectivityM(cut_set(1),cut_set(2))==0)
%     disp('局部实现调整...');
%     num1 = sum(subgraph(1,:));
%     num2 = sum(subgraph(2,:));
%     if num1 > num2
%         disp('---左边节点多，把 cut_set 距离赋给右边，让右边重新算');
%         % 左边节点多，把 cut_set 距离赋给右边，让右边重新算
%         e_ij_1 = (i_l(1)-j_l(1))^2 + (i_l(2)-j_l(2))^2;
%         e_ij_1 = sqrt(e_ij_1);
%         subDis2(cut1_ID_in_graph2,cut2_ID_in_graph2) = e_ij_1;
%         subDis2(cut2_ID_in_graph2,cut1_ID_in_graph2) = e_ij_1;
%         subCon2(cut1_ID_in_graph2,cut2_ID_in_graph2) = 1;
%         subCon2(cut2_ID_in_graph2,cut1_ID_in_graph2) = 1;
%         pos_2 = Alg_ARAP(subPP2,subDis2, subCon2);
%         i_r = pos_2(cut1_ID_in_graph2,:);
%         j_r = pos_2(cut2_ID_in_graph2,:);
%     else
%         disp('---右边节点多，把 cut_set 距离赋给左边，让左边重新算');
%         % 右边节点多，把 cut_set 距离赋给左边，让左边重新算
%         e_ij_2 = (i_r(1)-j_r(1))^2 + (i_r(2)-j_r(2))^2;
%         e_ij_2 = sqrt(e_ij_2);
%         subDis1(cut1_ID_in_graph1,cut2_ID_in_graph1) = e_ij_2;
%         subDis1(cut2_ID_in_graph1,cut1_ID_in_graph1) = e_ij_2;
%         subCon1(cut1_ID_in_graph1,cut2_ID_in_graph1) = 1;
%         subCon1(cut2_ID_in_graph1,cut1_ID_in_graph1) = 1;
%         pos_1 = Alg_ARAP(subPP1,subDis1, subCon1);
%         i_l = pos_1(cut1_ID_in_graph1,:);
%         j_l = pos_1(cut2_ID_in_graph1,:);
%     end
% end
%% 刚性变换
[pos2_2_pos1] = Func_RigidTransform(i_l,j_l,i_r,j_r,pos_2);
% [~,Z1,~]=procrustes(subPP1, pos_1, 'Scaling', false);
% [~,Z2,~]=procrustes(subPP2, pos2_2_pos1, 'Scaling', false);
%figure;
%plotgraph([],pos_1'*100,subCon1);%绘00制拓扑图
%hold on
%plotgraph([],pos2_2_pos1'*100,subCon2);%绘00制拓扑图



%i_l
%j_l
%pos2_2_pos1(cut1_ID_in_graph2,:)
%pos2_2_pos1(cut2_ID_in_graph2,:)

%% MI check
% 如果P1没问题，那么检查P2，若P2有问题，则根据P1加边；若P2没问题，无法依据MI进行任何加边。
% 如果P1有问题，那么检查P2，若P2没问题，则根据P2加边；若P2也有问题，说明此拆分无解，应当舍弃此次拆分，将2个子图当成1个图实现；
flag_P1 = 0;
flag_P2 = 0;
% 检查P1
pos_11 = pos_1;
cut_in_g1 = [cut1_ID_in_graph1 cut2_ID_in_graph1];
pos_11(cut_in_g1,:) = [];
pos_22 = pos2_2_pos1;
cut_in_g2 = [cut1_ID_in_graph2 cut2_ID_in_graph2];
pos_22(cut_in_g2,:) = [];
[k1,~] = size (pos_11);
[k2,~] = size (pos_22);

% plot P1--------------------------------------
subPP = [subPP1;subPP2];
subPos = [pos_1;pos2_2_pos1];
[a,~] = size(subPos);
Y=pdist(subPos);
trueDis=squareform(Y);
% 生成邻接矩阵
subCon=zeros(a, a);
for i=1:a
    for j=1:a
        if trueDis(i,j)>radius
            subCon(i,j)=0;
        elseif i==j
            subCon(i,j)=0;
        else
            subCon(i,j)=1;
        end
    end
end
% plotpositionsaa(subPP',subPos',subCon,100,'P1');
% hold on;
% plotsubgraph(points'*100,ConnectivityM,subgraph);
%--------------------------------------

for i = 1:k1
    for j = 1:k2
        distance_ij = (pos_22(j,1)-pos_11(i,1))^2 + (pos_22(j,2)-pos_11(i,2))^2;
        distance_ij = sqrt(distance_ij);
        %[i,j,distance_ij]
        if (distance_ij < radius)
            flag_P1 = 1;
            disp('P1有问题,MI成立');
            break;
        end
    end
    if (flag_P1 == 1)
        break;
    end
end
% 检查P2
pos2_flip = Func_Flip(pos2_2_pos1,cut1_ID_in_graph2,cut2_ID_in_graph2);
%[~,Z4,~]=procrustes(subPP1, pos_1, 'Scaling', false);
%[~,Z3,~]=procrustes(subPP2, pos2_flip, 'Scaling', false);
%figure;
% plotgraph([],pos_1'*100,subCon1);%绘00制拓扑图
% hold on
% plotgraph([],pos2_flip'*100,subCon2);%绘00制拓扑图

% plot P2--------------------------------------
% subPP = [subPP1;subPP2];
subPos_P2 = [pos_1;pos2_flip];
% [a,~] = size(subPos);
Y=pdist(subPos_P2);
trueDis=squareform(Y);
% 生成邻接矩阵
subCon=zeros(a, a);
for i=1:a
    for j=1:a
        if trueDis(i,j)>radius
            subCon(i,j)=0;
        elseif i==j
            subCon(i,j)=0;
        else
            subCon(i,j)=1;
        end
    end
end

% pause(0.5)
% plotpositionsaa(subPP',subPos_P2',subCon,100,'P2');
% hold on;
% plotsubgraph(points'*100,ConnectivityM,subgraph);
%--------------------------------------


pos22_flip = pos2_flip;
pos22_flip(cut_in_g2,:) = [];
for i = 1:k1
    for j = 1:k2
        distance_ij = (pos22_flip(j,1)-pos_11(i,1))^2 + (pos22_flip(j,2)-pos_11(i,2))^2;
        distance_ij = sqrt(distance_ij);
        %[i,j,distance_ij]
        if (distance_ij < radius)
            flag_P2 = 1;%说明P2有问题
            disp('P2有问题，MI成立');
            break;
        end
    end
    if (flag_P2 == 1)
        break;
    end
end

if (flag_P1 == 0) % P1没问题
    %disp('P1没问题,开始检查P2');
    if (flag_P2 == 0)% P2也没问题，MI没有加边
        disp('CASE 1: P1，P2都没问题，二者各自保留且MI不加边。');
    else % P2有问题,根据P1加边
        flag_MI = 1;
        disp('CASE 2: P1没问题，P2有问题，二合一并根据P1加边。');
        % 暂定分别从g1和g2中，各取cut1和cut2的1个邻居，进行加边
        % 加cut1的边
        neiIDs_cut1_in_g1 = find(subCon1(cut1_ID_in_graph1,:)>0);% cut1在g1的全部邻居
        neiID_cut1_in_g1 = neiIDs_cut1_in_g1(1);% 取第一个
        neiID_cut1_in_g1_global = ID1(neiID_cut1_in_g1);% 换算到global ID
        neiIDs_cut1_in_g2 = find(subCon2(cut1_ID_in_graph2,:)>0);% cut1在g2的全部邻居
        neiID_cut1_in_g2 = neiIDs_cut1_in_g2(1);% 取第一个
        neiID_cut1_in_g2_global = ID2(neiID_cut1_in_g2);% 换算到global ID
        dis1 = (pos_1(neiID_cut1_in_g1,1) - pos2_2_pos1(neiID_cut1_in_g2,1))^2+(pos_1(neiID_cut1_in_g1,2) - pos2_2_pos1(neiID_cut1_in_g2,2))^2;
        dis1 = sqrt(dis1);% 根据局部实现算距离
        distMatrix(neiID_cut1_in_g1_global,neiID_cut1_in_g2_global) = dis1;
        distMatrix(neiID_cut1_in_g2_global,neiID_cut1_in_g1_global) = dis1;% 更新 global distMatrix
        ConnectivityM (neiID_cut1_in_g1_global,neiID_cut1_in_g2_global) = 1;
        ConnectivityM (neiID_cut1_in_g2_global,neiID_cut1_in_g1_global) = 1;% 更新 global ConnectivityM
        % 加cut2的边
        neiIDs_cut2_in_g1 = find(subCon1(cut2_ID_in_graph1,:)>0);% cut2在g1的全部邻居
        neiID_cut2_in_g1 = neiIDs_cut2_in_g1(1);% 取第一个
        neiID_cut2_in_g1_global = ID1(neiID_cut2_in_g1);% 换算到global ID
        neiIDs_cut2_in_g2 = find(subCon2(cut2_ID_in_graph2,:)>0);% cut1在g2的全部邻居
        neiID_cut2_in_g2 = neiIDs_cut2_in_g2(1);% 取第一个
        neiID_cut2_in_g2_global = ID2(neiID_cut2_in_g2);% 换算到global ID
        dis2 = (pos_1(neiID_cut2_in_g1,1) - pos2_2_pos1(neiID_cut2_in_g2,1))^2+(pos_1(neiID_cut2_in_g1,2) - pos2_2_pos1(neiID_cut2_in_g2,2))^2;
        dis2 = sqrt(dis2);% 根据局部实现算距离
        distMatrix(neiID_cut2_in_g1_global,neiID_cut2_in_g2_global) = dis2;
        distMatrix(neiID_cut2_in_g2_global,neiID_cut2_in_g1_global) = dis2;% 更新 global distMatrix
        ConnectivityM (neiID_cut2_in_g1_global,neiID_cut2_in_g2_global) = 1;
        ConnectivityM (neiID_cut2_in_g2_global,neiID_cut2_in_g1_global) = 1;% 更新 global ConnectivityM
        
        %-------输出加边效果
        disp('加边：【起点ID 终点ID GroundTruth 实现距离】');
        d1 = (points(neiID_cut1_in_g1_global,1) - points(neiID_cut1_in_g2_global,1))^2 + (points(neiID_cut1_in_g1_global,2) - points(neiID_cut1_in_g2_global,2))^2;
        d1 = sqrt(d1);
        d2 = (points(neiID_cut2_in_g1_global,1) - points(neiID_cut2_in_g2_global,1))^2 + (points(neiID_cut2_in_g1_global,2) - points(neiID_cut2_in_g2_global,2))^2;
        d2 = sqrt(d2);
        [neiID_cut1_in_g1_global neiID_cut1_in_g2_global d1 dis1]
        [neiID_cut2_in_g1_global neiID_cut2_in_g2_global d2 dis2]
        
    end
else % P1有问题
    flag_MI = 1;
    if (flag_P2 == 0)% P2没问题，根据P2加边（pos_2_flip）
        disp('CASE 3: P1有问题，P2没问题，二合一并根据P2（P1_flip）加边。');
        % 暂定分别从g1和g2中，各取cut1和cut2的1个邻居，进行加边
        % 加cut1的边
        neiIDs_cut1_in_g1 = find(subCon1(cut1_ID_in_graph1,:)>0);% cut1在g1的全部邻居
        neiID_cut1_in_g1 = neiIDs_cut1_in_g1(1);% 取第一个
        neiID_cut1_in_g1_global = ID1(neiID_cut1_in_g1);% 换算到global ID
        neiIDs_cut1_in_g2 = find(subCon2(cut1_ID_in_graph2,:)>0);% cut1在g2的全部邻居
        neiID_cut1_in_g2 = neiIDs_cut1_in_g2(1);% 取第一个
        neiID_cut1_in_g2_global = ID2(neiID_cut1_in_g2);% 换算到global ID
        dis1 = (pos_1(neiID_cut1_in_g1,1) - pos2_flip(neiID_cut1_in_g2,1))^2+(pos_1(neiID_cut1_in_g1,2) - pos2_flip(neiID_cut1_in_g2,2))^2;
        dis1 = sqrt(dis1);% 根据局部实现算距离
        distMatrix(neiID_cut1_in_g1_global,neiID_cut1_in_g2_global) = dis1;
        distMatrix(neiID_cut1_in_g2_global,neiID_cut1_in_g1_global) = dis1;% 更新 global distMatrix
        ConnectivityM (neiID_cut1_in_g1_global,neiID_cut1_in_g2_global) = 1;
        ConnectivityM (neiID_cut1_in_g2_global,neiID_cut1_in_g1_global) = 1;% 更新 global ConnectivityM
        % 加cut2的边
        neiIDs_cut2_in_g1 = find(subCon1(cut2_ID_in_graph1,:)>0);% cut2在g1的全部邻居
        neiID_cut2_in_g1 = neiIDs_cut2_in_g1(1);% 取第一个
        neiID_cut2_in_g1_global = ID1(neiID_cut2_in_g1);% 换算到global ID
        neiIDs_cut2_in_g2 = find(subCon2(cut2_ID_in_graph2,:)>0);% cut1在g2的全部邻居
        neiID_cut2_in_g2 = neiIDs_cut2_in_g2(1);% 取第一个
        neiID_cut2_in_g2_global = ID2(neiID_cut2_in_g2);% 换算到global ID
        dis2 = (pos_1(neiID_cut2_in_g1,1) - pos2_flip(neiID_cut2_in_g2,1))^2+(pos_1(neiID_cut2_in_g1,2) - pos2_flip(neiID_cut2_in_g2,2))^2;
        dis2 = sqrt(dis2);% 根据局部实现算距离
        distMatrix(neiID_cut2_in_g1_global,neiID_cut2_in_g2_global) = dis2;
        distMatrix(neiID_cut2_in_g2_global,neiID_cut2_in_g1_global) = dis2;% 更新 global distMatrix
        ConnectivityM (neiID_cut2_in_g1_global,neiID_cut2_in_g2_global) = 1;
        ConnectivityM (neiID_cut2_in_g2_global,neiID_cut2_in_g1_global) = 1;% 更新 global ConnectivityM
        
        %-------输出加边效果
        disp('加边：【起点ID 终点ID GroundTruth 实现距离】');
        d1 = (points(neiID_cut1_in_g1_global,1) - points(neiID_cut1_in_g2_global,1))^2 + (points(neiID_cut1_in_g1_global,2) - points(neiID_cut1_in_g2_global,2))^2;
        d1 = sqrt(d1);
        d2 = (points(neiID_cut2_in_g1_global,1) - points(neiID_cut2_in_g2_global,1))^2 + (points(neiID_cut2_in_g1_global,2) - points(neiID_cut2_in_g2_global,2))^2;
        d2 = sqrt(d2);
        [neiID_cut1_in_g1_global neiID_cut1_in_g2_global d1 dis1]
        [neiID_cut2_in_g1_global neiID_cut2_in_g2_global d2 dis2]
        
    else% P2也有问题，舍弃此次拆分，二合一但不加边
        disp('CASE 4: P1，P2都有问题，二合一但MI不加边。');
    end
end

disp('Func_Check_MI End');

end

