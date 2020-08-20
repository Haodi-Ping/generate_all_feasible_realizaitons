function [points, ConnectivityM, distMatrix,nClean] = Func_FlipNetwork...
    (randstate, radius, nPoints,nf,nf_model)
%Func_FlipNetwork 通过修补，生成有flip，但无关节节点的拓扑。
%   %1. 对于单个的关节节点，将图分成2块，那么从2块中各取1个邻居，然后在此邻居中点放置一个新节点；
%   % 无需此步-2. 对于桥节点，将图分成2块，那么桥的2个端点各取1个邻居，然后在此邻居中间放置一个新节点，若距离大于2r，可放置多个；
%   3. 更新 ConnectivityM 和 distMatrix，因为新加的节点可能和别的点也有边。
% Boxscale = 100;
% disp('Begin Func_SuperFlipNetwork')
if isempty(nf_model)
    nf_model = 'additive';
end
%% 1. 随机生成初始网络
% disp('------ Func_SuperFlipNetwork: I. 随机生成初始网络 ------------------------------');
rand('seed',randstate);
points=rand(2,nPoints)-0.5;
points=points';

%%
% hole_index = 0;
% hole_pnts = [];
% thrd = 0.15;
% for i = 1:nPoints
%     xx = points(i,1);
%     yy = points(i,2);
%     if (xx>-thrd&&xx<thrd&&yy>-thrd&&yy<thrd)
%         hole_index = hole_index+1;
%         hole_pnts(hole_index) = i
%     end
% end
% [~,nn] = size(hole_pnts)
% nPoints = nPoints-nn;
% points(hole_pnts,:) = [];
%%
%points = points*0.6;
Y=pdist(points);
trueDistMatrix=squareform(Y);
% 生成邻接矩阵
ConnectivityM=zeros(nPoints, nPoints);
for i=1:nPoints
    for j=1:nPoints
        if trueDistMatrix(i,j)>radius
            ConnectivityM(i,j)=0;
        elseif i==j
            ConnectivityM(i,j)=0;
        else
            ConnectivityM(i,j)=1;
        end
    end
end
% sprintf('随机生成初始网络,节点数： %d ',nPoints)
% plotgraph([],points'*Boxscale,ConnectivityM,'1.初始网络');
%% 2. 删除悬挂点和孤岛点
% disp('删除悬挂点')
% disp('------ Func_SuperFlipNetwork: II. 删除悬挂点和孤岛点 ------------------------------');
degree = zeros(1,nPoints);
flag = 1;
while flag
    flag = 0;
    [m,~] = size(ConnectivityM);
    for i = 1: m
        [~,degree(i)] = size(find(ConnectivityM(i,:)>0));
        %degree(i) = degree(i) - 1;
    end
    island = find(degree < 2);% 度为0是孤岛，度为1是悬挂
    [~,nIsland] = size(island);
    if nIsland > 0
        flag = 1;
        ConnectivityM(island,:) = [];
        ConnectivityM(:,island) = [];
        trueDistMatrix(island,:) = [];
        trueDistMatrix(:,island) = [];
        %distMatrix(island,:) = [];
        %distMatrix(:,island) = [];
        points(island,:) = [];
        degree(island) =[];
        %points(:,island) = [];
    end
end
nClean = m;
island_num = nPoints-nClean;
% sprintf('删除悬挂点,节点数： %d ',island_num)
%figure
%plotgraph([],points'*Boxscale,ConnectivityM);

% disp("nIsland:"+(nPoints-nClean));

%plotgraph([],points'*Boxscale,ConnectivityM,'2.删除悬挂点');

%% 3. 如果初始网络不连通，取最大连通子图
% disp('如果初始网络不连通，取最大连通子图')
% disp('------ Func_SuperFlipNetwork: III. 如果初始网络不连通，取最大连通子图 ------------------------------');
sub_main=Func_SuperBFS(ConnectivityM);  %计算 去掉当前 cut set 后 连通子图
[sub_num,~]=size(sub_main);%
disjoint_num = 0;
if(sub_num>1)
    max = 0;
    max_index = 0;
    for i = 1:sub_num
        [~,current_num] = size(find(sub_main(i,:)>0));
        if (current_num>max)
            max = current_num;
            max_index = i;
        end
    end
    max_sub = sub_main(max_index,:);
    max_nodes = find(max_sub>0);
    points = points(max_nodes,:);
    ConnectivityM = ConnectivityM(max_nodes,max_nodes);
    trueDistMatrix = trueDistMatrix(max_nodes,max_nodes);
    disjoint_num = nClean-max;
    nClean = max;
end
% sprintf('删除孤岛子图,节点数： %d ',disjoint_num)
%plotgraph([],points'*Boxscale,ConnectivityM,'3.删除孤岛子图');

%% 4. 解决关节节点
% 找出关节节点
%arthro = []; % 记录找到的关节节点
% disp('解决关节节点')
% disp('------ Func_SuperFlipNetwork: IV. 解决关节节点 ------------------------------');
arthro_num = 0;
graph = ones(1,nClean);
for i = 1:nClean
    %i
    tempID = graph; % 操作针对备份的 graph
    tempID(i)=0; % 从 graph 中删除当前点
    subtempID = find(tempID>0); % subtempID 是去掉当前 "关节点" 后剩下的点
    subtempM = ConnectivityM(subtempID,subtempID); %从原图中取出 去掉当前点后 的子图关系 --- 子图节点数 × 子图节点数
    tempresult=Func_SuperBFS(subtempM);  %计算 去掉当前 cut set 后 连通子图
    [tempa,~]=size(tempresult);% a 为其个数，此处经测试，全部为 2，1个separator不会有 >2 个模块拆出来
    if (tempa>1)
        neis_i = find(ConnectivityM(i,:)>0);
        graph_1 = find(tempresult(1,:)>0); % 切出来的第一个子图
        graph_2 = find(tempresult(2,:)>0);
        [~,graph_1_num] = size(graph_1);
        [~,graph_2_num] = size(graph_2);
        % 删除完第i个点后，BFS结果中，得到的subgraph节点编号>i的点，相对于原图来说，是减了1的，所以要加回来
        for j = 1:graph_1_num
            if(graph_1(j)>=i)
                graph_1(j) = graph_1(j)+1;
            end
        end
        for k = 1:graph_2_num
            if(graph_2(k)>=i)
                graph_2(k) = graph_2(k)+1;
            end
        end
        neis_1 = intersect(neis_i,graph_1);
        neis_2 = intersect(neis_i,graph_2);
        nei1 = points(neis_1(1),:);
        nei2 = points(neis_2(1),:);
        % pnt =[ (nei1(1)+nei2(1))/2  (nei1(2)+nei2(2))/2];
        % disp('当前节点是关节节点 添加辅助点')
        arthro_num = arthro_num+1;
        points_auxiliary(arthro_num,1) = (nei1(1)+nei2(1))/2 ;
        points_auxiliary(arthro_num,2) = (nei1(2)+nei2(2))/2 ;
        %arthro(arthro_num) = i;
    end
end
if(arthro_num>0)
    nClean = nClean+arthro_num;
    [points, ConnectivityM, trueDistMatrix] = update_mats(points,radius,points_auxiliary);
end
% plotgraph([],points'*100,ConnectivityM,'整体调整完毕');%绘00制拓扑图
%plotgraph([],points'*Boxscale,ConnectivityM,'4.调整关节节点');

% sprintf('解决关节节点,节点数： %d ',arthro_num)
%arthro


%% 6. 所有修补均针对Ground Truth，最后一步才加噪声
distMatrix = zeros(nClean,nClean);
for i=1:nClean
    for j=i:nClean
        if i==j
            distMatrix(i,j)=0;
            continue
        end
        if ConnectivityM(i,j)==1
            if strcmp(nf_model,'additive')
            dis=trueDistMatrix(i,j)+nf*randn();
            end
            if strcmp(nf_model,'multiplicative')
            dis=trueDistMatrix(i,j)*(1+nf*randn());
            end
            if dis>0
                distMatrix(i,j)=dis;
                distMatrix(j,i)=dis;
            else
                distMatrix(i,j)=trueDistMatrix(i,j);
                distMatrix(j,i)=trueDistMatrix(i,j);
            end
        else
            distMatrix(i,j)=NaN;
            distMatrix(j,i)=NaN;
        end
    end
end
%sprintf('初始网络:%d,删除悬挂点:%d,删除孤岛子图:%d,添加辅助点:%d,最终网络:%d.',nPoints,island_num,disjoint_num,arthro_num,nClean)
% disp('End Func_SuperFlipNetwork')


end


%% update conn dis
function [points, ConnectivityM, trueDistMatrix] = update_mats(points,radius,points_auxiliary)
[num,~] = size(points);
[num1,~] = size(points_auxiliary);
nPoints = num+num1;
points(num+1:num+num1,:) = points_auxiliary;
Y=pdist(points);
trueDistMatrix=squareform(Y);
% 生成邻接矩阵
ConnectivityM=zeros(nPoints, nPoints);
for i=1:nPoints
    for j=1:nPoints
        if trueDistMatrix(i,j)>radius
            ConnectivityM(i,j)=0;
        elseif i==j
            ConnectivityM(i,j)=0;
        else
            ConnectivityM(i,j)=1;
        end
    end
end

end
