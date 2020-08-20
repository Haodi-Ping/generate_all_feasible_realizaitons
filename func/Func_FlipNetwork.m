function [points, ConnectivityM, distMatrix,nClean] = Func_FlipNetwork...
    (randstate, radius, nPoints,nf,nf_model)
%Func_FlipNetwork ͨ���޲���������flip�����޹ؽڽڵ�����ˡ�
%   %1. ���ڵ����Ĺؽڽڵ㣬��ͼ�ֳ�2�飬��ô��2���и�ȡ1���ھӣ�Ȼ���ڴ��ھ��е����һ���½ڵ㣻
%   % ����˲�-2. �����Žڵ㣬��ͼ�ֳ�2�飬��ô�ŵ�2���˵��ȡ1���ھӣ�Ȼ���ڴ��ھ��м����һ���½ڵ㣬���������2r���ɷ��ö����
%   3. ���� ConnectivityM �� distMatrix����Ϊ�¼ӵĽڵ���ܺͱ�ĵ�Ҳ�бߡ�
% Boxscale = 100;
% disp('Begin Func_SuperFlipNetwork')
if isempty(nf_model)
    nf_model = 'additive';
end
%% 1. ������ɳ�ʼ����
% disp('------ Func_SuperFlipNetwork: I. ������ɳ�ʼ���� ------------------------------');
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
% �����ڽӾ���
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
% sprintf('������ɳ�ʼ����,�ڵ����� %d ',nPoints)
% plotgraph([],points'*Boxscale,ConnectivityM,'1.��ʼ����');
%% 2. ɾ�����ҵ�͹µ���
% disp('ɾ�����ҵ�')
% disp('------ Func_SuperFlipNetwork: II. ɾ�����ҵ�͹µ��� ------------------------------');
degree = zeros(1,nPoints);
flag = 1;
while flag
    flag = 0;
    [m,~] = size(ConnectivityM);
    for i = 1: m
        [~,degree(i)] = size(find(ConnectivityM(i,:)>0));
        %degree(i) = degree(i) - 1;
    end
    island = find(degree < 2);% ��Ϊ0�ǹµ�����Ϊ1������
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
% sprintf('ɾ�����ҵ�,�ڵ����� %d ',island_num)
%figure
%plotgraph([],points'*Boxscale,ConnectivityM);

% disp("nIsland:"+(nPoints-nClean));

%plotgraph([],points'*Boxscale,ConnectivityM,'2.ɾ�����ҵ�');

%% 3. �����ʼ���粻��ͨ��ȡ�����ͨ��ͼ
% disp('�����ʼ���粻��ͨ��ȡ�����ͨ��ͼ')
% disp('------ Func_SuperFlipNetwork: III. �����ʼ���粻��ͨ��ȡ�����ͨ��ͼ ------------------------------');
sub_main=Func_SuperBFS(ConnectivityM);  %���� ȥ����ǰ cut set �� ��ͨ��ͼ
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
% sprintf('ɾ���µ���ͼ,�ڵ����� %d ',disjoint_num)
%plotgraph([],points'*Boxscale,ConnectivityM,'3.ɾ���µ���ͼ');

%% 4. ����ؽڽڵ�
% �ҳ��ؽڽڵ�
%arthro = []; % ��¼�ҵ��Ĺؽڽڵ�
% disp('����ؽڽڵ�')
% disp('------ Func_SuperFlipNetwork: IV. ����ؽڽڵ� ------------------------------');
arthro_num = 0;
graph = ones(1,nClean);
for i = 1:nClean
    %i
    tempID = graph; % ������Ա��ݵ� graph
    tempID(i)=0; % �� graph ��ɾ����ǰ��
    subtempID = find(tempID>0); % subtempID ��ȥ����ǰ "�ؽڵ�" ��ʣ�µĵ�
    subtempM = ConnectivityM(subtempID,subtempID); %��ԭͼ��ȡ�� ȥ����ǰ��� ����ͼ��ϵ --- ��ͼ�ڵ��� �� ��ͼ�ڵ���
    tempresult=Func_SuperBFS(subtempM);  %���� ȥ����ǰ cut set �� ��ͨ��ͼ
    [tempa,~]=size(tempresult);% a Ϊ��������˴������ԣ�ȫ��Ϊ 2��1��separator������ >2 ��ģ������
    if (tempa>1)
        neis_i = find(ConnectivityM(i,:)>0);
        graph_1 = find(tempresult(1,:)>0); % �г����ĵ�һ����ͼ
        graph_2 = find(tempresult(2,:)>0);
        [~,graph_1_num] = size(graph_1);
        [~,graph_2_num] = size(graph_2);
        % ɾ�����i�����BFS����У��õ���subgraph�ڵ���>i�ĵ㣬�����ԭͼ��˵���Ǽ���1�ģ�����Ҫ�ӻ���
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
        % disp('��ǰ�ڵ��ǹؽڽڵ� ��Ӹ�����')
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
% plotgraph([],points'*100,ConnectivityM,'����������');%��00������ͼ
%plotgraph([],points'*Boxscale,ConnectivityM,'4.�����ؽڽڵ�');

% sprintf('����ؽڽڵ�,�ڵ����� %d ',arthro_num)
%arthro


%% 6. �����޲������Ground Truth�����һ���ż�����
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
%sprintf('��ʼ����:%d,ɾ�����ҵ�:%d,ɾ���µ���ͼ:%d,��Ӹ�����:%d,��������:%d.',nPoints,island_num,disjoint_num,arthro_num,nClean)
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
% �����ڽӾ���
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
