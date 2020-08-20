function [separator,QI_connt,subgraph,cut_set,ConnectivityM,distMatrix] = ...
    Func_Partition(graph,radius,points,ConnectivityM,distMatrix)
%Func_Partition graph 1行n列，为1的元素代表此图包含此节点
%   points 只有 MI check 需要
%% 应该将 QI 和 MI 加的边存上，返回一个新的 ConnectivityM 和 distMatrix 矩阵
%% 声明输出结果
separator = 0; % 记录separator的个数
QI_connt = 0;
%MI_count = 0;
subgraph = []; % 存放分割结果
cut_set = []; % 存放此分割结果的 cut set
% disp('Func_Partition')
if(sum(graph) == 0)
    % 空图 do nothing
    % disp('------ Func_Partition: 空图 do nothing ------------------------------');
else
    %% 预处理，找到关节节点，是不能用来当cut set的
    % disp('------ Func_Partition: 预处理，关节节点不当 cut set ------------------------------')
    arthro = []; % 记录找到的关节节点
    arthro_num = 0;
    % [~,num] = size(ConnectivityM);
    ID = find (graph>0);
    [~,nodes_num] = size(ID);
    for i = 1:nodes_num
        tempID = graph; % 操作针对备份的 graph
        tempID(ID(i))=0; % 从 graph 中删除当前点
        subtempID = find(tempID>0); % subtempID 是去掉当前 "关节点" 后剩下的点
        subtempM = ConnectivityM(subtempID,subtempID); %从原图中取出 去掉当前 cut set 后 的子图关系 --- 子图节点数 × 子图节点数
        tempresult=Func_SuperBFS(subtempM);  %计算 去掉当前 cut set 后 连通子图
        [tempa,~]=size(tempresult);% a 为其个数，此处经测试，全部为 2，1个separator不会有 >2 个模块拆出来
        if (tempa>1)
            arthro_num = arthro_num+1;
            arthro(arthro_num) = ID(i);
        end
    end
    %arthro
    %% 开始分割
    %disp('------ Func_Partition: 开始分割 ------------------------------');
    state=0; % 标记是否成功
    [~,num] = size(ConnectivityM);
    ID = find (graph>0);
    [~,nodes_num] = size(ID);
    for i=1:nodes_num
        % 如果 i 是关节节点，跳过
        is_i_arthro = size(find(arthro == ID(i)));
        if (is_i_arthro>0)
            continue;
        end
        
        % 如果 j 是关节节点，跳过
        for j=i+1:nodes_num
            is_j_arthro = size(find(arthro == ID(j)));
            if (is_j_arthro>0)
                continue;
            end
            % 开始分割
            
            markID=graph;% 首先markID 初始化为传入的图 1行npoints列 --- 1 × 100
            markID(ID(i))=0; %
            markID(ID(j))=0; % 尝试删除ID(i)与ID（j）个节点，即尝试 cut set = {v_ID(i),v_ID(j) }
            if ConnectivityM(ID(i),ID(j))==0
                continue
            end
            subID=find(markID>0); % subID 是去掉当前 cut set 后剩下的点
            subM=ConnectivityM(subID,subID); %从原图中取出 去掉当前 cut set 后 的子图关系 --- 子图节点数 × 子图节点数
            result=Func_SuperBFS(subM);  %计算 去掉当前 cut set 后 连通子图
            [a,~]=size(result); % a 为其个数，此处经测试，全部为 2，1个separator不会有 >2 个模块拆出来
            
            if (a>1)
                %Pending_Cut_Set=[ID(i) ID(j)]
                % 找到一个 separator
                % 检查QI，解决不了，才往最终结果加
                [ConnectivityM,distMatrix,QI_flag] = ...
                    Func_Check_QI(ID(i),ID(j),ConnectivityM,distMatrix,radius,points);
                %flag_arthro = Func_Find_Arthro(result,ConnectivityM);
                %if (QI_flag == 0 && flag_arthro == 0)
                separator = separator+1;
                if (QI_flag == 0 )
                    % QI 失败,划分结束
                    % MI_flag = Func_Check_MI(cut_set,subgraph,radius,points,distMatrix,ConnectivityM)
                    subgraph=zeros(a,num); % a 为切出来的子图数量，num为节点数
                    cut_set = [ID(i),ID(j)];
                    for k=1:a % 对于切出来的每个子图,存到 subgraph 中
                        resultID=find(result(k,:)>0); % 将当前子图包含的节点存到 resultID
                        subgraph(k,subID(resultID))=1; % 将当前子图包含的节点存到 todoM
                        subgraph(k,ID(i))=1;
                        subgraph(k,ID(j))=1; % 将  cut set 存到 todoM
                        %subG = find(subgraph(k,:)>0)
                    end
%                     for ii = 1:a
%                         sub_ii = find(subgraph(ii,:)>0)
%                     end
                    state=1;
                    break;
                else
                    % 此 separator 被 QI 解决
                    % 更新 dismat 进行加边
                    QI_connt = QI_connt+1;
                    continue;
                end
                %
                
                
                %             num1 = sum(subgraph(1,:));
                %             num2 = sum(subgraph(2,:));
                %             if (num1<=3 || num2<=3) % 有3个点的图，直接存
                %                 state=1;
                %                 break;
                %             end
                
                
                
                %
                
                
            end
        end
        if(state==1)% 只要成功找到一对  cut set 就停止找 i
            break;
        end
    end
    
end

end

