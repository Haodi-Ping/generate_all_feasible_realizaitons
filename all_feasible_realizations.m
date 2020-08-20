%% 生成全部 possible 实现， 并根据 CFFC 找出 feasible 实现
clc
clear all
close all
warning('off');
currentFolder = pwd;
addpath(genpath(currentFolder));
%% 初始化
switch_debug = 0;
switch_plots = 1;
randstate = 16;
npoints = 50;
radius = 0.16;
switch_CFFC = 0;
nf = 0;
Boxscale = 100;
[points, ConnectivityM, distMatrix,nClean] = ...
    Func_FlipNetwork(randstate, radius, npoints,nf,'');
average_degree = sum(sum(ConnectivityM))/nClean;

ground_truth = points; % ground truth
realization_init = zeros(nClean,2,1);
realization_flip = zeros(nClean,2,1); % realization_flip 与 realization_init 一一对应
realization_possible = zeros(nClean,2,1); % 全部 possible
realization_feasible = zeros(nClean,2,1); % realization_init 与 realization_flip 的并集

if switch_plots
    plotgraph([],points'*Boxscale,ConnectivityM,'Network Topo');
end
nbrNums = sum(ConnectivityM,2);
maxNbrNum = max(nbrNums);
neighbors = zeros(nClean, maxNbrNum);
neighbors(:,1)=nbrNums;
for i=1:nClean
    cnt = 1;
    for j=1:nClean
        if ConnectivityM(i,j)==1
            cnt=cnt+1;
            neighbors(i,cnt)=j;
        end
    end
end
%% records
ConnectivityM_cpy = ConnectivityM;
distMatrix_cpy = distMatrix;
num_separator = 0;
num_LFFC = 0;

%% LFFC and NC-Tree Construction
tic
graph = ones(1,nClean); % 原图为 1 行 nClean 列的全1数组
NC_Tree(1,:) = graph;
NC_Tree_level = 0;
is_current_level_partition_able = 1;
while (is_current_level_partition_able)
    is_current_level_partition_able = 0;
    for i = 2^NC_Tree_level:(2^(NC_Tree_level+1)-1) % 遍历当前level
        [separator,QI_count,temp_subgraph,temp_cut_set,ConnectivityM,distMatrix] = ...
            Func_Partition(NC_Tree(i,:),radius,points,ConnectivityM,distMatrix);
        num_LFFC = num_LFFC + QI_count;
        num_separator = num_separator + separator;
        [a,~] = size(temp_subgraph);
        NC_leaf_marks(i) = 0;
        if(a>1)
            NC_cut_sets(i,:) = temp_cut_set;
            NC_Tree(2*i,:)= temp_subgraph(1,:);
            NC_Tree(2*i+1,:)= temp_subgraph(2,:);
            is_current_level_partition_able = 1;
        else
            NC_Tree(2*i,:)= zeros(1,nClean);
            NC_Tree(2*i+1,:)= zeros(1,nClean);
            if (sum(NC_Tree(i,:))>0)
                NC_leaf_marks(i) = 1;
            end
        end
    end
    if(is_current_level_partition_able == 1)
        NC_Tree_level = NC_Tree_level+1;
    else
        %sprintf('终止partition',NC_Tree_level+1)
    end
end

if NC_Tree_level==0
    if QI_count == 0
        disp('The input graph is already 3-connected ...')
    else
        disp('All separators are solved by LFFC.')
    end
    %return
end
NC_Tree(2^(NC_Tree_level+1):2^(NC_Tree_level+2)-1,:) = [];
time_LFFC = toc;
if switch_plots
    try
        plot_NC_tree(NC_Tree,NC_Tree_level,points,ConnectivityM,...
            NC_cut_sets,NC_leaf_marks,'NC-Tree');
    catch
        disp('I can not plot NC-Tree ...')
    end
end
%%
% [ARAP_Pos] = Alg_ARAP(ConnectivityM,distMatrix,points);
% [~,Z,~]=procrustes(points, ARAP_Pos, 'Scaling', false);
% ARAP_Pos=Z;
% plotpositionsaa(points',ARAP_Pos',ConnectivityM,Boxscale,'ARAP-LFFC');
%% record leaf nodes
tic
% subgraph = zeros(1,nClean);
% sub_index = 0;
% for i = 1:(2^(NC_Tree_level+1)-1)
%     if(NC_leaf_marks(i) == 1)
%         sub_index = sub_index+1;
%         subgraph(sub_index,:) = NC_Tree(i,:);
%     end
% end
%% local realization of leaf nodes
NC_num = 2^(NC_Tree_level+1)-1;
NC_realizations = cell(1,NC_num);
Poss_realizations = cell(1,NC_num);
for j = 1:NC_num
    if(NC_leaf_marks(j) == 1)
        id = find(NC_Tree(j,:));
        sub_num = sum(NC_Tree(j,:));
        sub_dist = distMatrix(id,id);
        sub_conn = ConnectivityM(id,id);
        sub_points = points(id,:);
        coord_arap = Alg_ARAP(sub_conn,sub_dist,sub_points);
        [~,coord_arap,~ ]=procrustes(sub_points, coord_arap, 'Scaling', false);
        [sub_neighbors] = connectivity2neighbors(sub_conn);
        [coord_wcs] = Alg_WCS...
            (sub_num,sub_points,coord_arap,sub_neighbors,sub_conn,sub_dist,Boxscale);
        %SubgraphCoord(2*j-1:2*j,id)  = coord_wcs';
        NC_realizations{j} = coord_wcs;
        Poss_realizations{j} = coord_wcs;
    else
        NC_realizations{j} = NaN;
        Poss_realizations{j} = NaN;
    end
end
time_localozation = toc;
% if switch_plots
%     plot_NC_tree_diff...
%         (NC_Tree,NC_Tree_level,points,distMatrix,ConnectivityM,NC_cut_sets,NC_leaf_marks,NC_realizations);
% end
%% CFFC checking
tic
pair_possible = 0;
num_CFFC = 0;
index = NC_Tree_level;
is_current_level_merge_able = 1;
while (switch_CFFC == 1 && index>=0 && is_current_level_merge_able)
    if(NC_leaf_marks(1) == 1)
        disp('NC-Tree only contains a root node and node.graph is 3-connected.');
        break;
    end
    sprintf('merging %d level...',index+1)
    %is_current_level_merge_able = 0; % level都不能合并，那么level-1就不用尝试合并了
    for i = 2^index:2:(2^(index+1)-1) % 遍历当前level
        if (sum(NC_Tree(i,:))==0) % 空图不处理
            continue;
        end
        if (NC_leaf_marks(i)==0 || NC_leaf_marks(i+1)==0) % 非叶子节点不处理,仅当2个叶子节点才处理；
            continue;
        end
        disp('try to merge...');
        subgraph_merging(1,:) = NC_Tree(i,:);
        subgraph_merging(2,:) = NC_Tree(i+1,:);
        this_fartherID = floor(i/2);
        cut_set = NC_cut_sets(this_fartherID,:);
        pnt_1 = find(subgraph_merging(1,:)>0);
        pnt_2 = find(subgraph_merging(2,:)>0);
        r_l = NC_realizations{i};
        r_r = NC_realizations{i+1};
        pnt_merge = union(pnt_1,pnt_2);
        points_merge = points(pnt_merge,:);
        conn_merge = ConnectivityM(pnt_merge,pnt_merge);
        %         [r_merge, flag_CFFC, ID_merge] = Func_CFFC_iterative_merge...
        %             (cut_set,subgraph_merging,radius,r_l,r_r,ConnectivityM,distMatrix);
        [r_merge, flag_CFFC, ID_merge] = Func_CFFC...
            (cut_set,subgraph_merging,radius,r_l,r_r);
        if switch_debug
            plotpositionsaa(points_merge',r_merge{1}',conn_merge,Boxscale,...
                ['merge-init ' num2str(flag_CFFC(1))]);
            plotpositionsaa(points_merge',r_merge{2}',conn_merge,Boxscale,...
                ['merge-flip ' num2str(flag_CFFC(2))]);
        end
        %Poss_realizations{this_fartherID} = [r_merge{1} r_merge{2}];
        %flag_CFFC
        switch sum(flag_CFFC)
            case 0 % both of the realizations are feasible
                NC_realizations{this_fartherID} = [r_merge{1} r_merge{2}];
                %NC_Tree(i,:) = zeros(1,nClean);
                %NC_Tree(i+1,:) = zeros(1,nClean);
                NC_Tree(this_fartherID,ID_merge) = 1;
                % NC_leaf_marks(i) = 1;
                % NC_leaf_marks(i+1) = 1;
                % NC_leaf_marks(this_fartherID) = 2;
                is_current_level_merge_able = 0;
            case 1 % one of the realizations is feasible
                num_CFFC = num_CFFC+1;
                pos_merge = (1-flag_CFFC(1))*r_merge{1}+(1-flag_CFFC(2))*r_merge{2};
                NC_realizations{this_fartherID} = pos_merge;
                NC_Tree(i,:) = zeros(1,nClean);
                NC_Tree(i+1,:) = zeros(1,nClean);
                NC_Tree(this_fartherID,ID_merge) = 1;
                NC_leaf_marks(i) = 0;
                NC_leaf_marks(i+1) = 0;
                NC_leaf_marks(this_fartherID) = 1;
                NC_cut_sets(this_fartherID,:) = [0 0];
            case 2 % none of the realizations is feasible
                is_current_level_merge_able = 0;
                disp('no realization of the current components is feasible...')
        end
    end
    if(is_current_level_merge_able == 1)
        index = index - 1;
    end
end
if switch_plots
    try
        plot_NC_tree(NC_Tree,NC_Tree_level,points,ConnectivityM,...
            NC_cut_sets,NC_leaf_marks,'Simplified NC-Tree');
    catch
        disp('I can not plot Simplified NC-Tree ...')
    end
end
time_CFFC = toc;
%% record leaf graph
subgraph = zeros(1,nClean);
sub_index = 0;
for i = 1:(2^(NC_Tree_level+1)-1)
    if(NC_leaf_marks(i) == 1)
        sub_index = sub_index+1;
        subgraph(sub_index,:) = NC_Tree(i,:);
        this_coord = NC_realizations{i};
        this_id = find(NC_Tree(i,:));
        SubgraphCoord(2*sub_index-1:2*sub_index,this_id)  = this_coord';
    end
end
%% all feasible realizaitons
%NC_Node_num = sum(NC_Tree,2);
index_count = 0;
for i = NC_Tree_level:-1:1
    for j = 2^i:2:(2^(i+1)-1)
        if NC_leaf_marks(j)*NC_leaf_marks(j+1)>0
            %['merging ' num2str(j) ' and ' num2str(j+1)]
            r_ls = NC_realizations{j};
            r_rs = NC_realizations{j+1};
            this_fartherID = floor(j/2);
            cut_set = NC_cut_sets(this_fartherID,:);
            subgraph_merging(1,:) = NC_Tree(j,:);
            subgraph_merging(2,:) = NC_Tree(j+1,:);
            flag_both_unfeasible = 0;
            for ii = 1:NC_leaf_marks(j)
                for jj = 1:NC_leaf_marks(j+1)
                    r_l = r_ls(:,1+2*(ii-1):2+2*(ii-1));
                    r_r = r_rs(:,1+2*(jj-1):2+2*(jj-1));
                    [r_merge, flag_CFFC, ID_merge] = Func_CFFC...
                        (cut_set,subgraph_merging,radius,r_l,r_r);
                    index_count = index_count+1;
                    ID_l = find(NC_Tree(j,:));
                    ID_r = find(NC_Tree(j+1,:));
                    points_l = points(ID_l,:);
                    conn_l = ConnectivityM(ID_l,ID_l);
                    points_r = points(ID_r,:);
                    conn_r = ConnectivityM(ID_r,ID_r);
                    if switch_debug
                        plotpositionsaa(points_l',r_l',conn_l,Boxscale,...
                            ['r-l ' num2str(index_count)]);
                        plotpositionsaa(points_r',r_r',conn_r,Boxscale,...
                            ['r-r ' num2str(index_count)]);
                    end
                    %flag_CFFC
                    flag = sum(flag_CFFC);
                    points_merge = points(ID_merge,:);
                    conn_merge = ConnectivityM(ID_merge,ID_merge);
                    if switch_debug
                        plotpositionsaa(points_merge',r_merge{1}',conn_merge,Boxscale,...
                            ['init ' num2str(index_count) ';' num2str(flag_CFFC(1))]);
                        plotpositionsaa(points_merge',r_merge{2}',conn_merge,Boxscale,...
                            ['flip ' num2str(index_count) ';' num2str(flag_CFFC(2))]);
                        
                    end
                    %                     pair_possible = pair_possible+1;
                    %                     realization_init(:,:,pair_possible) = r_merge{1};
                    %                     realization_init(:,:,pair_possible) = r_merge{2};
                    %                     realization_possible(:,:,2*pair_possible-1) = r_merge{1};
                    %                     realization_possible(:,:,2*pair_possible) = r_merge{2};
                    if isnan(NC_realizations{this_fartherID})
                        Poss_realizations{this_fartherID} = [r_merge{1} r_merge{2}];
                    else
                        Poss_realizations{this_fartherID} = ...
                            [Poss_realizations{this_fartherID} r_merge{1} r_merge{2}];
                    end
                    switch flag
                        case 0 % both feasible
                            NC_leaf_marks(this_fartherID) =  ...
                                NC_leaf_marks(this_fartherID)+2;
                            if isnan(NC_realizations{this_fartherID})
                                NC_realizations{this_fartherID} = ...
                                    [r_merge{1} r_merge{2}];
                            else
                                NC_realizations{this_fartherID} = ...
                                    [NC_realizations{this_fartherID} r_merge{1} r_merge{2}];
                            end
                        case 1 % alternative feasible
                            NC_leaf_marks(this_fartherID) =  ...
                                NC_leaf_marks(this_fartherID)+1;
                            if isnan(NC_realizations{this_fartherID})
                                NC_realizations{this_fartherID} = ...
                                    (1-flag_CFFC(1))*r_merge{1}+(1-flag_CFFC(2))*r_merge{2};
                            else
                                NC_realizations{this_fartherID} = ...
                                    [NC_realizations{this_fartherID} ((1-flag_CFFC(1))*r_merge{1}+(1-flag_CFFC(2))*r_merge{2})];
                            end
                        case 2 % both unfeasible 节点merge，直接使用merge节点实现
                            %disp('both unfeasible!!!')
                            if flag_both_unfeasible==1
                                continue
                            end
                            flag_both_unfeasible=1;
                            merge_num = size(ID_merge,2);
                            merge_dist = distMatrix_cpy(ID_merge,ID_merge);
                            merge_conn = ConnectivityM_cpy(ID_merge,ID_merge);
                            merge_points = points(ID_merge,:);
                            merge_arap = Alg_ARAP(merge_conn,merge_dist,merge_points);
                            [~,merge_arap,~ ]=procrustes(merge_points, merge_arap, 'Scaling', false);
                            [merge_neighbors] = connectivity2neighbors(merge_conn);
                            [merge_wcs] = Alg_WCS(merge_num,merge_points,merge_arap,...
                                merge_neighbors,merge_conn,merge_dist,Boxscale);
                            
                            NC_leaf_marks(this_fartherID) =  ...
                                NC_leaf_marks(this_fartherID)+1;
                            if isnan(NC_realizations{this_fartherID})
                                NC_realizations{this_fartherID} = merge_wcs;
                            else
                                NC_realizations{this_fartherID} = ...
                                    [NC_realizations{this_fartherID} merge_wcs];
                            end
                    end
                end
            end
            %NC_leaf_marks
        end
    end
end
mean_residual_distribution = zeros(1,NC_leaf_marks(1));
root_P = NC_realizations{1};
for i = 1: NC_leaf_marks(1)
    this_P = root_P(:,2*i-1:2*i);
    realization_feasible(:,:,i) = this_P;
    [~,temp] = F_residual(this_P,points,Boxscale);
    mean_residual_distribution(i) = temp;
    if switch_plots
        plotpositionsaa(points',this_P',ConnectivityM,Boxscale,...
            ['Feasible:' num2str(i) ' of ' num2str(NC_leaf_marks(1))]);
    end
end
root_poss = Poss_realizations{1};
[~,num_poss] = size(root_poss);
num_poss = num_poss/2;
for i = 1:num_poss
    realization_possible(:,:,i) = root_poss(:,2*i-1:2*i);
end