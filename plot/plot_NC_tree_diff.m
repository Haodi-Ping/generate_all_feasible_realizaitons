function plot_NC_tree_diff...
    (NC_Tree,NC_Tree_level,points,distMatrix,ConnectivityM,cut_sets,leaf_marks,leaf_realizations)
%PLOT_NC_TREE show NC-Tree
%   此处显示详细说明
Boxscale = 100;
fig_width = 30;
% fig_height = 15;
% node_gap = 1; % 节点之间的间距
% sub_size = fig_width/(2^(NC_Tree_level-1));
colNum = 2^NC_Tree_level;
rowNum = NC_Tree_level+1;
figure
for i = 0:NC_Tree_level%NC_Tree_level:-1:0
    this_level_num = 2^i; % the node number of i_th level of NC-Tree
    for j = 1:this_level_num
        %hold on
        this_tree_index = 2^i + j - 1;
        this_num = sum(NC_Tree(this_tree_index,:));
        if this_num == 0
            continue
        end        
        %this_plot_index = i * colNum + j;
        this_occupation_num = colNum/2^i;
        this_plot_index = i * colNum + (j-1) * this_occupation_num + 1;
        subplot(rowNum,colNum,...
            this_plot_index:(this_plot_index + this_occupation_num - 1));
        %this_subgraph = find(NC_Tree(this_tree_index,:));  
        leaf_mark = 0;
        if (leaf_marks(this_tree_index) == 1)
            leaf_mark = 1;
        end
        try 
            cut_set = cut_sets(this_tree_index,:);
        catch
            cut_set = [0 0];
        end
        id = find(NC_Tree(this_tree_index,:));
        sub_num = sum(NC_Tree(this_tree_index,:));
        sub_dist = distMatrix(id,id);
        sub_conn = ConnectivityM(id,id);
        sub_points = points(id,:);
        coord_arap = Alg_ARAP(sub_conn,sub_dist,sub_points);
        [sub_neighbors] = connectivity2neighbors(sub_conn);
        [coord_wcs] = Alg_WCS...
            (sub_num,sub_points,coord_arap,sub_neighbors,sub_conn,sub_dist,Boxscale);
        %this_coord = leaf_realizations{}
         [~,coord_wcs,~] = procrustes(sub_points, coord_wcs, 'Scaling', false);
        plot_NC_node_diff(sub_points',coord_wcs',sub_conn,Boxscale,'');
    end
end

%set(gcf,'unit','centimeters','position',[2 4 fig_width fig_height]);
end
