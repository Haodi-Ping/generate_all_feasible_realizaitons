function [r_merge, flag, ID_merge] = Func_CFFC...
    (cut_set,subgraph_merging,radius,r_l,r_r)
%Func_CFFC 此处显示有关此函数的摘要
%   此处显示详细说明

graph_1 = subgraph_merging(1,:);
graph_2 = subgraph_merging(2,:);
num_1 = sum(graph_1);
num_2 = sum(graph_2);
ID_1 = find(graph_1);
ID_2 = find(graph_2);

ID_merge = union(ID_1,ID_2);
[~,ind_1] = ismember(ID_1,ID_merge);
[~,ind_2] = ismember(ID_2,ID_merge);
[~,ind_cut] = ismember(cut_set,ID_merge);

cut1_ID_in_graph1 = find (ID_1 == cut_set(1));
cut2_ID_in_graph1 = find (ID_1 == cut_set(2));
cut1_ID_in_graph2 = find (ID_2 == cut_set(1));
cut2_ID_in_graph2 = find (ID_2 == cut_set(2));

cut_realization_g1 = r_l([cut1_ID_in_graph1,cut2_ID_in_graph1],:);
cut_realization_g2 = r_r([cut1_ID_in_graph2,cut2_ID_in_graph2],:);
%% 根据 cut set 的 realization 计算变换矩阵
[d,cut_realization_g2_transform,transform]=procrustes...
    (cut_realization_g1, cut_realization_g2, 'Scaling', false);%, 'Scaling', false
c = repmat(transform.c(1,:),num_2,1);
T = transform.T;
b = transform.b;
r_r_transform = b*r_r*T + c;
%% 方式 1 对齐子图的 realization
r_merge_1(ind_1,:) = r_l;
r_merge_1(ind_2,:) = r_r_transform;
r_merge_1(ind_cut,:) = (cut_realization_g1+cut_realization_g2_transform)/2;
%% 方式 2 对齐子图的 realization - 方式 1 的 flip
[r_r_flip] = Func_Flip(r_r_transform,cut1_ID_in_graph2,cut2_ID_in_graph2);
r_merge_2(ind_1,:) = r_l;
r_merge_2(ind_2,:) = r_r_flip;
r_merge_2(ind_cut,:) = (cut_realization_g1+cut_realization_g2_transform)/2;
%% check CFFC 
r_l([cut1_ID_in_graph1,cut2_ID_in_graph1],:) = [];
r_r_transform([cut1_ID_in_graph2,cut2_ID_in_graph2],:) = [];
r_r_flip([cut1_ID_in_graph2,cut2_ID_in_graph2],:) = [];

expand_1 = repelem(r_l, (num_2 - 2),1);
expand_2_init = repmat(r_r_transform, (num_1 - 2), 1);
expand_2_flip = repmat(r_r_flip, (num_1 - 2), 1);
vect_init = expand_1 - expand_2_init;
vect_flip = expand_1 - expand_2_flip;
dist_init = sum(abs(vect_init).^2,2).^(1/2);
dist_flip = sum(abs(vect_flip).^2,2).^(1/2);
% ind_cut_in_merge_1 = num_2*(cut1_ID_in_graph1-1) + cut1_ID_in_graph2;
% ind_cut_in_merge_2 = num_2*(cut1_ID_in_graph1-1) + cut2_ID_in_graph2;
% ind_cut_in_merge_3 = num_2*(cut2_ID_in_graph1-1) + cut1_ID_in_graph2;
% ind_cut_in_merge_4 = num_2*(cut2_ID_in_graph1-1) + cut2_ID_in_graph2;
% ind_cut_in_merge = ...
%     [ind_cut_in_merge_1 ind_cut_in_merge_2 ind_cut_in_merge_3 ind_cut_in_merge_4];
% dist_init(ind_cut_in_merge) = [];
% dist_flip(ind_cut_in_merge) = [];
flag_1 = 0;
flag_2 = 0;
if min(dist_init) <= radius
    flag_1 = 1;
end
if min(dist_flip) <= radius
    flag_2 = 1;
end
%% save result
r_merge{1} = r_merge_1;
r_merge{2} = r_merge_2;
flag = [flag_1;flag_2];
end