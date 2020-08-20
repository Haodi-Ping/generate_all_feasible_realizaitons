function plot_NC_tree...
    (NC_Tree,NC_Tree_level,points,ConnectivityM,cut_sets,leaf_marks,graphlabel)
%PLOT_NC_TREE show NC-Tree
%   此处显示详细说明
Boxscale = 100;
% fig_width = 30;
% fig_height = 15;
% node_gap = 1; % 节点之间的间距
% sub_size = fig_width/(2^(NC_Tree_level-1));
colNum = 2^NC_Tree_level;
rowNum = NC_Tree_level+1;
h = figure('color','w');
set(h,'name',graphlabel,'Numbertitle','off')

for i = 0:NC_Tree_level % NC_Tree_level:-1:0
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
        plot_tree_node(NC_Tree(this_tree_index,:),points'*Boxscale,...
            ConnectivityM,num2str(this_num),cut_set,leaf_mark);
    end
end
%suptitle(graphlabel)
pause(0.5)
%set(gcf,'unit','centimeters','position',[2 4 fig_width fig_height]);
end

%% plot_tree_node
function plot_tree_node(subgraph,points,ConnectivityM,label,cut_set,leaf_mark)
% node = find(subgraph);
% sub_points = points(:,node);
% sub_ConnectivityM = ConnectivityM(node,node);
[~,npts] = size(ConnectivityM);
% plot points
r1 = 1;
r2 = 2;
% plot edges
for j = 2:npts
    idx = find(ConnectivityM(1:j,j));
    if ~isempty(idx)
        len = length(idx);
        Pj = points(:,j)*ones(1,len);
        hold on
        plot([Pj(r1,:); points(r1,idx)],[Pj(r2,:); points(r2,idx)],'g');  
    end
end
for i=1:npts
    if subgraph(i) == 1
        if ismember(i,cut_set)
            h = plot(points(r1,i),points(r2,i),'or','markersize',3);
        else
            h = plot(points(r1,i),points(r2,i),'ok','markersize',3);
        end
        %text(points(r1,i)+1,points(r2,i)+1,num2str(i));
        %set(h,'linewidth',2);
        hold on
    end
end
xlim([-60 60]);
ylim([-60 60]);
axis('square');
set(gca,'XTick',[])
set(gca,'YTick',[])
%set(gca,'Visible','off');

%set(gcf,'unit','centimeters','position',[2 4 sub_size sub_size]);
if sum(cut_set)>0
    label = [label ': ' num2str(cut_set)];
end
if leaf_mark == 1
    title(label,'Color','r');
    return
end
title(label);
end
