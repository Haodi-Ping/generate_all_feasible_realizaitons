function [ratio] = Func_Ratio(subgraph,distMatrix,ConnectivityM)
%Func_Ratio 此处显示有关此函数的摘要
%   此处显示详细说明
[sub_num,~] = size(subgraph);
ratio = zeros(1,sub_num);
for i=1:sub_num
    r=find(subgraph(i,:));
    [~,nnn]=size(r);
    edge=sum(sum(ConnectivityM(r,r)))/2;
    ratio(i)=(edge-2*nnn+3)*2/nnn/(nnn-1);
end

% 
% [sub_num,~] = size(subgraph);
% residual_error = zeros(1,sub_num);
% ratio = zeros(1,sub_num);
% % 暂定为
% max_residual_error = -1;
% for i = 1:sub_num
%     graph_sub = find(subgraph(i,:)>0);
%     [~,node_num] = size(graph_sub);
%     points_sub = points(graph_sub,:);
%     dist_sub = distMatrix(graph_sub,graph_sub);
%     ConnectivityM_sub = ConnectivityM(graph_sub,graph_sub);
%     edge_num = sum(sum(ConnectivityM_sub))/2;
%     pos_sub = Alg_ARAP( ConnectivityM_sub,dist_sub,points_sub);
%     %[~,regist_pos_sub,~]=procrustes(points_sub, pos_sub, 'Scaling', false);
%     %pos_super(graph_sub,:) = regist_pos_sub;
%     % plotpositionsaa(points_sub',regist_pos_sub',ConnectivityM_sub,100,'SubGraph');
%     Y = pdist(pos_sub);
%     dist_pos = squareform(Y);
%     residual_error(i) = sum(sum(abs(dist_pos - dist_sub)))/edge_num;
%     %residual_error(i) = sum(sum(abs(dist_pos - dist_sub)));
%     if (residual_error(i)>max_residual_error)
%         max_residual_error = residual_error(i);
%     end
% end
% residual_error
% for i = 1: sub_num
%     if(residual_error(i) == 0)
%         residual_error(i) = residual_error(i)+0.1;
%     end
%     ratio(i) = exp(max_residual_error/residual_error(i));
% end

end

