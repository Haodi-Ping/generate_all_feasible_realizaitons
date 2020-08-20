function [residual_lst,mean_residual] = F_residual(pos,points,Boxscale)
%F_residual ����ʵ�ֽ���� ground truth �Ĳв� 2019-9-18 16:08:21
%   ÿ����Ĳв�浽������
num_points = size(points,1);
residual_lst = zeros(1,num_points);
for i=1:num_points
    residual_lst(i) = abs(norm(pos(i,:) - points(i,:))*Boxscale);
end
mean_residual = mean(residual_lst);
end

