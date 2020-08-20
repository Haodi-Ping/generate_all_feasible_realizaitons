function [Pos_Flip] = Func_Flip(Pos_Origin,axis_start,axis_end)
%Func_Flip 将 Pos_Origin 沿着 （axis_start,axis_end） Flip
%   此处显示详细说明

% 对称轴2端点
% pnt_i = [0 0]
% pnt_j = [1 1]

% 转换到三维
pnt_i = [Pos_Origin(axis_start,:) 0];
pnt_j = [Pos_Origin(axis_end,:) 0];
% 平移到原点
pnt_ii = pnt_i - pnt_i;
pnt_jj = pnt_j - pnt_i;
% 计算对称轴的单位向量
n = [pnt_jj(1)-pnt_ii(1) pnt_jj(2)-pnt_ii(2)];
mod = (pnt_ii(1)-pnt_jj(1))^2+(pnt_ii(2)-pnt_jj(2))^2;
mod = sqrt(mod);
n = n/mod;
% 计算R矩阵
x = n(1);
y = n(2);
%z = 0;
R = [2*x*x-1 2*x*y 0;2*x*y 2*y*y-1 0;0 0 -1];
% 逐个翻折
[m,~] = size(Pos_Origin);
Pos_Flip = zeros(m,2);
for index = 1:m
    %temp = ARAP_Pos(i);
    %index
    %temp = [];
    temp = [Pos_Origin(index,:) 0];
    temp = temp - pnt_i;
    flip = temp*R;
    flip = flip + pnt_i;
    Pos_Flip(index,:) =  flip(:,1:2);
end
end

