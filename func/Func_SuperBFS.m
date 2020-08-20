function [subResult] = Func_SuperBFS(Cmartix)
%Cmartix
%Func_SuperBFS 猜测BFS会卡在有孤岛点的图(度为0的节点)
%   故先检测孤岛点，没有孤岛点，再去调用BFS
[total_num,~] = size(Cmartix);

for i = 1: total_num
    %find(Cmartix(i,:)>0)
    [~,degree(i)] = size(find(Cmartix(i,:)>0));
end
LonelyNode = find(degree == 0);% 度为0是孤岛，度为1是悬挂
[~,LonelyNum] = size(LonelyNode);
if (LonelyNum == 0)
    subResult = BFS(Cmartix);
else
    % disp('使用BFS_Refined')
    % 有孤岛点的处理方法，先拿走孤岛点，然后调用BFS，再将孤岛点放到BFS结果
    %LonelyNode
    %LonelyNum
    Cmartix(LonelyNode,:) = [];
    Cmartix(:,LonelyNode) = [];
    subResult = BFS( Cmartix );
    [subNum,totalNum] = size(subResult);
    % 孤岛点放到BFS结果  b=[a(:,1),zeros(5,1),a(:,2:5)]
    for i = 1:LonelyNum
        subResult = [subResult(:,1:LonelyNode(i)-1) zeros(subNum,1) subResult(:,LonelyNode(i):totalNum)];
        totalNum = totalNum+1;
    end
    for i = 1:LonelyNum
        a = zeros(1,totalNum);
        a(LonelyNode(i)) = 1;
        subResult = [subResult;a];
    end    
end
[subNum,totalNum] = size(subResult);
if subNum>2
    subResult_new = zeros(2,totalNum);
    sub_nodes = sum(subResult,2);
    [~,max_index] = max(sub_nodes);
    subResult_new(1,:) = subResult(max_index,:);
    subResult(max_index,:) = [];
    subResult_new(2,:) = sum(subResult,1);
    subResult = subResult_new;
end
% for i = 1: total_num
%     find(Cmartix(i,:)>0)
% end

end

