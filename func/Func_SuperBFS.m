function [subResult] = Func_SuperBFS(Cmartix)
%Cmartix
%Func_SuperBFS �²�BFS�Ῠ���йµ����ͼ(��Ϊ0�Ľڵ�)
%   ���ȼ��µ��㣬û�йµ��㣬��ȥ����BFS
[total_num,~] = size(Cmartix);

for i = 1: total_num
    %find(Cmartix(i,:)>0)
    [~,degree(i)] = size(find(Cmartix(i,:)>0));
end
LonelyNode = find(degree == 0);% ��Ϊ0�ǹµ�����Ϊ1������
[~,LonelyNum] = size(LonelyNode);
if (LonelyNum == 0)
    subResult = BFS(Cmartix);
else
    % disp('ʹ��BFS_Refined')
    % �йµ���Ĵ������������߹µ��㣬Ȼ�����BFS���ٽ��µ���ŵ�BFS���
    %LonelyNode
    %LonelyNum
    Cmartix(LonelyNode,:) = [];
    Cmartix(:,LonelyNode) = [];
    subResult = BFS( Cmartix );
    [subNum,totalNum] = size(subResult);
    % �µ���ŵ�BFS���  b=[a(:,1),zeros(5,1),a(:,2:5)]
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

