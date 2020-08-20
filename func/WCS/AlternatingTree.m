function [re, level, G]=AlternatingTree(B,M,root)
% WCS 用
b=B.Edges.EndNodes;
M=[M; M(:,2), M(:,1)];
m=max(b(:));                %压缩表中最大值就是邻接矩阵的宽与高
A=compresstable2matrix(b);  %从邻接压缩表构造图的矩阵表示
%G=graph(b(:,1),b(:,2))
%figure
%plot(G);
%type=[0,1];

head=1;             %队列头
tail=1;             %队列尾，开始队列为空，tail==head
queue(head)=root;      %向头中加入图第一个节点
level(head)=0;
head=head+1;        %队列扩展

flag=root;             %标记某个节点是否访问过了
re=[];              %最终结果
while tail~=head    %判断队列是否为空
    i=queue(tail);  %取队尾节点
    for j=1:m
        if A(i,j)==1 &&isempty(find(flag==j,1) )    %如果节点相连并且没有访问过
            if(mod(level(tail),2)==0)               %如果是0,2,4,6,8层
                if(ismember([i,j],M,'rows')==0)     %要求搜索未匹配的边
                    queue(head)=j;                          %新节点入列
                    level(head)=level(tail)+1;
                    head=head+1;                            %扩展队列
                    flag=[flag j];                          %对新节点进行标记
                    re=[re;i j];                            %将边存入结果
                end
            elseif(mod(level(tail),2)==1)           %如果是1,3,5,7层
                if(ismember([i,j],M,'rows')==1)     %要求搜索匹配的边
                    queue(head)=j;                          %新节点入列
                    level(head)=level(tail)+1;
                    head=head+1;                            %扩展队列
                    flag=[flag j];                          %对新节点进行标记
                    re=[re;i j];                            %将边存入结果
                end 
            end
        end
    end
    tail=tail+1;            
end

A=compresstable2matrix(re);
G=graph(A);


end