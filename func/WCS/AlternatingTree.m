function [re, level, G]=AlternatingTree(B,M,root)
% WCS ��
b=B.Edges.EndNodes;
M=[M; M(:,2), M(:,1)];
m=max(b(:));                %ѹ���������ֵ�����ڽӾ���Ŀ����
A=compresstable2matrix(b);  %���ڽ�ѹ������ͼ�ľ����ʾ
%G=graph(b(:,1),b(:,2))
%figure
%plot(G);
%type=[0,1];

head=1;             %����ͷ
tail=1;             %����β����ʼ����Ϊ�գ�tail==head
queue(head)=root;      %��ͷ�м���ͼ��һ���ڵ�
level(head)=0;
head=head+1;        %������չ

flag=root;             %���ĳ���ڵ��Ƿ���ʹ���
re=[];              %���ս��
while tail~=head    %�ж϶����Ƿ�Ϊ��
    i=queue(tail);  %ȡ��β�ڵ�
    for j=1:m
        if A(i,j)==1 &&isempty(find(flag==j,1) )    %����ڵ���������û�з��ʹ�
            if(mod(level(tail),2)==0)               %�����0,2,4,6,8��
                if(ismember([i,j],M,'rows')==0)     %Ҫ������δƥ��ı�
                    queue(head)=j;                          %�½ڵ�����
                    level(head)=level(tail)+1;
                    head=head+1;                            %��չ����
                    flag=[flag j];                          %���½ڵ���б��
                    re=[re;i j];                            %���ߴ�����
                end
            elseif(mod(level(tail),2)==1)           %�����1,3,5,7��
                if(ismember([i,j],M,'rows')==1)     %Ҫ������ƥ��ı�
                    queue(head)=j;                          %�½ڵ�����
                    level(head)=level(tail)+1;
                    head=head+1;                            %��չ����
                    flag=[flag j];                          %���½ڵ���б��
                    re=[re;i j];                            %���ߴ�����
                end 
            end
        end
    end
    tail=tail+1;            
end

A=compresstable2matrix(re);
G=graph(A);


end