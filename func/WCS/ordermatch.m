function  [matching,subgraph]=ordermatch(G, A, n, m, points)
% G �Ǹ��ݶ� ConnectivityM �õ���G��ʽ��Edges��Nodes��ͼ
% A �� ConnectivityM ��������û�õ�
% n �Ǳ���
% m �ǽڵ���
% points Ҳû��
num=0;
subgraph(1,:)=zeros(1,n);
Redge=[];
matching(:,1)=1:n+3;
matching(:,2)=0;
B=zeros(n+3+2*m);
C=B+inf;
for i=1:n
    for k=1:2
        B(i, n+3+G.Edges.EndNodes(i,k)*2-1)=1;
        B(n+3+G.Edges.EndNodes(i,k)*2-1,i)=1;
        B(i, n+3+G.Edges.EndNodes(i,k)*2)=1;
        B(n+3+G.Edges.EndNodes(i,k)*2,i)=1;
        
        C(i, n+3+G.Edges.EndNodes(i,k)*2-1)=1;
        C(n+3+G.Edges.EndNodes(i,k)*2-1,i)=1;
        C(i, n+3+G.Edges.EndNodes(i,k)*2)=1;
        C(n+3+G.Edges.EndNodes(i,k)*2,i)=1;
    end
    B(n+1:n+3,:)=[1,1,1]'*B(i,:);
    C(n+1:n+3,:)=[1,1,1]'*C(i,:);
    %%%%%%%%%%%%%�ڶ���ͼ�ϼ������ƥ��%%%%%%%%%%%%
    [assign,cost]=assignmentoptimal(C);
    %%%%%%%%%%%%% matching ������ƥ��%%%%%%%%%%%%%
    submatch=[(1:n+3)' assign(1:n+3)];
    submatch(Redge,2)=1;
    test1 = find(submatch(1:i,2)==0);
    submatch(Redge,2)=0;
    test2 = find(submatch(n+1:n+3,2)==0);
    if(~isempty(test1)||~isempty(test2))
        BG=graph(B+B');
        if(isempty(test1))
            opt=i;
        else
            opt=test1;
        end;
        [edge, level, AltTree]=AlternatingTree(BG, submatch, opt);
        ID = find(edge(:,2)<n+1);
        result=zeros(1,n);
        resultID = edge(ID,2)';
        result(resultID)=1;
        result(opt)=1;
        num=num+1;
        Redge(num)=i;
        subgraph(num,:)=result;
        
        connectSub = zeros(m);
        for j=find(result)
            connectSub(G.Edges.EndNodes(j,1),G.Edges.EndNodes(j,2))=1;
        end
        %plotgraph([],points,connectSub);
        
        B(i,:)=0;
        C(i,:)=inf;
    end
    %i
end
matching = submatch;
end

