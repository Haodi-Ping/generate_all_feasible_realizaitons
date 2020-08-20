function [separator,QI_connt,subgraph,cut_set,ConnectivityM,distMatrix] = ...
    Func_Partition(graph,radius,points,ConnectivityM,distMatrix)
%Func_Partition graph 1��n�У�Ϊ1��Ԫ�ش����ͼ�����˽ڵ�
%   points ֻ�� MI check ��Ҫ
%% Ӧ�ý� QI �� MI �ӵıߴ��ϣ�����һ���µ� ConnectivityM �� distMatrix ����
%% ����������
separator = 0; % ��¼separator�ĸ���
QI_connt = 0;
%MI_count = 0;
subgraph = []; % ��ŷָ���
cut_set = []; % ��Ŵ˷ָ����� cut set
% disp('Func_Partition')
if(sum(graph) == 0)
    % ��ͼ do nothing
    % disp('------ Func_Partition: ��ͼ do nothing ------------------------------');
else
    %% Ԥ�����ҵ��ؽڽڵ㣬�ǲ���������cut set��
    % disp('------ Func_Partition: Ԥ�����ؽڽڵ㲻�� cut set ------------------------------')
    arthro = []; % ��¼�ҵ��Ĺؽڽڵ�
    arthro_num = 0;
    % [~,num] = size(ConnectivityM);
    ID = find (graph>0);
    [~,nodes_num] = size(ID);
    for i = 1:nodes_num
        tempID = graph; % ������Ա��ݵ� graph
        tempID(ID(i))=0; % �� graph ��ɾ����ǰ��
        subtempID = find(tempID>0); % subtempID ��ȥ����ǰ "�ؽڵ�" ��ʣ�µĵ�
        subtempM = ConnectivityM(subtempID,subtempID); %��ԭͼ��ȡ�� ȥ����ǰ cut set �� ����ͼ��ϵ --- ��ͼ�ڵ��� �� ��ͼ�ڵ���
        tempresult=Func_SuperBFS(subtempM);  %���� ȥ����ǰ cut set �� ��ͨ��ͼ
        [tempa,~]=size(tempresult);% a Ϊ��������˴������ԣ�ȫ��Ϊ 2��1��separator������ >2 ��ģ������
        if (tempa>1)
            arthro_num = arthro_num+1;
            arthro(arthro_num) = ID(i);
        end
    end
    %arthro
    %% ��ʼ�ָ�
    %disp('------ Func_Partition: ��ʼ�ָ� ------------------------------');
    state=0; % ����Ƿ�ɹ�
    [~,num] = size(ConnectivityM);
    ID = find (graph>0);
    [~,nodes_num] = size(ID);
    for i=1:nodes_num
        % ��� i �ǹؽڽڵ㣬����
        is_i_arthro = size(find(arthro == ID(i)));
        if (is_i_arthro>0)
            continue;
        end
        
        % ��� j �ǹؽڽڵ㣬����
        for j=i+1:nodes_num
            is_j_arthro = size(find(arthro == ID(j)));
            if (is_j_arthro>0)
                continue;
            end
            % ��ʼ�ָ�
            
            markID=graph;% ����markID ��ʼ��Ϊ�����ͼ 1��npoints�� --- 1 �� 100
            markID(ID(i))=0; %
            markID(ID(j))=0; % ����ɾ��ID(i)��ID��j�����ڵ㣬������ cut set = {v_ID(i),v_ID(j) }
            if ConnectivityM(ID(i),ID(j))==0
                continue
            end
            subID=find(markID>0); % subID ��ȥ����ǰ cut set ��ʣ�µĵ�
            subM=ConnectivityM(subID,subID); %��ԭͼ��ȡ�� ȥ����ǰ cut set �� ����ͼ��ϵ --- ��ͼ�ڵ��� �� ��ͼ�ڵ���
            result=Func_SuperBFS(subM);  %���� ȥ����ǰ cut set �� ��ͨ��ͼ
            [a,~]=size(result); % a Ϊ��������˴������ԣ�ȫ��Ϊ 2��1��separator������ >2 ��ģ������
            
            if (a>1)
                %Pending_Cut_Set=[ID(i) ID(j)]
                % �ҵ�һ�� separator
                % ���QI��������ˣ��������ս����
                [ConnectivityM,distMatrix,QI_flag] = ...
                    Func_Check_QI(ID(i),ID(j),ConnectivityM,distMatrix,radius,points);
                %flag_arthro = Func_Find_Arthro(result,ConnectivityM);
                %if (QI_flag == 0 && flag_arthro == 0)
                separator = separator+1;
                if (QI_flag == 0 )
                    % QI ʧ��,���ֽ���
                    % MI_flag = Func_Check_MI(cut_set,subgraph,radius,points,distMatrix,ConnectivityM)
                    subgraph=zeros(a,num); % a Ϊ�г�������ͼ������numΪ�ڵ���
                    cut_set = [ID(i),ID(j)];
                    for k=1:a % �����г�����ÿ����ͼ,�浽 subgraph ��
                        resultID=find(result(k,:)>0); % ����ǰ��ͼ�����Ľڵ�浽 resultID
                        subgraph(k,subID(resultID))=1; % ����ǰ��ͼ�����Ľڵ�浽 todoM
                        subgraph(k,ID(i))=1;
                        subgraph(k,ID(j))=1; % ��  cut set �浽 todoM
                        %subG = find(subgraph(k,:)>0)
                    end
%                     for ii = 1:a
%                         sub_ii = find(subgraph(ii,:)>0)
%                     end
                    state=1;
                    break;
                else
                    % �� separator �� QI ���
                    % ���� dismat ���мӱ�
                    QI_connt = QI_connt+1;
                    continue;
                end
                %
                
                
                %             num1 = sum(subgraph(1,:));
                %             num2 = sum(subgraph(2,:));
                %             if (num1<=3 || num2<=3) % ��3�����ͼ��ֱ�Ӵ�
                %                 state=1;
                %                 break;
                %             end
                
                
                
                %
                
                
            end
        end
        if(state==1)% ֻҪ�ɹ��ҵ�һ��  cut set ��ֹͣ�� i
            break;
        end
    end
    
end

end

