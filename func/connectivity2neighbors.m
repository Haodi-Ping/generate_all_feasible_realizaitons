function [neighbors] = connectivity2neighbors(ConnectivityM)
%CONNECTIVITY2NEIGHBORS �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
nClean = size(ConnectivityM,1);
nbrNums = sum(ConnectivityM,2);
maxNbrNum = max(nbrNums);
neighbors = zeros(nClean, maxNbrNum);% neighbors ����Ϊ���ڵ���������Ϊ������ھ���
neighbors(:,1)=nbrNums;% neighbors ��һ�д���ǰ�ڵ��ھ���
for i=1:nClean
    cnt = 1;
    for j=1:nClean
        if ConnectivityM(i,j)==1
            cnt=cnt+1;
            neighbors(i,cnt)=j;
        end
    end
end

end

