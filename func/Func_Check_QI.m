function [connectivityMat,distMat,QI_flag] = Func_Check_QI(cut1, cut2,connectivityMat,distMat,radius,points)
%Func_Check_QI 
%   传入QI仅为在控制台输出ground truth，检查QI计算结果
%disp('Func_Check_QI')
QI_flag = 0; % 为 0 代表此 separator 不能被 QI 解决
if(connectivityMat(cut1,cut2)==0)
    % cut set无边 --- do nothing, QI 只能解决 cut set 有边的情形
    %disp('QI 未加边：cut set间无连接  [cut1 cut2]')
    %[cut1 cut2]
else
    % cut set有边 --- 构造 r-Quad 并检查 QI
    [~,pntNum] = size(connectivityMat); % pntNum 为节点数
    comNei = []; % cut1, cut2 的公共邻居
    count_comNei = 0; % cut1, cut2 的公共邻居 数量
    % 构造公共邻居
    for i = 1:pntNum
        if(connectivityMat(i,cut1)==1 && connectivityMat(i,cut2)==1 && i~=cut1 && i~=cut2)
            count_comNei = count_comNei+1;
            comNei(count_comNei) = i;
        end
    end
    if(count_comNei<2)
        % 不够2个公共邻居，无法构造 r-Quad --- do nothing
        %disp('QI 未加边：cut set无多余2个公共邻居  [cut1 cut2]')
        % [cut1 cut2]
    else
        % 构造 r-Quad 并检查QI
        for i = 1:count_comNei
            for j = i+1:count_comNei
                if(connectivityMat(comNei(i),comNei(j))==0) % 公共邻居没有边
                    e = distMat(cut1,cut2);
                    a = distMat(cut1,comNei(i));
                    b = distMat(cut2,comNei(i));
                    c = distMat(cut1,comNei(j));
                    d = distMat(cut2,comNei(j));
                    ang_cal_1 = (a*a+e*e-b*b)/(2*a*e);
                    ang_cal_2 = (c*c+e*e-d*d)/(2*c*e);
                    ang_1 = acos(ang_cal_1);
                    ang_2 = acos(ang_cal_2);
                    x_1 = sqrt(a*a + c*c - 2*a*c*cos(ang_1+ang_2));
                    x_2 = sqrt(a*a + c*c - 2*a*c*cos(abs(ang_1-ang_2)));
                    if (x_1 >radius && x_2 <= radius)
                        %disp("---QI check 成功:"+cut1+","+cut2);
                        QI_flag = 1;
                        connectivityMat(comNei(i),comNei(j))=1;
                        connectivityMat(comNei(j),comNei(i))=1;
                        distMat(comNei(i),comNei(j))=x_1;
                        distMat(comNei(j),comNei(i))=x_1;
                        %disp('QI 加边：cut1 cut2 comNei(i) comNei(j) QI_Length GroundTruth')
                        %gt = (points(comNei(i),1)-points(comNei(j),1))^2 + (points(comNei(i),2)-points(comNei(j),2))^2;
                        %gt = sqrt(gt);
                        %[cut1 cut2 comNei(i) comNei(j) x_1 gt]
                        %if(abs(distMat(comNei(i),comNei(j))-x_1)>0.05)
                            %disp("---difference found:");
                        %end
                        break;
                    else
                        %flag = 0;  % do nothing
                        %disp("---QI check 不成立:"+cut1+","+cut2);
                        %disp('QI 未加边：基于公共邻居QI条件不成立  [cut1 cut2 comNei(i) comNei(j)]')
                        % [cut1 cut2 comNei(i) comNei(j)]
                    end
                end
            end
            
            if (QI_flag == 1)
                break;
            end
            
        end
        
    end
end

