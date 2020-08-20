function plot_iteration(POS,POS_itera,npoints,iteration_num)
% PP = POS{}
BoxScale = 100;
%figure
idx1 = 1; idx2 = 2;
markersize = 3;
%% plot lines
for i = 1:npoints
    hold on
    PP = POS_itera{i};
%    this_itera_start = PP(1:end-1,:);
%    this_itera_end = PP(2:end,:);
%     hh = plot([this_itera_start(:,idx1); this_itera_end(:,idx1)],...
%             [this_itera_start(:,idx2); this_itera_end(:,idx2)],...
%             'color',[.5,.5,.5]);
    for j = 1:iteration_num - 1
       hh =  plot([PP(j,idx1); PP(j+1,idx1)],...
             [PP(j,idx2); PP(j+1,idx2)],...
             'color',[.5,.5,.5]);
        set(hh,'linewidth',1.5);
    end
end
%% plot points
for i = 2:iteration_num-1
    hold on
    PP = POS{i};
    h = plot(PP(:,idx1),PP(:,idx2),'ob','markersize',markersize,'linewidth',3.5); % 画ground truth点
    set(h,'linewidth',1.5);
end
%% plot initial state
PP = POS{1};
hold on
h=plot(PP(:,idx1),PP(:,idx2),'d','markersize',3.8,'linewidth',2,'color',[48/255,128/255,20/255]); % 画ground truth点
%% plot final state
PP = POS{end};
hold on
h=plot(PP(:,idx1),PP(:,idx2),'+r','markersize',8,'linewidth',2); % 画ground truth点

%% set axis
axis('square'); 
size_font = 10;
set(gca,'YTick',-50:20:50,'XTick',-50:20:50);%,'Fontname', 'Times New Roman'
set(gca,'FontSize',8) %%设置横纵坐标字体的大小

xlabel('X','fontsize',size_font,'color','k');
ylabel('Y','fontsize',size_font,'color','k');
grid on;
box on;
% 
% axis(0.6*BoxScale*[-1,1,-1,1]);
% set(gca,'ytick',(-60:20:60))
% set(gca,'xtick',(-60:20:60))

%xlabel('X','FontSize',20),ylabel('Y','FontSize',20);
% xlabel('X'),ylabel('Y');
% grid on;
% box on;
%title([graphlabel ': ' num2str(mean_residual)]);
pause(0.1); 
hold off

end

