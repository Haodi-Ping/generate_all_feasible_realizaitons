%%************************************************************************
%% Plots anchor points, estimations of unknown points and actual locations
%% and the discrepancies between them
%% Blue Diamond : Anchor
%% Red Star : Estimated position of unknown point
%% Green Circle : Actual position of unknown point
%%
%% Input: 
%% P0: Anchor positions. If there is no anchor, input P0 = [];
%% PP: Actual positions of unknown points
%% X_opt: Estimated positions of unknown points
%% Plane 'xy','yz','xz' : Desired 2-D plane if points are 3-D
%%       'xyz'          : 3-D 
%%************************************************************************
% plotpositionsaa(points',pos_g2o',ConnectivityM,Boxscale,['G2O: ' num2str(mean_residual_g2o)]);
   function plotpositionsaa(PP,Xopt,ConnectivityM,BoxScale,graphlabel)

   [~,mean_residual] = F_residual(Xopt',PP',BoxScale);
   h = figure;
   set(h,'name',graphlabel,'Numbertitle','off')
   %axes('FontSize',18,'FontWeight','bold');
   
   markersize = 9; 
   PP=PP*BoxScale;
   Xopt=Xopt*BoxScale;
   dim = size(Xopt,1);
   plane = 'xy';
   if (dim == 2); plane = 'xy'; end
   if (dim == 3); plane = 'xyz'; end
 
   if strcmp(plane,'xy')
      if strcmp(plane,'xy')
         idx1 = 1; idx2 = 2; 
      elseif strcmp(plane,'xz')
         idx1 = 1; idx2 = 3;
      elseif strcmp(plane,'yz')
         idx1 = 2; idx2 = 3;
      end
      h=plot(Xopt(idx1,:),Xopt(idx2,:),'+','markersize',markersize,'linewidth',2,'color',[.2 .7 .3]); % 画实现点 'color',[.5 .5 .5]
      hold on; %grid on
      h=plot([Xopt(idx1,:); PP(idx1,:)],[Xopt(idx2,:); PP(idx2,:)],'b');
      set(h,'linewidth',2);
      h=plot(PP(idx1,:),PP(idx2,:),'or','markersize',markersize,'linewidth',2); % 画ground truth点
      set(h,'linewidth',1.5);
      axis('square'); axis(0.6*BoxScale*[-1,1,-1,1]);
      %xlabel('X','FontSize',20),ylabel('Y','FontSize',20);
      xlabel('X'),ylabel('Y');
      grid on;
      title([graphlabel ': ' num2str(mean_residual)]);
      pause(0.1); hold off
   end
   %% 画原图
   hold on;
  %[dim,nfix] = size(P0);  
  [dim,npts] = size(PP); 
    DD = ConnectivityM;
  Ds = DD(1:npts,1:npts); 
  %Da = DD(1:npts,npts+(1:nfix)); 

  r1 = 1; r2 = 2; 

  hold on;      
  for j = 2:npts           
     idx = find(Ds(1:j,j)); 
     if ~isempty(idx)
        len = length(idx); 
        Pj = PP(:,j)*ones(1,len); 
        hh = plot([Pj(r1,:); PP(r1,idx)],[Pj(r2,:); PP(r2,idx)],'color',[.5,.5,.5],'linewidth',1.1);
        %set(hh,'linewidth',0.5);
     end    
  end
  h=plot(PP(idx1,:),PP(idx2,:),'or','markersize',markersize,'linewidth',2); % 画ground truth点
  for i=1:npts
    %h = plot(PP(r1,i),PP(r2,i),'or','markersize',6);  
    text(PP(r1,i)+1.5,PP(r2,i)+1.5,num2str(i));
    set(h,'linewidth',2);
    hold on
  end
  

  pause(0.2)
%   xlim([-60 60]);
%   ylim([-60 60]);
  
  end
%%************************************************************************


