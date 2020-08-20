%%**************************************************************
%% Plot connectivity graph
%%
%%**************************************************************

  function plotgraph(P0,PP,DD,graphlabel)
  h = figure;
  set(h,'name',graphlabel,'Numbertitle','off')

  %axes('FontSize',14,'FontWeight','bold');
  markersize = 9; 


  [dim,nfix] = size(P0);  
  [dim,npts] = size(PP); 

  Ds = DD(1:npts,1:npts); 
  Da = DD(1:npts,npts+(1:nfix)); 

  r1 = 1; r2 = 2; 

  hold on;      
  for j = 2:npts           
     idx = find(Ds(1:j,j)); 
     if ~isempty(idx)
        len = length(idx); 
        Pj = PP(:,j)*ones(1,len); 
        plot([Pj(r1,:); PP(r1,idx)],[Pj(r2,:); PP(r2,idx)],'k'); % »­Ïß
     end    
  end
  if ~isempty(P0)
     for j = 1:nfix         
        idx = find(Da(:,j)); 
        len = length(idx); 
        Pj = P0(:,j)*ones(1,len); 
        h = plot([Pj(r1,:); PP(r1,idx)],[Pj(r2,:); PP(r2,idx)],'b');
        set(h,'linewidth',2);    
        h = plot(P0(r1,j),P0(r2,j),'bd','markersize',markersize);   
        set(h,'linewidth',3);      
     end
  end
  for i=1:npts
    h = plot(PP(r1,i),PP(r2,i),'or','markersize',markersize);    
    set(h,'linewidth',2);
    hold on
  end
  
  for i=1:npts
    h = plot(PP(r1,i),PP(r2,i),'or','markersize',markersize);  
    text(PP(r1,i)+1.5,PP(r2,i)+1.5,num2str(i));
    set(h,'linewidth',2);
    hold on
  end
  
  xlim([-60 60]);
  ylim([-60 60]);
  grid on
  box on
  axis('square');
  hold off
%%**************************************************************
